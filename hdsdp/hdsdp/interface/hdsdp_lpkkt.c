#ifdef HEADERPATH
#include "interface/hdsdp_lpkkt.h"
#else
#include "hdsdp_lpkkt.h"
#endif

/* Given scaling matrix D. Solve linear system
 
    A D^2 A' y = r
 
 using
 
 1. Cholesky factorization
 2. Preconditioned conjugate residual or conjugate gradient
 
*/

static hdsdp_retcode HLpKKTIAnalyze( hdsdp_lp_kkt *kkt, int iNormalEqn ) {
    
    /* Analyze the sparsity pattern of normal equation */
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    assert( iNormalEqn );
    int *iNzBuffer = NULL;
    
    if ( iNormalEqn ) {
        
        /* Normal equation A * D^2 * A' */
        int nKKTCol = kkt->nLpRow;
        int nKKTNz = 0;
        int nColNz = 0;
        
        int nKKTAlloc = nKKTCol * 100;
        
        HDSDP_INIT(kkt->KKTMatBeg, int, nKKTCol + 1);
        HDSDP_INIT(kkt->KKTMatIdx, int, nKKTAlloc);
        HDSDP_INIT(iNzBuffer, int, nKKTCol);
        
        HDSDP_MEMCHECK(kkt->KKTMatBeg);
        HDSDP_MEMCHECK(iNzBuffer);
        
        for ( int iKKTCol = 0; iKKTCol < nKKTCol; ++iKKTCol ) {
            
            nColNz = 0;
            HDSDP_ZERO(iNzBuffer, int, nKKTCol);
            for ( int iTransElem = kkt->colMatTransBeg[iKKTCol]; iTransElem < kkt->colMatTransBeg[iKKTCol + 1]; ++iTransElem ) {
                int iTransRow = kkt->colMatTransIdx[iTransElem];
                for ( int iLpElem = kkt->colMatBeg[iTransRow]; iLpElem < kkt->colMatBeg[iTransRow + 1]; ++iLpElem ) {
                    iNzBuffer[kkt->colMatIdx[iLpElem]] = 1;
                }
            }
            
            /* Count number of nonzeros in each KKT column 
               Only consider the lower triangular part of the matrix */
            for ( int iRow = iKKTCol; iRow < nKKTCol; ++iRow ) {
                if ( iNzBuffer[iRow] ) {
                    kkt->KKTMatIdx[nKKTNz] = iRow;
                    nColNz += iNzBuffer[iRow];
                    nKKTNz += 1;
                }
            }
            
            /* Allocate more memory in case it runs out */
            kkt->KKTMatBeg[iKKTCol + 1] = nKKTNz;
            
            if ( nKKTNz >= nKKTAlloc - nKKTCol - 1 ) {
                nKKTAlloc = nKKTAlloc * 2;
                int *KKTMatIdxTmp = kkt->KKTMatIdx;
                HDSDP_REALLOC(kkt->KKTMatIdx, int, nKKTAlloc);
                
                if ( !kkt->KKTMatIdx ) {
                    kkt->KKTMatIdx = KKTMatIdxTmp;
                    retcode = HDSDP_RETCODE_FAILED;
                    goto exit_cleanup;
                }
            }
        }
        
        /* Done with analysis */
        HDSDP_REALLOC(kkt->KKTMatIdx, int, nKKTNz);
        
        /* Allocate memory for matrix elements */
        HDSDP_INIT(kkt->KKTMatElem, double, nKKTNz);
        HDSDP_MEMCHECK(kkt->KKTMatElem);
        
        /* Symbolic analysis */
        HDSDP_CALL(HFpLinsysCreate(&kkt->KKTDirect, nKKTCol, HDSDP_LINSYS_SPARSE_DIRECT));
        HDSDP_CALL(HFpLinsysSymbolic(kkt->KKTDirect, kkt->KKTMatBeg, kkt->KKTMatIdx));
        
    } else {
        
        assert( 0 );
        /* Augmented system
            [ D      A'  ]
            [ A   -r * I ] */
        
        int nRow = kkt->nLpRow;
        int nCol = kkt->nLpCol;
        int nElem = kkt->nLpElem;
        
        const int *colMatBeg = kkt->colMatBeg;
        const int *colMatIdx = kkt->colMatIdx;
        const double *colMatElem = kkt->colMatElem;
        
        kkt->nKKTAugCol = nRow + nCol;
        kkt->nKKTAugElem = nRow + nCol + kkt->nLpElem;
        
        HDSDP_INIT(kkt->KKTMatBeg, int, kkt->nKKTAugCol + 1);
        HDSDP_INIT(kkt->KKTMatIdx, int, kkt->nKKTAugElem);
        HDSDP_INIT(kkt->KKTMatElem, double, kkt->nKKTAugElem);
        
        /* Set up symbolic structure of KKT system */
        /* First n columns */
        for ( int iCol = 0; iCol < nCol; ++iCol ) {
            kkt->KKTMatBeg[iCol] = colMatBeg[iCol] + iCol;
            kkt->KKTMatIdx[kkt->KKTMatBeg[iCol]] = iCol;
            for ( int iElem = colMatBeg[iCol]; iElem < colMatBeg[iCol + 1]; ++iElem ) {
                kkt->KKTMatBeg[iCol + iElem + 1] = colMatIdx[iElem] + nCol;
                kkt->KKTMatElem[iCol + iElem + 1] = colMatElem[iElem];
            }
        }
        
        kkt->KKTMatBeg[nCol] = kkt->nLpElem + nCol;
        
        /* Last m columns */
        for ( int iCol = 0; iCol < nRow; ++iCol ) {
            kkt->KKTMatBeg[iCol + nCol + 1] = kkt->KKTMatBeg[iCol + nCol] + 1;
            kkt->KKTMatIdx[iCol + nCol + 1] = nCol + iCol;
        }
        
        HDSDP_CALL(HFpLinsysCreate(&kkt->KKTDirect, nRow + nCol, HDSDP_LINSYS_SPARSE_INDEFINITE));
        HDSDP_CALL(HFpLinsysSymbolic(kkt->KKTDirect, kkt->KKTMatBeg, kkt->KKTMatIdx));
    }
    
exit_cleanup:
    
    HDSDP_FREE(iNzBuffer);
    
    return retcode;
}

static hdsdp_retcode HLpKKTISetupAugmented( hdsdp_lp_kkt *kkt, double *dScalingMatrix, double dPrimalReg, double dDualReg ) {

    /* Method should not be invoked unless an augmented solver is ready */
    assert(0);
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    for ( int iCol = 0; iCol < kkt->nLpCol; ++iCol ) {
        kkt->KKTMatElem[kkt->KKTMatBeg[iCol]] = 1.0;
        for ( int iElem = kkt->colMatBeg[iCol]; iElem < kkt->colMatBeg[iCol + 1]; ++iElem ) {
            kkt->KKTMatElem[iCol + iElem + 1] = kkt->colMatElem[iElem] * dScalingMatrix[iCol];
        }
    }
    
    HDSDP_CALL(HFpLinsysNumeric(kkt->KKTDirect, kkt->KKTMatBeg, kkt->KKTMatIdx, kkt->KKTMatElem));
 
exit_cleanup:
    return retcode;
}

static hdsdp_retcode HLpKKTISetupNormal( hdsdp_lp_kkt *kkt, double *dScalingMatrix, double dPrimalReg, double dDualReg ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    double dTStart = HUtilGetTimeStamp();
    
    /* Regularize scaling matrix */
    for ( int iCol = 0; iCol < kkt->nLpCol; ++iCol ) {
        kkt->dScalingBuffer[iCol] = dScalingMatrix[iCol] * dScalingMatrix[iCol] + dPrimalReg;
    }
    
    /* Clear the normal matrix */
    HDSDP_ZERO(kkt->KKTMatElem, double, kkt->KKTMatBeg[kkt->nLpRow]);
    /* Add dual regularization to the diagonal of the normal matrix
       Since we only store the lower triangular part of the normal matrix, the diagonal must be on the first
       row of each column.
     */
    for ( int iKKTCol = 0; iKKTCol < kkt->nLpRow; ++iKKTCol ) {
        kkt->KKTMatElem[kkt->KKTMatBeg[iKKTCol]] = dDualReg;
    }
    
    double dNormalEqnElem = 0.0;
    for ( int iKKTCol = 0; iKKTCol < kkt->nLpRow; ++iKKTCol ) {
        HDSDP_ZERO(kkt->dScaledColBuffer, double, kkt->nLpCol);
        /* Scale the row iKKTCol by scaling matrix */
        for ( int iColElem = kkt->colMatTransBeg[iKKTCol]; iColElem < kkt->colMatTransBeg[iKKTCol + 1]; ++iColElem ) {
            int iCol = kkt->colMatTransIdx[iColElem];
            kkt->dScaledColBuffer[iCol] = kkt->dScalingBuffer[iCol] * kkt->colMatTransElem[iColElem];
        }
        
        for ( int iRowElem = kkt->KKTMatBeg[iKKTCol]; iRowElem < kkt->KKTMatBeg[iKKTCol + 1]; ++iRowElem ) {
            dNormalEqnElem = 0.0;
            int iRow = kkt->KKTMatIdx[iRowElem];
            /* Compute inner product between row [iRow] and [iKKTCol] of A */
            for ( int iElem = kkt->colMatTransBeg[iRow]; iElem < kkt->colMatTransBeg[iRow + 1]; ++iElem ) {
                dNormalEqnElem += kkt->dScaledColBuffer[kkt->colMatTransIdx[iElem]] * kkt->colMatTransElem[iElem];
            }
            kkt->KKTMatElem[iRowElem] += dNormalEqnElem;
        }
    }
    
    /* Set up the normal equation with scaling matrix */
    HDSDP_CALL(HFpLinsysNumeric(kkt->KKTDirect, kkt->KKTMatBeg, kkt->KKTMatIdx, kkt->KKTMatElem));
    kkt->dFactorTime += HUtilGetTimeStamp() - dTStart;
 
exit_cleanup:
    return retcode;
}

extern hdsdp_retcode HLpKKTCreate( hdsdp_lp_kkt **pkkt ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    if ( !pkkt ) {
        retcode = HDSDP_RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    hdsdp_lp_kkt *kkt = NULL;
    HDSDP_INIT(kkt, hdsdp_lp_kkt, 1);
    HDSDP_MEMCHECK(kkt);
    HDSDP_ZERO(kkt, hdsdp_lp_kkt, 1);
    
    *pkkt = kkt;
    
exit_cleanup:
    return retcode;
}

extern hdsdp_retcode HLpKKTInit( hdsdp_lp_kkt *kkt, int nRow, int nCol, int *colMatBeg, int *colMatIdx, double *colMatElem, 
                                 int *colMatTransBeg, int *colMatTransIdx, double *colMatTransElem ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    kkt->nLpRow = nRow;
    kkt->nLpCol = nCol;
    kkt->nLpElem = colMatBeg[nCol];
    
    kkt->colMatBeg = colMatBeg;
    kkt->colMatIdx = colMatIdx;
    kkt->colMatElem = colMatElem;
    kkt->colMatTransBeg = colMatTransBeg;
    kkt->colMatTransIdx = colMatTransIdx;
    kkt->colMatTransElem = colMatTransElem;
    
    HDSDP_INIT(kkt->dScalingBuffer, double, nCol);
    HDSDP_MEMCHECK(kkt->dScalingBuffer);
    
    HDSDP_INIT(kkt->dScaledColBuffer, double, nCol);
    HDSDP_MEMCHECK(kkt->dScaledColBuffer);
    
    /* Start the analysis of normal equation */
    HDSDP_CALL(HLpKKTIAnalyze(kkt, 1));
    
    /* Create sparse iterative solver */
    HDSDP_CALL(HFpLinsysCreate(&kkt->KKTIterative, nRow, HDSDP_LINSYS_SPARSE_ITERATIVE));
    HDSDP_CALL(HFpLinsysSymbolic(kkt->KKTIterative, kkt->KKTMatBeg, kkt->KKTMatIdx));
    
    kkt->nFactor = 0;
    kkt->nSolve = 0;
    
exit_cleanup:
    return retcode;
}

extern hdsdp_retcode HLpKKTSetup( hdsdp_lp_kkt *kkt, int LpMethod, double *dScalingMatrix, double dPrimalReg, double dDualReg ) {
    
    /* Set up the KKT matrix for different types of iterations
       If LpMethod is LP_ITER_MEHROTRA, then the scaling matrix is X^{-0.5} S^{0.5}
       The augmented system is directly formed, factorized and solved by kkt->KKTDirect
        
       [ X^{-1} S    A' ] [ dx ] = [ r1 ]
       [     A       0  ] [ dy ] = [ r2 ]
     
       If LpMethod is LP_ITER_PRIMAL, then the scaling matrix is X^{-1}
       The augmented system is formed to be
       
            [ I    XA']
            [ AX    0 ]
     
       and iterative solver is solved normal equation with preconditioned iterative solver
     
          AX^2A' * dy = r1
     */
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    kkt->LpMethod = LpMethod;
    
    HDSDP_CALL(HLpKKTISetupNormal(kkt, dScalingMatrix, dPrimalReg, dDualReg));
    kkt->nFactor += 1;
    
exit_cleanup:
    return retcode;
}

extern hdsdp_retcode HLpKKTSolveAugmented( hdsdp_lp_kkt *kkt, double *dLhsVec, double *dRhsVec ) {

    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    assert( 0 );
    
    /* Now only primal-dual method can invoke augmented system */
    if ( kkt->LpMethod != LP_ITER_MEHROTRA || kkt->LpMethod != LP_ITER_HSD ) {
        retcode = HDSDP_RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    HDSDP_CALL(HFpLinsysSolve(kkt->KKTDirect, 1, dLhsVec, dRhsVec));
    kkt->nSolve += 1;
    
exit_cleanup:
    return retcode;
}

extern hdsdp_retcode HLpKKTSolveNormalEqn( hdsdp_lp_kkt *kkt, int nRhs, double *dLhsVec, double *dRhsVec ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    /* Now only primal method can invoke normal equation */
    if ( kkt->LpMethod == LP_ITER_PRIMAL ) {
        HDSDP_CALL(HFpLinsysSolve(kkt->KKTIterative, nRhs, dLhsVec, dRhsVec));
    } else {
        double dTStart = HUtilGetTimeStamp();
        HDSDP_CALL(HFpLinsysSolve(kkt->KKTDirect, nRhs, dLhsVec, dRhsVec));
        kkt->dSolveTime += HUtilGetTimeStamp() - dTStart;
        kkt->nSolve += 1;
    }
    
exit_cleanup:
    return retcode;
}

extern double HLpKKTGetFactorSolveTimeRatio( hdsdp_lp_kkt *kkt ) {
    
    if ( kkt->nFactor == 0 ) {
        return 1.0;
    }
    
    double dFactorTime = kkt->dFactorTime / kkt->nFactor;
    double dSolveTime =  kkt->dSolveTime / kkt->nSolve;
    
    return dFactorTime / dSolveTime;
}

extern void HLpKKTClear( hdsdp_lp_kkt *kkt ) {
    
    if ( !kkt ) {
        return;
    }
    
    HDSDP_FREE(kkt->KKTMatBeg);
    HDSDP_FREE(kkt->KKTMatIdx);
    HDSDP_FREE(kkt->KKTMatElem);
    
    HDSDP_FREE(kkt->dScalingBuffer);
    HDSDP_FREE(kkt->dScaledColBuffer);
    
    HFpLinsysDestroy(&kkt->KKTDirect);
    HFpLinsysDestroy(&kkt->KKTIterative);
    
    HDSDP_ZERO(kkt, hdsdp_lp_kkt, 1);
    
    return;
}

extern void HLpKKTDestroy( hdsdp_lp_kkt **pkkt ) {
    
    if ( !pkkt ) {
        return;
    }
    
    HLpKKTClear(*pkkt);
    HDSDP_FREE(*pkkt);
    
    return;
}
