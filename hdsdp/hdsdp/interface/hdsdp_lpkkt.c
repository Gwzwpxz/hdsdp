#ifdef HEADERPATH
#include "interface/hdsdp_lpkkt.h"
#else
#include "hdsdp_lpkkt.h"
#endif

/* Given scaling matrix D. Solve linear system
 
    [ D^-2 + Rp     A'   ] [ dy ] = [ r1 ]
    [ A           0 - Rd ] [ dx ] = [ r2 ]
 using
 
 1. Augmented system direct solver
 2. Preconditioned normal equation iterative solver
*/

static hdsdp_retcode HLpKKTISetupPrimal( hdsdp_lp_kkt *kkt, double *dScalingMatrix ) {
       
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

static hdsdp_retcode HLpKKTISetupPrimalDual( hdsdp_lp_kkt *kkt, double *dScalingMatrix ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    for ( int iCol = 0; iCol < kkt->nLpCol; ++iCol ) {
        kkt->KKTMatElem[kkt->KKTMatBeg[iCol]] = dScalingMatrix[iCol] * dScalingMatrix[iCol];
    }
    
    HDSDP_CALL(HFpLinsysNumeric(kkt->KKTDirect, kkt->KKTMatBeg, kkt->KKTMatIdx, kkt->KKTMatElem));
 
exit_cleanup:
    return retcode;
}

static void HLpKKTIRegularize( hdsdp_lp_kkt *kkt, double dPrimalReg, double dDualReg ) {
    
    if ( dPrimalReg > 0.0 ) {
        for ( int iCol = 0; iCol < kkt->nLpCol; ++iCol ) {
            kkt->KKTMatElem[kkt->KKTMatBeg[iCol]] = dPrimalReg;
        }
    }
    
    if ( dDualReg > 0.0 ) {
        for ( int iCol = 0; iCol < kkt->nLpRow; ++iCol ) {
            kkt->KKTMatElem[kkt->KKTMatBeg[kkt->nKKTAugCol - iCol - 1]] = - dDualReg;
        }
    }
    
    return;
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
    
    /* Augmented system
        [ I      XA  ]
        [ AX  -r * I ]
     */
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
    HDSDP_CALL(HFpLinsysCreate(&kkt->KKTIterative, nRow, HDSDP_LINSYS_SPARSE_ITERATIVE));
    HDSDP_CALL(HFpLinsysSymbolic(kkt->KKTDirect, kkt->KKTMatBeg, kkt->KKTMatIdx));
    HDSDP_CALL(HFpLinsysSymbolic(kkt->KKTIterative, kkt->KKTMatBeg, kkt->KKTMatIdx));
    
    kkt->nFactor = 0;
    kkt->nSolve = 0;
    
exit_cleanup:
    return retcode;
}

extern hdsdp_retcode HLpKKTSetup( hdsdp_lp_kkt *kkt, lp_method LpMethod, double *dScalingMatrix, double dPrimalReg, double dDualReg ) {
    
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
    
    if ( LpMethod == LP_ITER_PRIMAL ) {
        HDSDP_CALL(HLpKKTISetupPrimal(kkt, dScalingMatrix));
    } else {
        HDSDP_CALL(HLpKKTISetupPrimalDual(kkt, dScalingMatrix));
        HLpKKTIRegularize(kkt, dPrimalReg, dDualReg);
    }
    
exit_cleanup:
    return retcode;
}

extern hdsdp_retcode HLpKKTSolveAugmented( hdsdp_lp_kkt *kkt, double *dLhsVec, double *dRhsVec ) {

    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    /* Now only primal-dual method can invoke augmented system */
    if ( kkt->LpMethod != LP_ITER_MEHROTRA || kkt->LpMethod != LP_ITER_HSD ) {
        retcode = HDSDP_RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    HDSDP_CALL(HFpLinsysSolve(kkt->KKTDirect, 1, dLhsVec, dRhsVec));
    
exit_cleanup:
    return retcode;
}

extern hdsdp_retcode HLpKKTSolveNormalEqn( hdsdp_lp_kkt *kkt, double *dLhsVec, double *dRhsVec ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    /* Now only primal method can invoke normal equation */
    if ( kkt->LpMethod != LP_ITER_PRIMAL ) {
        retcode = HDSDP_RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    HDSDP_CALL(HFpLinsysSolve(kkt->KKTIterative, 1, dLhsVec, dRhsVec));
    
exit_cleanup:
    return retcode;
}

extern void HLpKKTClear( hdsdp_lp_kkt *kkt ) {
    
    if ( !kkt ) {
        return;
    }
    
    HDSDP_FREE(kkt->KKTMatBeg);
    HDSDP_FREE(kkt->KKTMatIdx);
    HDSDP_FREE(kkt->KKTMatElem);
    
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
