/** @file hdsdp\_lpsolve.c
 *  @brief Specialized LP solver
 *
 */

#ifdef HEADERPATH
#include "interface/hdsdp_utils.h"
#include "interface/hdsdp_linsolver.h"
#include "interface/hdsdp_lpsolve.h"
#include "interface/hdsdp_lpkkt.h"
#include "linalg/vec_opts.h"
#include "linalg/sparse_opts.h"
#else
#include "hdsdp_utils.h"
#include "hdsdp_linsolver.h"
#include "hdsdp_lpsolve.h"
#include "hdsdp_lpkkt.h"
#include "vec_opts.h"
#include "sparse_opts.h"
#endif

#include <math.h>

#define PRIMAL_HISTORY  (3)


/* Implement primal statistics for primal IPM. The following statistics are collected to
 determine the convergence of primal IPM
 
 1. Primal iteration sequence over history window size PRIMAL_HISTORY
 2. Primal barrier parameter sequence over the whole history
 3. Condition number estimate
 4. Euclidean distance between consecutive iterations
 5. Thresholded scaled distance between consecutive iterations
 
 */
static hdsdp_retcode HPrimalStatsCreate( hdsdp_primal_stats **pStats ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    if ( !pStats ) {
        retcode = HDSDP_RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    hdsdp_primal_stats *stats = NULL;
    HDSDP_INIT(stats, hdsdp_primal_stats, 1);
    HDSDP_MEMCHECK(stats);
    HDSDP_ZERO(stats, hdsdp_primal_stats, 1);
    
    *pStats = stats;
 
exit_cleanup:
    return retcode;
}

static hdsdp_retcode HPrimalStatsInit( hdsdp_primal_stats *stats, int nCol, hdsdp_lpsolver_params *params ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    stats->params = params;
    stats->nCol = nCol;
    HDSDP_INIT(stats->dPrimalMuHistory, double, params->nMaxIter);
    HDSDP_INIT(stats->dPrimalIterHistory, double, PRIMAL_HISTORY * nCol);
    HDSDP_MEMCHECK(stats->dPrimalIterHistory);
    
exit_cleanup:
    return retcode;
}

static void HPrimalStatsUpdate( hdsdp_primal_stats *stats, int iIter, double *dColVal, double dMu ) {
    
    HDSDP_MEMCPY(stats->dPrimalIterHistory + stats->nCol,
                 stats->dPrimalIterHistory, double, stats->nCol);
    HDSDP_MEMCPY(stats->dPrimalIterHistory, dColVal, double, stats->nCol);
    stats->dPrimalMuHistory[iIter] = dMu;
    
    /* Compute convergence statistics */
    
    return;
}

static int HPrimalStatsSuperlinerTest( hdsdp_primal_stats *stats ) {
    
    int isSuperLin = 0;
    
    return isSuperLin;
}

static void HPrimalStatsClear( hdsdp_primal_stats *stats ) {
    
    if ( !stats ) {
        return;
    }
    
    HDSDP_FREE(stats->dPrimalIterHistory);
    HDSDP_FREE(stats->dPrimalMuHistory);
    HDSDP_ZERO(stats, hdsdp_primal_stats, 1);

    return;
}

static void HPrimalStatsDestroy( hdsdp_primal_stats **pStats ) {
    
    if ( !pStats ) {
        return;
    }
    
    HPrimalStatsClear(*pStats);
    HDSDP_FREE(*pStats);
    
    return;
}

static hdsdp_lpsolver_params HLpSolverIGetDefaultParams(void) {
    
    hdsdp_lpsolver_params params;
    
    /* Optimization tolerance */
    params.dAbsOptTol = 1e-06;
    params.dAbsFeasTol = 1e-06;
    params.dRelOptTol = 1e-10;
    params.dRelFeasTol = 1e-10;
    
    params.dKKTPrimalReg = 1e-14;
    params.dKKTDualReg = 1e-12;
    
    params.dPotentialRho = 2.0;
    params.dPrimalUpdateStep = 0.95;
    params.dDualUpdateStep = 0.95;
    params.dIterativeTol = 1e-12;
    params.dScalingThreshTol = 1e-03;
    params.dBarrierLowerBnd = 1e-16;
    
    params.nThreads = 8;
    params.nScalIter = 1;
    params.nMaxIter = 1000;
    
    params.iStartMethod = MEHROTRA_START;
    params.iScalMethod = SCAL_GEOMETRIC;
    params.iPrimalMethod = 1;
    
    params.LpMethod = LP_ITER_MEHROTRA;
    
    return params;
}

static hdsdp_retcode HLpSolverIScaleData( hdsdp_lpsolver *HLp ) {
    /* Perform ruiz and geometric scaling on the IPM data */
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    if ( HLp->params.iScalMethod == SCAL_RUIZ ) {
        csp_ruizscal(HLp->nRow, HLp->nCol, HLp->colMatBeg, HLp->colMatIdx, HLp->colMatElem,
                     HLp->dRowScal, HLp->dColScal, HLp->params.nScalIter);
    } else if ( HLp->params.iScalMethod == SCAL_GEOMETRIC ) {
        csp_geoscal(HLp->nRow, HLp->nCol, HLp->colMatBeg, HLp->colMatIdx, HLp->colMatElem, HLp->dRowScal, HLp->dColScal);
    } else {
        csp_l2scal(HLp->nRow, HLp->nCol, HLp->colMatBeg, HLp->colMatIdx, HLp->colMatElem, HLp->dRowScal, HLp->dColScal);
    }

    /* Scale transpose matrix */
    csp_colscal(HLp->nRow, HLp->colMatTransBeg, HLp->colMatTransIdx, HLp->colMatTransElem, HLp->dRowScal);
    csp_rowscal(HLp->nRow, HLp->colMatTransBeg, HLp->colMatTransIdx, HLp->colMatTransElem, HLp->dColScal);
    
    /* Scale b and c */
    vvrscl(&HLp->nRow, HLp->dRowScal, HLp->dRowRhs);
    vvrscl(&HLp->nCol, HLp->dColScal, HLp->dColObj);
    
exit_cleanup:
    return retcode;
}

static hdsdp_retcode HLpSolverIComputeMehrotraStartingPoint( hdsdp_lpsolver *HLp ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    /* Factorize A * A' */
    for ( int iCol = 0; iCol < HLp->nCol; ++iCol ) {
        HLp->dAuxiColVector1[iCol] = 1.0;
    }
    
    HDSDP_CALL(HLpKKTSetup(HLp->Hkkt, LP_ITER_MEHROTRA, HLp->dAuxiColVector1, 0.0, 0.0));
    
    /* Compute A * c */
    csp_Axpy(HLp->nCol, HLp->colMatBeg, HLp->colMatIdx, HLp->colMatElem, 1.0, HLp->dColObj, HLp->dAuxiColVector2);
    
    /* Solve least squares to get x and y */
    /* x = A' * (AAT \ b) */
    HDSDP_CALL(HLpKKTSolveNormalEqn(HLp->Hkkt, 1, HLp->dRowRhs, HLp->dAuxiRowVector1));
    csp_ATxpy(HLp->nCol, HLp->colMatBeg, HLp->colMatIdx, HLp->colMatElem,
               1.0, HLp->dAuxiRowVector1, HLp->dColVal);
    /* y = AAT \ (A * c) */
    HDSDP_CALL(HLpKKTSolveNormalEqn(HLp->Hkkt, 1, HLp->dAuxiColVector2, HLp->dRowDual));
    /* s = c - A' * y */
    HDSDP_MEMCPY(HLp->dColDual, HLp->dColObj, double, HLp->nCol);
    csp_ATxpy(HLp->nCol, HLp->colMatBeg, HLp->colMatIdx, HLp->colMatElem, -1.0, HLp->dRowDual, HLp->dColDual);
    
    /* Compute coordinate shift */
    int iXmin = idmin(&HLp->nCol, HLp->dColVal, &HIntConstantOne);
    int iSmin = idmin(&HLp->nCol, HLp->dColDual, &HIntConstantOne);
    
    double dDeltaX = HDSDP_MAX(-1.5 * HLp->dColVal[iXmin], 0.0);
    double dDeltaS = HDSDP_MAX(-1.5 * HLp->dColDual[iSmin], 0.0);
    
    double dDeltaXS = 0.0;
    double dDeltaXC = 0.0;
    double dDeltaSC = 0.0;
    double dSumX = 0.0;
    double dSumS = 0.0;
    
    for ( int iCol = 0; iCol < HLp->nCol; ++iCol ) {
        dDeltaXS += (HLp->dColVal[iCol] + dDeltaX) * (HLp->dColDual[iCol] + dDeltaS);
        dSumX += HLp->dColVal[iCol];
        dSumS += HLp->dColDual[iCol];
    }
    
    dDeltaXS = dDeltaXS * 0.5;
    dDeltaXC = dDeltaX + dDeltaXS / ( dSumS + HLp->nCol * dDeltaS );
    dDeltaSC = dDeltaS + dDeltaXS / ( dSumX + HLp->nCol * dDeltaX );
    
    for ( int iCol = 0; iCol < HLp->nCol; ++iCol ) {
        HLp->dColVal[iCol] += dDeltaXC;
        HLp->dColDual[iCol] += dDeltaSC;
    }
    
    /* Done with Mehrotra's starting point */
exit_cleanup:
    return retcode;
}

static hdsdp_retcode HLpSolverIInitialize( hdsdp_lpsolver *HLp, int iStartMethod ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    if ( iStartMethod == TRIVIAL_START ) {
        for ( int iCol = 0; iCol < HLp->nCol; ++iCol ) {
            HLp->dColVal[iCol] = 1e+02;
            HLp->dColDual[iCol] = 1e+02;
        }
    } else {
        /* Compute Mehrotra starting point */
        HDSDP_CALL(HLpSolverIComputeMehrotraStartingPoint(HLp));
    }
    
exit_cleanup:
    return retcode;
}

static void HLpSolverISetupResidual( hdsdp_lpsolver *HLp ) {
    
    /* Get A * x - b * tau */
    csp_Axpy(HLp->nCol, HLp->colMatBeg, HLp->colMatIdx, HLp->colMatElem,
              1.0, HLp->dColVal, HLp->dPrimalInfeasVec);
    
    double dMinusTau = - HLp->dTau;
    axpy(&HLp->nCol, &dMinusTau, HLp->dRowRhs, &HIntConstantOne,
         HLp->dPrimalInfeasVec, &HIntConstantOne);
    
    /* Get A' * y + s - c * tau */
    HDSDP_MEMCPY(HLp->dDualInfeasVec, HLp->dColDual, double, HLp->nCol);
    csp_ATxpy(HLp->nCol, HLp->colMatBeg, HLp->colMatIdx, HLp->colMatElem,
               1.0, HLp->dRowDual, HLp->dDualInfeasVec);
    axpy(&HLp->nCol, &dMinusTau, HLp->dColObj, &HIntConstantOne, HLp->dDualInfeasVec, &HIntConstantOne);
    
    /* Get complementarity XSe */
    for ( int iCol = 0; iCol < HLp->nCol; ++iCol ) {
        HLp->dComplVec[iCol] = HLp->dColVal[iCol] * HLp->dColDual[iCol];
    }

    return;
}

static int HLpSolverICollectStats( hdsdp_lpsolver *HLp ) {
    
    return 1;
}

static hdsdp_retcode HLpSolverITakePrimalDualStep( hdsdp_lpsolver *HLp ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    
    
exit_cleanup:
    return retcode;
}

static hdsdp_retcode HLpSolverITakeHSDStep( hdsdp_lpsolver *HLp ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    assert( 0 );
    
    
exit_cleanup:
    return retcode;
}

static hdsdp_retcode HLpSolverITakePrimalStep( hdsdp_lpsolver *HLp ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
exit_cleanup:
    return retcode;
}

extern hdsdp_retcode HLpSolverCreate( hdsdp_lpsolver **pHLp ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    if ( !pHLp ) {
        retcode = HDSDP_RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    hdsdp_lpsolver *HLp = NULL;
    HDSDP_INIT(HLp, hdsdp_lpsolver, 1);
    HDSDP_MEMCHECK(HLp);
    HDSDP_ZERO(HLp, hdsdp_lpsolver, 1);
    *pHLp = HLp;
    
exit_cleanup:
    return retcode;
}

extern hdsdp_retcode HLpSolverInit( hdsdp_lpsolver *HLp, int nRow, int nCol ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    HLp->nRow = nRow;
    HLp->nCol = nCol;
    
    HDSDP_INIT(HLp->dRowRhs, double, nRow);
    HDSDP_INIT(HLp->dColObj, double, nCol);
    HDSDP_INIT(HLp->dRowScal, double, nRow);
    HDSDP_INIT(HLp->dColScal, double, nCol);
    
    HDSDP_MEMCHECK(HLp->dRowRhs);
    HDSDP_MEMCHECK(HLp->dColObj);
    HDSDP_MEMCHECK(HLp->dRowScal);
    HDSDP_MEMCHECK(HLp->dColScal);
    
    HLp->params = HLpSolverIGetDefaultParams();
    HDSDP_CALL(HPrimalStatsCreate(&HLp->stats));
    HDSDP_CALL(HPrimalStatsInit(HLp->stats, nCol, &HLp->params));
    
    HLp->LpStatus = HDSDP_UNKNOWN;
    
    HDSDP_INIT(HLp->dColVal, double, nCol);
    HDSDP_INIT(HLp->dColDual, double, nCol);
    HDSDP_INIT(HLp->dRowDual, double, nRow);
    
    HDSDP_MEMCHECK(HLp->dColVal);
    HDSDP_MEMCHECK(HLp->dColDual);
    HDSDP_MEMCHECK(HLp->dRowDual);
    
    HDSDP_INIT(HLp->dPrimalInfeasVec, double, nRow);
    HDSDP_INIT(HLp->dDualInfeasVec, double, nCol);
    HDSDP_INIT(HLp->dComplVec, double, nCol);
    HDSDP_INIT(HLp->dScalingMatrix, double, nCol);
    
    HDSDP_MEMCHECK(HLp->dPrimalInfeasVec);
    HDSDP_MEMCHECK(HLp->dDualInfeasVec);
    HDSDP_MEMCHECK(HLp->dComplVec);
    HDSDP_MEMCHECK(HLp->dScalingMatrix);
    
    HDSDP_INIT(HLp->dColValDirection, double, nCol);
    HDSDP_INIT(HLp->dColDualDirection, double, nCol);
    HDSDP_INIT(HLp->dRowDualDirection, double, nRow);
    
    HDSDP_MEMCHECK(HLp->dColValDirection);
    HDSDP_MEMCHECK(HLp->dColDualDirection);
    HDSDP_MEMCHECK(HLp->dRowDualDirection);
    
    HDSDP_INIT(HLp->dAuxiColVector1, double, nCol);
    HDSDP_INIT(HLp->dAuxiColVector2, double, nCol);
    HDSDP_INIT(HLp->dAuxiColVector3, double, nCol);
    
    HDSDP_MEMCHECK(HLp->dAuxiColVector1);
    HDSDP_MEMCHECK(HLp->dAuxiColVector2);
    HDSDP_MEMCHECK(HLp->dAuxiColVector3);
    
    HDSDP_INIT(HLp->dAuxiRowVector1, double, nRow);
    HDSDP_INIT(HLp->dAuxiRowVector2, double, nRow);
    
    HDSDP_MEMCHECK(HLp->dAuxiRowVector1);
    HDSDP_MEMCHECK(HLp->dAuxiRowVector2);
    
    HDSDP_CALL(HLpKKTCreate(&HLp->Hkkt));
    
    HLp->dBarrierMu = HDSDP_INFINITY;
    
    HLp->pObjVal = HDSDP_INFINITY;
    HLp->dObjVal = - HDSDP_INFINITY;
    HLp->dPrimalDualGap = HDSDP_INFINITY;
    
    HLp->dKappa = 1.0;
    HLp->dTau = 1.0;
    HLp->dKappaDirection = 0.0;
    HLp->dTauDirection = 0.0;
    
exit_cleanup:
    return retcode;
}

extern hdsdp_retcode HLpSolverSetData( hdsdp_lpsolver *HLp, int *colMatBeg, int *colMatIdx, double *colMatElem,
                                       int *colMatTransBeg, int *colMatTransIdx, double *colMatTransElem, double *rowRHS,
                                       double *colObj ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    HLp->nElem = colMatBeg[HLp->nCol];
    
    HDSDP_INIT(HLp->colMatBeg, int, HLp->nCol + 1);
    HDSDP_INIT(HLp->colMatIdx, int, HLp->nElem);
    HDSDP_INIT(HLp->colMatElem, double, HLp->nElem);
    
    HDSDP_MEMCHECK(HLp->colMatBeg);
    HDSDP_MEMCHECK(HLp->colMatIdx);
    HDSDP_MEMCHECK(HLp->colMatElem);
    
    HDSDP_INIT(HLp->colMatTransBeg, int, HLp->nRow + 1);
    HDSDP_INIT(HLp->colMatTransIdx, int, HLp->nElem);
    HDSDP_INIT(HLp->colMatTransElem, double, HLp->nElem);
    
    HDSDP_MEMCHECK(HLp->colMatTransBeg);
    HDSDP_MEMCHECK(HLp->colMatTransIdx);
    HDSDP_MEMCHECK(HLp->colMatTransElem);
    
    /* Copy problem data */
    HDSDP_MEMCPY(HLp->colMatBeg, colMatBeg, int, HLp->nCol + 1);
    HDSDP_MEMCPY(HLp->colMatIdx, colMatIdx, int, HLp->nElem);
    HDSDP_MEMCPY(HLp->colMatElem, colMatElem, double, HLp->nElem);
    
    HDSDP_MEMCPY(HLp->colMatTransBeg, colMatTransBeg, int, HLp->nRow + 1);
    HDSDP_MEMCPY(HLp->colMatTransIdx, colMatTransIdx, int, HLp->nElem);
    HDSDP_MEMCPY(HLp->colMatTransElem, colMatTransElem, double, HLp->nElem);
    
    HDSDP_MEMCPY(HLp->dRowRhs, rowRHS, double, HLp->nRow);
    HDSDP_MEMCPY(HLp->dColObj, colObj, double, HLp->nCol);
    
    /* Initialize KKT solver */
    HDSDP_CALL(HLpKKTInit(HLp->Hkkt, HLp->nRow, HLp->nCol,
                          HLp->colMatBeg, HLp->colMatIdx, HLp->colMatElem,
                          HLp->colMatTransBeg, HLp->colMatTransIdx, HLp->colMatTransElem));
    
exit_cleanup:
    return retcode;
}

extern hdsdp_retcode HLpSolverOptimize( hdsdp_lpsolver *HLp ) {
    /* Invoke hybrid primal-dual and primal solver */
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    /* Compute starting point */
    HDSDP_CALL(HLpSolverIScaleData(HLp));
    HDSDP_CALL(HLpSolverIInitialize(HLp, HLp->params.iStartMethod));
    
    int nMaxIter = HLp->params.nMaxIter;
    int goOn = 0;
    
    for ( int nIter = 0; nIter < nMaxIter; ++nIter ) {
        
        if ( HLp->params.LpMethod == LP_ITER_PRIMAL ) {
            HDSDP_CALL(HLpSolverITakePrimalStep(HLp));
        } else if ( HLp->params.LpMethod == LP_ITER_HSD ) {
            HDSDP_CALL(HLpSolverITakeHSDStep(HLp));
        } else {
            HDSDP_CALL(HLpSolverITakePrimalDualStep(HLp));
        }
        
        goOn = HLpSolverICollectStats(HLp);
        
        if ( !goOn ) {
            break;
        }
    }
    
exit_cleanup:
    return retcode;
}

extern void HLpSolverClear( hdsdp_lpsolver *HLp ) {
    
    if ( !HLp ) {
        return;
    }
    
    HDSDP_FREE(HLp->colMatBeg);
    HDSDP_FREE(HLp->colMatIdx);
    HDSDP_FREE(HLp->colMatElem);
    
    HDSDP_FREE(HLp->colMatTransBeg);
    HDSDP_FREE(HLp->colMatTransIdx);
    HDSDP_FREE(HLp->colMatTransElem);
    
    HDSDP_FREE(HLp->dRowRhs);
    HDSDP_FREE(HLp->dColObj);
    HDSDP_FREE(HLp->dRowScal);
    HDSDP_FREE(HLp->dColScal);
    
    HPrimalStatsDestroy(&HLp->stats);
    
    HDSDP_FREE(HLp->dColVal);
    HDSDP_FREE(HLp->dColDual);
    HDSDP_FREE(HLp->dRowDual);
    
    HDSDP_FREE(HLp->dPrimalInfeasVec);
    HDSDP_FREE(HLp->dDualInfeasVec);
    HDSDP_FREE(HLp->dComplVec);
    HDSDP_FREE(HLp->dScalingMatrix);
    
    HDSDP_FREE(HLp->dColValDirection);
    HDSDP_FREE(HLp->dColDualDirection);
    HDSDP_FREE(HLp->dRowDualDirection);
    
    HDSDP_FREE(HLp->dAuxiColVector1);
    HDSDP_FREE(HLp->dAuxiColVector2);
    HDSDP_FREE(HLp->dAuxiColVector3);
    
    HDSDP_FREE(HLp->dAuxiRowVector1);
    HDSDP_FREE(HLp->dAuxiRowVector2);
    
    HLpKKTDestroy(&HLp->Hkkt);
    
    HDSDP_ZERO(HLp, hdsdp_lpsolver, 1);
    
    return;
}

extern void HLpSolverDestroy( hdsdp_lpsolver **pHLp ) {
    
    if ( !pHLp ) {
        return;
    }
    
    HLpSolverClear(*pHLp);
    HDSDP_FREE(*pHLp);
}
