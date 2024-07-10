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
    stats->dCondNumberEst = HDSDP_INFINITY;
    stats->dIterDiffMetric = HDSDP_INFINITY;
    stats->dIterDiffMetricScal = HDSDP_INFINITY;
    stats->dIterDiffMetricThresh = HDSDP_INFINITY;
    stats->dIterDiffMetricAggressive = HDSDP_INFINITY;
    
    HDSDP_INIT(stats->dPrimalMuHistory, double, params->nMaxIter + 1);
    HDSDP_INIT(stats->dPrimalIterHistory, double, nCol);
    HDSDP_MEMCHECK(stats->dPrimalIterHistory);
    
exit_cleanup:
    return retcode;
}

static void HPrimalStatsUpdate( hdsdp_primal_stats *stats, int iIter, double *dColVal, double dMu ) {
    
    stats->dPrimalMuHistory[iIter] = dMu;
    stats->iIterLast = iIter;
    
    if ( iIter <= 1 ) {
        HDSDP_MEMCPY(stats->dPrimalIterHistory, dColVal, double, stats->nCol);
        return;
    }
    
    /* Compute convergence statistics */
    stats->dIterDiffMetric = 0.0;
    stats->dIterDiffMetricScal = 0.0;
    stats->dIterDiffMetricThresh = 0.0;
    
    /* Aggressive statistic ignores difference in the small-magnitude elements*/
    stats->dIterDiffMetricAggressive = 0.0;
    
    /* Compute distance between consecutive iterates. Use infinity norm instead of 2 norm */
    for ( int iCol = 0; iCol < stats->nCol; ++iCol ) {
        
        /* More aggressive when a history point is used */
        double dColValElem = stats->dPrimalIterHistory[iCol];
        double dColValDiff = fabs(dColVal[iCol] - stats->dPrimalIterHistory[iCol]);
        double dColValScalDiff = dColValDiff / dColValElem;
        
        stats->dIterDiffMetric = HDSDP_MAX(dColValDiff, stats->dIterDiffMetric);
        stats->dIterDiffMetricScal = HDSDP_MAX(dColValScalDiff, stats->dIterDiffMetricScal);
        
        /* Compute thresholded difference */
        if ( dColVal[iCol] > stats->params->dScalingThreshTol ) {
            stats->dIterDiffMetricThresh = HDSDP_MAX(dColValScalDiff, stats->dIterDiffMetricThresh);
            stats->dIterDiffMetricAggressive = HDSDP_MAX(dColValScalDiff, stats->dIterDiffMetricAggressive);
            if ( stats->dIterDiffMetricAggressive >= 0.8 ) {
                // printf("Here \n");
            }
            
        } else {
            stats->dIterDiffMetricThresh = HDSDP_MAX(dColValDiff, stats->dIterDiffMetricThresh);
        }
    }
    
    if ( stats->dIterDiffMetricThresh < 1.0 ) {
        stats->dCondNumberEst = (1 + stats->dIterDiffMetricThresh) / (1 - stats->dIterDiffMetricThresh);
        stats->dCondNumberEst = stats->dCondNumberEst * stats->dCondNumberEst;
    } else {
        stats->dCondNumberEst = HDSDP_INFINITY;
    }
    
    /* Overwrite history */
    HDSDP_MEMCPY(stats->dPrimalIterHistory, dColVal, double, stats->nCol);
    
    return;
}

#ifndef LSUPERTEST
#define LSUPERTEST (5)
#endif
static void HPrimalStatsSuperlinerTest( hdsdp_primal_stats *stats ) {
    
    /* Take most recent 5 iterates */
    int isSuperLin = 0;
    int nSuperTest = LSUPERTEST;
    double *dMuHist = stats->dPrimalMuHistory;
    
    if ( stats->iIterLast < 3 ) {
        stats->isSuperLin = 0;
        return;
    }
    
    if ( stats->iIterLast <= nSuperTest * 2 ) {
        nSuperTest = stats->iIterLast / 2;
    }
    
    double dRecentAvg = 0.0;
    double dHistAvg = 0.0;
    
    for ( int iElem = 0; iElem < nSuperTest; ++iElem ) {
        dRecentAvg += log(dMuHist[stats->iIterLast - iElem]) - log(dMuHist[stats->iIterLast - iElem - 1]);
        dHistAvg += log(dMuHist[stats->iIterLast - nSuperTest - iElem]) - log(dMuHist[stats->iIterLast - nSuperTest - iElem - 1]);
    }
    
    if ( dRecentAvg < dHistAvg ) {
        isSuperLin = 1;
    }
    
    stats->isSuperLin = isSuperLin;
    return;
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
    params.dAbsOptTol = 1.0;
    params.dAbsFeasTol = 1.0;
    params.dRelOptTol = 1e-12;
    params.dRelFeasTol = 1e-12;
    
    params.dKKTPrimalReg = 1e-14;
    params.dKKTDualReg = 1e-12;
    
    params.dPotentialRho = 100.0;
    params.dPrimalUpdateStep = 0.995;
    params.dDualUpdateStep = 0.995;
    params.dIterativeTol = 1e-12;
    params.dScalingThreshTol = 1e-04;
    params.dBarrierLowerBndCoeff = 1e-03;
    params.dTimeLimit = 3600.0;
    
    params.nThreads = 8;
    params.nScalIter = 10;
    params.nMaxIter = 100;
    
    params.iStartMethod = MEHROTRA_START;
    params.iScalMethod = SCAL_GEOMETRIC;
    params.iPrimalMethod = 0;
    
    params.LpMethod = LP_ITER_MEHROTRA;
    
    return params;
}

static void HLpSolverICollectLpStats( hdsdp_lpsolver *HLp ) {
    
    /* Get statistics */
    HLp->lpstats.dObjInfNorm = 0.0;
    HLp->lpstats.dObjTwoNorm = 0.0;
    HLp->lpstats.dObjOneNorm = 0.0;
    
    HLp->lpstats.dRhsInfNorm = 0.0;
    HLp->lpstats.dRhsTwoNorm = 0.0;
    HLp->lpstats.dRhsOneNorm = 0.0;
    
    for ( int iCol = 0; iCol < HLp->nCol; ++iCol ) {
        double dAbsElem = fabs(HLp->dColObj[iCol]);
        HLp->lpstats.dObjInfNorm = HDSDP_MAX(dAbsElem, HLp->lpstats.dObjInfNorm);
        HLp->lpstats.dObjOneNorm += dAbsElem;
        HLp->lpstats.dObjTwoNorm += dAbsElem * dAbsElem;
    }
    
    HLp->lpstats.dObjTwoNorm = sqrt(HLp->lpstats.dObjTwoNorm);
    
    for ( int iRow = 0; iRow < HLp->nRow; ++iRow ) {
        double dAbsElem = fabs(HLp->dRowRhs[iRow]);
        HLp->lpstats.dRhsInfNorm = HDSDP_MAX(dAbsElem, HLp->lpstats.dRhsInfNorm);
        HLp->lpstats.dRhsOneNorm += dAbsElem;
        HLp->lpstats.dRhsTwoNorm += dAbsElem * dAbsElem;
    }
    
    HLp->lpstats.dRhsTwoNorm = sqrt(HLp->lpstats.dRhsTwoNorm);
    HLp->lpstats.dAMatFNorm = 0.0;
    
    for ( int iElem = 0; iElem < HLp->nElem; ++iElem ) {
        HLp->lpstats.dAMatAbsNorm += fabs(HLp->colMatElem[iElem]);
        HLp->lpstats.dAMatFNorm += HLp->colMatElem[iElem] * HLp->colMatElem[iElem];
    }
    
    HLp->lpstats.dAMatFNorm = sqrt(HLp->lpstats.dAMatFNorm);
    HLp->lpstats.nAMatNz = HLp->colMatBeg[HLp->nCol];
    
    if ( HLp->lpstats.dObjInfNorm < 1e-08 ) {
        HLp->lpstats.isNoObj = 1;
    }
    
    if ( HLp->lpstats.dRhsInfNorm < 1e-08 ) {
        HLp->lpstats.isNoRhs = 1;
    }
    
    HLp->params.dScalingThreshTol = 1e-03 / HLp->lpstats.dAMatFNorm;

    /* Adjust solver parameters */
    if ( HLp->lpstats.isNoObj ) {
        HLp->params.dScalingThreshTol = 1e-03 / sqrt(HLp->nCol);
        HLp->params.dBarrierLowerBndCoeff = 1e-05;
    }
    
    return;
}

static hdsdp_retcode HLpSolverIScaleData( hdsdp_lpsolver *HLp ) {
    /* Perform ruiz and geometric scaling on the IPM data */
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    for ( int iElem = 0; iElem < HLp->nCol; ++iElem ) {
        HLp->dColScal[iElem] = 1.0;
    }
    for ( int iElem = 0; iElem < HLp->nRow; ++iElem ) {
        HLp->dRowScal[iElem] = 1.0;
    }
    
    if ( HLp->params.iScalMethod == SCAL_RUIZ ) {
        HDSDP_CALL(csp_ruizscal(HLp->nRow, HLp->nCol, HLp->colMatBeg, HLp->colMatIdx, HLp->colMatElem,
                                HLp->dRowScal, HLp->dColScal, HLp->params.nScalIter));
    } else if ( HLp->params.iScalMethod == SCAL_GEOMETRIC ) {
        csp_geoscal(HLp->nRow, HLp->nCol, HLp->colMatBeg, HLp->colMatIdx, HLp->colMatElem, HLp->dRowScal, HLp->dColScal);
    } else if ( HLp->params.iScalMethod == SCAL_L2NORM ) {
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
    
    double dsNorm = nrm1(&HLp->nCol, HLp->dColDual, &HIntConstantOne);
    if ( dsNorm < 1e-08 ) {
        for ( int iElem = 0; iElem < HLp->nCol; ++iElem ) {
            HLp->dColDual[iElem] += 1.0;
        }
    }
    
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
    
    HLp->dBarrierMu = dot(&HLp->nCol, HLp->dColVal, &HIntConstantOne, HLp->dColDual, &HIntConstantOne);
    HLp->dBarrierMu = HLp->dBarrierMu / HLp->nCol;

exit_cleanup:
    return retcode;
}

static hdsdp_retcode HLpSolverIComputeStartingPoint( hdsdp_lpsolver *HLp, int iStartMethod ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    /* Compute starting point */
    if ( iStartMethod == TRIVIAL_START ) {
        for ( int iCol = 0; iCol < HLp->nCol; ++iCol ) {
            HLp->dColVal[iCol] = 1e+02;
            HLp->dColDual[iCol] = 1e+02;
        }
        HLp->dBarrierMu = 1e+02;
    } else {
        /* Compute Mehrotra starting point */
        HDSDP_CALL(HLpSolverIComputeMehrotraStartingPoint(HLp));
    }
    
exit_cleanup:
    return retcode;
}

static void HLpSolverISetupResidual( hdsdp_lpsolver *HLp ) {
    
    /* Get A * x - b * tau */
    HDSDP_ZERO(HLp->dPrimalInfeasVec, double, HLp->nRow);
    csp_Axpy(HLp->nCol, HLp->colMatBeg, HLp->colMatIdx, HLp->colMatElem,
              1.0, HLp->dColVal, HLp->dPrimalInfeasVec);
    
    double dMinusTau = - HLp->dTau;
    axpy(&HLp->nRow, &dMinusTau, HLp->dRowRhs, &HIntConstantOne,
         HLp->dPrimalInfeasVec, &HIntConstantOne);
    
    /* Get A' * y + s - c * tau */
    HDSDP_MEMCPY(HLp->dDualInfeasVec, HLp->dColDual, double, HLp->nCol);
    csp_ATxpy(HLp->nCol, HLp->colMatBeg, HLp->colMatIdx, HLp->colMatElem,
               1.0, HLp->dRowDual, HLp->dDualInfeasVec);
    axpy(&HLp->nCol, &dMinusTau, HLp->dColObj, &HIntConstantOne, HLp->dDualInfeasVec, &HIntConstantOne);

    return;
}

static int HLpSolverIComputeSolutionStats( hdsdp_lpsolver *HLp, int iIter ) {
    
    int goOn = 1;
    
    if ( iIter == 0 ) {
        hdsdp_printf("Using Hybrid Primal-Primal-Dual solver \n\n");
        hdsdp_printf("    %5s  %15s  %15s  %8s  %8s  %8s  %10s   %5s \n",
                     "nIter", "pObj", "dObj", "pInf", "dInf", "Mu", "P/D Step", "T [PD]");
    }
    
    HLpSolverISetupResidual(HLp);
    
    /* Primal infeasibility */
    HLp->dPrimalInfeas = 0.0;
    for ( int iRow = 0; iRow < HLp->nRow; ++iRow ) {
        double dScaledInfeas = HLp->dPrimalInfeasVec[iRow] * HLp->dRowScal[iRow];
        HLp->dPrimalInfeas +=  dScaledInfeas * dScaledInfeas;
    }
    
    HLp->dPrimalInfeas = sqrt(HLp->dPrimalInfeas);
    HLp->dPrimalInfeasRel = HLp->dPrimalInfeas / (HLp->lpstats.dRhsTwoNorm + 1.0);
    
    /* Dual infeasibility */
    HLp->dDualInfeas = 0.0;
    for ( int iCol = 0; iCol < HLp->nCol; ++iCol ) {
        double dScaledInfeas = HLp->dDualInfeasVec[iCol] * HLp->dColScal[iCol];
        HLp->dDualInfeas += dScaledInfeas * dScaledInfeas;
    }
    
    HLp->dDualInfeas = sqrt(HLp->dDualInfeas);
    HLp->dDualInfeasRel = HLp->dDualInfeas / (HLp->lpstats.dObjTwoNorm + 1.0);
    
    /* Duality gap */
    HLp->pObjVal = dot(&HLp->nCol, HLp->dColObj, &HIntConstantOne, HLp->dColVal, &HIntConstantOne);
    HLp->dObjVal = dot(&HLp->nRow, HLp->dRowRhs, &HIntConstantOne, HLp->dRowDual, &HIntConstantOne);
    HLp->dPrimalDualGap = fabs(HLp->pObjVal - HLp->dObjVal);
    HLp->dPrimalDualGapRel = fabs(HLp->pObjVal - HLp->dObjVal) / (fabs(HLp->pObjVal) + fabs(HLp->dObjVal) + 1.0);
    
    /* Compute current barrier parameter */
    double dCompl = dot(&HLp->nCol, HLp->dColVal, &HIntConstantOne, HLp->dColDual, &HIntConstantOne);
    dCompl = dCompl / HLp->nCol;
    
    if ( (HLp->dPrimalDualGapRel <= HLp->params.dRelOptTol)  &&
         (HLp->dPrimalInfeasRel  <= HLp->params.dRelFeasTol) &&
         (HLp->dDualInfeasRel    <= HLp->params.dRelFeasTol) &&
         (HLp->dPrimalDualGap    <= HLp->params.dAbsOptTol)  &&
         (HLp->dPrimalInfeas     <= HLp->params.dAbsFeasTol) &&
         (HLp->dDualInfeas       <= HLp->params.dAbsFeasTol) ) {
        
        goOn = 0;
    }
    
    /* Logging */
#ifdef HDSDP_PIPMMETRIC_DEBUG
    hdsdp_printf("    %5d  %+15.8e  %+15.8e  %8.2e  %8.2e  %8.2e  %5.2f %5.2f  %4.1f  | B: %8.2e | P: %5.1e | F/S: %5.1f"
                 " Cond: %6.1e | L2: %6.1e | SL2: %6.1e | Agg: %6.1e [%d] \n", iIter, HLp->pObjVal, HLp->dObjVal, HLp->dPrimalInfeasRel,
                 HLp->dDualInfeasRel, dCompl, HLp->pStep, HLp->dStep, HUtilGetTimeStamp() - HLp->dTStart, HLp->dBarrierMu,
                 HLp->dProxNorm, HLpKKTGetFactorSolveTimeRatio(HLp->Hkkt), HLp->pstats->dCondNumberEst, HLp->pstats->dIterDiffMetric,
                 HLp->pstats->dIterDiffMetricThresh, HLp->pstats->dIterDiffMetricAggressive, HLp->pstats->isSuperLin);
#else
    hdsdp_printf("    %5d  %+15.8e  %+15.8e  %8.2e  %8.2e  %5.2f %5.2f  %4.1f \n",
                 iIter, HLp->pObjVal, HLp->dObjVal, HLp->dPrimalInfeasRel, HLp->dDualInfeasRel,
           HLp->pStep, HLp->dStep, HUtilGetTimeStamp() - HLp->dTStart);
#endif
    
    return goOn;
}

static int HLpSolverICheckPrimalStats( hdsdp_lpsolver *HLp, int iIter) {
    
    int iPrimalStart = 0;
    
    HPrimalStatsUpdate(HLp->pstats, iIter, HLp->dColVal, HLp->dBarrierMu);
    
    if ( iIter == 0 ) {
        return 0;
    }
    
    HPrimalStatsSuperlinerTest(HLp->pstats);
    
    /* Determine whether to switch to primal IPM based on convergence statistics */
    double dCondUbEst = HLp->pstats->dCondNumberEst;
    double dEuclideanDist = HLp->pstats->dIterDiffMetric;
//    double dScalDist = HLp->pstats->dIterDiffMetricScal;
//    double dThreshDist = HLp->pstats->dIterDiffMetricThresh;
//    double dThreshAggDist = HLp->pstats->dIterDiffMetricAggressive;
    
    int iPrimalStartCond1 = 0;
    int iPrimalStartCond2 = 0;
    
    iPrimalStartCond1 = HLp->params.iPrimalMethod && (HLp->params.LpMethod != LP_ITER_PRIMAL);
    iPrimalStartCond2 = (dCondUbEst < 100.0 || dEuclideanDist < 1e-05) && \
                        (HLp->dPrimalDualGapRel < 1e-03 && HLp->dPrimalDualGapRel > HLp->params.dRelOptTol * 1e+02);
    
    /* Relaxed condition for debugging */
//    iPrimalStartCond2 = (dCondUbEst < 100.0 || dEuclideanDist < dAdaTol);
    
    if ( iPrimalStartCond1 && iPrimalStartCond2 ) {// && !HLp->pstats->isSuperLin ) {
        iPrimalStart = 1;
    }
    
    return iPrimalStart;
}

static double HLpSolverISingleRatioTest( int nElem, double *dVal, double *dStep ) {
    
    double dRatio = HDSDP_INFINITY;
    
    for ( int iElem = 0; iElem < nElem; ++iElem ) {
        double dTmp = dStep[iElem] / dVal[iElem];
        dRatio = HDSDP_MIN(dRatio, dTmp);
    }
    
    if ( dRatio > 0.0 ) {
        return 100.0;
    }
    
    return 1.0 / fabs(dRatio);
}

static void HLpSolverIRatioTest( hdsdp_lpsolver *HLp, double *dPrimalStep, double *dDualStep ) {
    
    *dPrimalStep = HLpSolverISingleRatioTest(HLp->nCol, HLp->dColVal, HLp->dColValDirection);
    *dDualStep = HLpSolverISingleRatioTest(HLp->nCol, HLp->dColDual, HLp->dColDualDirection);
    
    return;
}

/* Primal-dual interior point method main algorithm. Use Mehrotra's predictor-corrector step */
static hdsdp_retcode HLpSolverITakePrimalDualStep( hdsdp_lpsolver *HLp ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    double *dXSinv = HLp->dAuxiColVector1;
    double *dColTmp = HLp->dAuxiColVector2;
    double *dColTmp2 = HLp->dAuxiColVector3;
    double *dXSinvRd = HLp->dAuxiColVector4;
    double *dNormalEqnRhsTmp = HLp->dAuxiRowVector1;
    double pStepCorr = 0.0;
    double dStepCorr = 0.0;
    double dBarrierAffine = 0.0;
    double dBarrierTarget = 0.0;
    
    /* Set up the scaling matrix */
    for ( int iCol = 0; iCol < HLp->nCol; ++iCol ) {
        dColTmp[iCol] = sqrt(HLp->dColVal[iCol] / HLp->dColDual[iCol]);
        dXSinv[iCol] = dColTmp[iCol] * dColTmp[iCol];
    }
    
    /* Setup and factorize normal matrix */
    HDSDP_CALL(HLpKKTSetup(HLp->Hkkt, LP_ITER_MEHROTRA, dColTmp,
                HLp->params.dKKTPrimalReg, HLp->params.dKKTDualReg));
    
    /* Predictor step */
    // rtmp = A * x - rp - A * (xsinv .* rd) = b - A * (xsinv .* rd);
    HDSDP_MEMCPY(dNormalEqnRhsTmp, HLp->dRowRhs, double, HLp->nRow);
    HDSDP_MEMCPY(dXSinvRd, HLp->dDualInfeasVec, double, HLp->nCol);
    vvscl(&HLp->nCol, dXSinv, dXSinvRd);
    csp_Axpy(HLp->nCol, HLp->colMatBeg, HLp->colMatIdx, HLp->colMatElem, -1.0, dXSinvRd, dNormalEqnRhsTmp);
    
    /* dycorr = info.factor \ rtmp; */
    HLpKKTSolveNormalEqn(HLp->Hkkt, 1, dNormalEqnRhsTmp, HLp->dRowDualDirection);
    /* dxcorr = xsinv .* rd - x + xsinv .* (A' * dycorr); */
    /* xsinv .* rd */
    HDSDP_MEMCPY(HLp->dColValDirection, dXSinvRd, double, HLp->nCol);
    /* A' * dycorr */
    HDSDP_ZERO(dColTmp, double, HLp->nCol);
    csp_ATxpy(HLp->nCol, HLp->colMatBeg, HLp->colMatIdx, HLp->colMatElem, 1.0, HLp->dRowDualDirection, dColTmp);
    /* xsinv .* (A' * dycorr) */
    vvscl(&HLp->nCol, dXSinv, dColTmp);
    /*  xsinv .* rd - x  */
    axpy(&HLp->nCol, &HDblConstantMinusOne, HLp->dColVal, &HIntConstantOne, HLp->dColValDirection, &HIntConstantOne);
    /* xsinv .* rd - x + xsinv .* (A' * dycorr); */
    axpy(&HLp->nCol, &HDblConstantOne, dColTmp, &HIntConstantOne, HLp->dColValDirection, &HIntConstantOne);
    
    /* dscorr = - s - dxcorr ./ xsinv; */
    for ( int iCol = 0; iCol < HLp->nCol; ++iCol ) {
        HLp->dColDualDirection[iCol] = -HLp->dColDual[iCol] - HLp->dColValDirection[iCol] / dXSinv[iCol];
    }
    
    /* Done with predictor step */
    
    /* Ratio test */
    HLpSolverIRatioTest(HLp, &pStepCorr, &dStepCorr);
    pStepCorr = HDSDP_MIN(pStepCorr, 1.0);
    dStepCorr = HDSDP_MIN(dStepCorr, 1.0);
    
    /* Estimate next barrier parameter */
    for ( int iElem = 0; iElem < HLp->nCol; ++iElem ) {
        dBarrierAffine += (HLp->dColVal[iElem] + pStepCorr * HLp->dColValDirection[iElem]) \
                            * (HLp->dColDual[iElem] + dStepCorr * HLp->dColDualDirection[iElem]);
    }
    dBarrierAffine = dBarrierAffine / HLp->nCol;
    dBarrierTarget = (dBarrierAffine / HLp->dBarrierMu);
    dBarrierTarget = HLp->dBarrierMu * dBarrierTarget * dBarrierTarget * dBarrierTarget;
    double dBarrierLowerBnd = HLp->params.dRelFeasTol * HLp->params.dBarrierLowerBndCoeff;
    dBarrierTarget = HDSDP_MAX(dBarrierTarget, dBarrierLowerBnd);
    HLp->dBarrierMu = HDSDP_MIN(dBarrierTarget, HLp->dBarrierMu);
    
    /* Corrector step */
    /* rmu = x .* s + dxcorr .* dscorr - mu; */
    for ( int iElem = 0; iElem < HLp->nCol; ++iElem ) {
        dColTmp[iElem] = HLp->dColVal[iElem] * HLp->dColDual[iElem] + \
                            HLp->dColValDirection[iElem] * HLp->dColDualDirection[iElem] - HLp->dBarrierMu;
        dColTmp2[iElem] = dColTmp[iElem] / HLp->dColDual[iElem];
    }
    
    /* rtmp = A * (rmu ./ s) - rp - A * (xsinv .* rd); */
    HDSDP_ZERO(dNormalEqnRhsTmp, double, HLp->nRow);
    csp_Axpy(HLp->nCol, HLp->colMatBeg, HLp->colMatIdx, HLp->colMatElem, 1.0, dColTmp2, dNormalEqnRhsTmp);
    csp_Axpy(HLp->nCol, HLp->colMatBeg, HLp->colMatIdx, HLp->colMatElem, -1.0, dXSinvRd, dNormalEqnRhsTmp);
    axpy(&HLp->nRow, &HDblConstantMinusOne, HLp->dPrimalInfeasVec, &HIntConstantOne, dNormalEqnRhsTmp, &HIntConstantOne);
    
    /* dy = info.factor \ rtmp; */
    HDSDP_CALL(HLpKKTSolveNormalEqn(HLp->Hkkt, 1, dNormalEqnRhsTmp, HLp->dRowDualDirection));
    /* dx = xsinv .* rd - rmu ./ s + xsinv .* (A' * dy); */
    HDSDP_MEMCPY(HLp->dColValDirection, dXSinvRd, double, HLp->nCol);
    /* xsinv .* rd - rmu ./ s */
    for ( int iElem = 0; iElem < HLp->nCol; ++iElem ) {
        /* dx */
        HLp->dColValDirection[iElem] -= dColTmp[iElem] / HLp->dColDual[iElem];
        /* ds */
        HLp->dColDualDirection[iElem] = - dColTmp[iElem] / HLp->dColVal[iElem];
    }
    /* Set up A' * dy and dColTmp = rmu is reused */
    HDSDP_ZERO(dColTmp, double, HLp->nCol);
    csp_ATxpy(HLp->nCol, HLp->colMatBeg, HLp->colMatIdx, HLp->colMatElem, 1.0, HLp->dRowDualDirection, dColTmp);
    /* xsinv .* (A' * dy) */
    vvscl(&HLp->nCol, dXSinv, dColTmp);
    axpy(&HLp->nCol, &HDblConstantOne, dColTmp, &HIntConstantOne, HLp->dColValDirection, &HIntConstantOne);
    
    /* ds = - rmu ./ x - dx ./ xsinv; */
    for ( int iElem = 0; iElem < HLp->nCol; ++iElem ) {
        HLp->dColDualDirection[iElem] -= HLp->dColValDirection[iElem] / dXSinv[iElem];
    }
    
    /* Done with corrector step */
    HLpSolverIRatioTest(HLp, &HLp->pStep, &HLp->dStep);
    HLp->pStep = HDSDP_MIN(HLp->params.dPrimalUpdateStep * HLp->pStep, 1.0);
    HLp->dStep = HDSDP_MIN(HLp->params.dDualUpdateStep * HLp->dStep, 1.0);
    
    /* Take step into next iteration */
    axpy(&HLp->nCol, &HLp->pStep, HLp->dColValDirection, &HIntConstantOne, HLp->dColVal, &HIntConstantOne);
    axpy(&HLp->nCol, &HLp->dStep, HLp->dColDualDirection, &HIntConstantOne, HLp->dColDual, &HIntConstantOne);
    axpy(&HLp->nRow, &HLp->dStep, HLp->dRowDualDirection, &HIntConstantOne, HLp->dRowDual, &HIntConstantOne);
    
    HLp->dBarrierMu = dot(&HLp->nCol, HLp->dColVal, &HIntConstantOne, HLp->dColDual, &HIntConstantOne);
    HLp->dBarrierMu = HLp->dBarrierMu / HLp->nCol;
    HLp->dBarrierMu = HDSDP_MAX(HLp->dBarrierMu, dBarrierLowerBnd);
    
exit_cleanup:
    return retcode;
}

static hdsdp_retcode HLpSolverIPreparePrimal( hdsdp_lpsolver *HLp ) {
    /* Prepare the factorization in primal IPM */
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    /* Overwrite with debug data */
    
#if 0
#include "debug_data_lp.h"
    assert( ndbgCol == HLp->nCol && ndbgRow == HLp->nRow );
    HDSDP_MEMCPY(HLp->dColVal, dbgColVal, double, ndbgCol);
    HDSDP_MEMCPY(HLp->dRowDual, dbgRowDual, double, ndbgRow);
    HDSDP_MEMCPY(HLp->dColDual, dbgColDual, double, ndbgCol);
    HLp->dBarrierMu = dbgBarrierParam;
#endif
    /* Decide the scaling matrix */
    double *dScalingMatrix = HLp->dScalingMatrix;
    HDSDP_MEMCPY(dScalingMatrix, HLp->dColVal, double, HLp->nCol);
    
    /* Predict the elements close to 0 to avoid tiny stepsize caused by inexactness */
    if ( 0 ) {
        for ( int iElem = 0; iElem < HLp->nCol; ++iElem ) {
            if ( dScalingMatrix[iElem] < 1e-05 ) {
                dScalingMatrix[iElem] = dScalingMatrix[iElem] * 0.1;
            }
        }
    }
    
    /* Set up the normal matrix. In the current implementation, this is the only factorization
       in primal interior point method */
    HDSDP_CALL(HLpKKTSetup(HLp->Hkkt, LP_ITER_PRIMAL, dScalingMatrix,
                           HLp->params.dKKTPrimalReg, HLp->params.dKKTDualReg));
    
    /* Adjust primal parameters */
    HLp->params.dPrimalUpdateStep = 0.95;
    HLp->params.dDualUpdateStep = 0.95;
    HLp->params.LpMethod = LP_ITER_PRIMAL;
    
exit_cleanup:
    return retcode;
}

static double HLpSolverIExtrapolatePotentialParam( hdsdp_lpsolver *HLp ) {
    /* Use extrapolation to get the next barrier parameter.  */
    double dPotentialRho = HLp->params.dPotentialRho;
    int iIter = HLp->pstats->iIterLast;
    
    dPotentialRho = HDSDP_MAX(dPotentialRho, HLp->pstats->dPrimalMuHistory[iIter - 1] / HLp->pstats->dPrimalMuHistory[iIter]);
    
    return dPotentialRho;
}

static void HLpSolverIMatVec( hdsdp_lpsolver *HLp, double dXCoeff, double *dXVec, double *dYVec ) {
    
    if ( 0 ) {
        // HLpKKTMultiply(HLp->Hkkt, dXCoeff, dXVec, dYVec);
    } else {
        /* Compute y = y + coeff * A * X^2 * A' * x  */
        HDSDP_ZERO(HLp->dKrylovAuxVec6, double, HLp->nCol);
        csp_ATxpy(HLp->nCol, HLp->colMatBeg, HLp->colMatIdx, HLp->colMatElem, 1.0, dXVec, HLp->dKrylovAuxVec6);
        for ( int iElem = 0; iElem < HLp->nCol; ++iElem ) {
            HLp->dKrylovAuxVec6[iElem] *= HLp->dAuxiColVector4[iElem];
        }
        csp_Axpy(HLp->nCol, HLp->colMatBeg, HLp->colMatIdx, HLp->colMatElem, dXCoeff, HLp->dKrylovAuxVec6, dYVec);
    }
    
    for ( int iElem = 0; iElem < HLp->nRow; ++iElem ) {
        dYVec[iElem] += HLp->params.dKKTDualReg * dXVec[iElem];
    }
    
    return;
}

static int HLpSolverITestIterativeDir( hdsdp_lpsolver *HLp, double *dRowDualDirection ) {
    
    /* Shifted scaling matrix */
    double *dShiftedScalingMatrix = HLp->dAuxiColVector1;
    /* Shifted scaling error matrix */
    double *dErrScalMatrix = HLp->dAuxiColVector2;
    double *dColTmp = HLp->dAuxiColVector3;
    double *dShiftedScalingMatrixSqr = HLp->dAuxiColVector4;
    double *dTrialRowDualDirection = HLp->dAuxiRowVector2;
    double *dTrialPrimalInfeasVec = HLp->dPrimalInfeasVec;
    double dTrialPrimalInfeas = 0.0;
    double dTrialPrimalInfeasRel = 0.0;
    
    /* Copy solution to trial space */
    HDSDP_MEMCPY(dTrialRowDualDirection, dRowDualDirection, double, HLp->nRow);
    
    /* Get dy */
    scal(&HLp->nRow, &HLp->dBarrierMu, dTrialRowDualDirection, &HIntConstantOne);
    
    /* Get ds = -rd - A' * dy */
    HDSDP_ZERO(HLp->dColDualDirection, double, HLp->nCol);
    csp_ATxpy(HLp->nCol, HLp->colMatBeg, HLp->colMatIdx, HLp->colMatElem,
              -1.0, dTrialRowDualDirection, HLp->dColDualDirection);
    axpy(&HLp->nCol, &HDblConstantMinusOne, HLp->dDualInfeasVec, &HIntConstantOne,
         HLp->dColDualDirection, &HIntConstantOne);
    
    /* Get dx = vxinv .* v - (vsqr .* (s + ds)) / mu */
    for ( int iElem = 0; iElem < HLp->nCol; ++iElem ) {
        HLp->dColValDirection[iElem] = dErrScalMatrix[iElem] * dShiftedScalingMatrix[iElem] - \
        (dShiftedScalingMatrixSqr[iElem] * (HLp->dColDual[iElem] + HLp->dColDualDirection[iElem])) / HLp->dBarrierMu;
    }
    
    HLpSolverIRatioTest(HLp, &HLp->pStep, &HLp->dStep);
    
    if ( HLp->pStep < 1e-04 ) {
        hdsdp_printf("=> Iterative test ends with tiny stepsize \n");
        return 0;
    }
    
    /* Test for feasibility of primal */
    HDSDP_MEMCPY(dColTmp, HLp->dColVal, double, HLp->nCol);
    axpy(&HLp->nCol, &HLp->pStep, HLp->dColValDirection, &HIntConstantOne, dColTmp, &HIntConstantOne);
    HDSDP_MEMCPY(dTrialPrimalInfeasVec, HLp->dRowRhs, double, HLp->nRow);
    csp_Axpy(HLp->nCol, HLp->colMatBeg, HLp->colMatIdx, HLp->colMatElem, -1.0, dColTmp, dTrialPrimalInfeasVec);
    
    for ( int iRow = 0; iRow < HLp->nRow; ++iRow ) {
        double dScaledInfeas = dTrialPrimalInfeasVec[iRow] * HLp->dRowScal[iRow];
        dTrialPrimalInfeas += dScaledInfeas * dScaledInfeas;
    }
    
    dTrialPrimalInfeas = sqrt(dTrialPrimalInfeas);
    dTrialPrimalInfeasRel = dTrialPrimalInfeas / (HLp->lpstats.dRhsTwoNorm + 1.0);
    
//    hdsdp_printf("=> Iterative test: pStep %5.2f dStep %5.2f pInf: %5.1e \n",
//                 HLp->pStep, HLp->dStep, dTrialPrimalInfeasRel);
    
    if ( dTrialPrimalInfeasRel > HLp->dPrimalInfeasRel && dTrialPrimalInfeasRel > HLp->params.dRelOptTol ) {
        return 0;
    }
    
    return 1;
}

static hdsdp_retcode HLpSolverIConjGrad( hdsdp_lpsolver *HLp, int nMaxIter, double dTol, double *dRhs, double *dSol ) {
    /* Apply preconditioned conjugate gradient method to solve the normal matrix
     
     The iterative solver may interact with the LP problem. Specifically, in every inner iteration
     we recover the primal direction dx and stop if
     
        1. Ratio test gives sufficiently large stepsize
        2. || A * (x + alpha_x * dx ) - b || is reduced
     
     This idea was exploited in
     Zanetti, F., & Gondzio, J. (2023). A new stopping criterion for Krylov solvers applied in interior point methods.
     SIAM Journal on Scientific Computing, 45(2), A703-A728.
     
     */
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    int nCol = HLp->nRow;
    
    double dDotMTimesd = 0.0;
    double resiDotPreInvResi = 0.0;
    double alpha = 0.0;
    double beta = 0.0;
    double resiNorm = 0.0;
    double rhsNorm = 0.0;
    
    double *iterResi = HLp->dKrylovAuxVec1;
    double *iterResiNew = HLp->dKrylovAuxVec2;
    double *iterDirection = HLp->dKrylovAuxVec3;
    double *MTimesDirection = HLp->dKrylovAuxVec4;
    double *preInvResi = HLp->dKrylovAuxVec5;
    double *iterVec = dSol;
    double *rhsVec = dRhs;
    
    int iter = 0;
    int iPassTest = 0;
    
    rhsNorm = nrm2(&nCol, rhsVec, &HIntConstantOne);
    dTol = dTol * rhsNorm;
    
    /* Initial guess */
    HDSDP_CALL(HLpKKTSolveNormalEqn(HLp->Hkkt, 1, rhsVec, iterVec));
    
    HDSDP_MEMCPY(iterResi, rhsVec, double, nCol);
    HLpSolverIMatVec(HLp, -1.0, iterVec, iterResi);
    
    resiNorm = nrm2(&nCol, iterResi, &HIntConstantOne);
    
    iPassTest = HLpSolverITestIterativeDir(HLp, iterVec);
    
    if ( resiNorm < dTol ) {
        goto exit_cleanup;
    }
    
    HDSDP_MEMCPY(iterDirection, iterResi, double, nCol);
    HDSDP_CALL(HLpKKTSolveNormalEqn(HLp->Hkkt, 1, iterDirection, NULL));
    HDSDP_MEMCPY(preInvResi, iterDirection, double, nCol);
    
    HDSDP_ZERO(MTimesDirection, double, nCol);
    HLpSolverIMatVec(HLp, 1.0, iterDirection, MTimesDirection);
    
    for ( iter = 0; iter < nMaxIter; ++iter ) {
        
        resiDotPreInvResi = dot(&nCol, preInvResi, &HIntConstantOne, iterResi, &HIntConstantOne);
        dDotMTimesd = dot(&nCol, iterDirection, &HIntConstantOne, MTimesDirection, &HIntConstantOne);
        alpha = resiDotPreInvResi / dDotMTimesd;
        axpy(&nCol, &alpha, iterDirection, &HIntConstantOne, iterVec, &HIntConstantOne);
        
        HDSDP_MEMCPY(iterResiNew, iterResi, double, nCol);
        double minusAlpha = - alpha;
        
        axpy(&nCol, &minusAlpha, MTimesDirection, &HIntConstantOne, iterResiNew, &HIntConstantOne);
        HDSDP_MEMCPY(preInvResi, iterResiNew, double, nCol);
        HDSDP_CALL(HLpKKTSolveNormalEqn(HLp->Hkkt, 1, preInvResi, NULL));
        beta = dot(&nCol, iterResiNew, &HIntConstantOne, preInvResi, &HIntConstantOne);
        beta = beta / resiDotPreInvResi;
        axpby(&nCol, &HDblConstantOne, preInvResi, &HIntConstantOne, &beta, iterDirection, &HIntConstantOne);
        HDSDP_ZERO(MTimesDirection, double, nCol);
        HLpSolverIMatVec(HLp, 1.0, iterDirection, MTimesDirection);
        
        // HDSDP_MEMCPY(iterResi, iterResiNew, double, nCol);
        HDSDP_MEMCPY(iterResi, rhsVec, double, nCol);
        HLpSolverIMatVec(HLp, -1.0, iterVec, iterResi);
        resiNorm = nrm2(&nCol, iterResi, &HIntConstantOne);
        iPassTest = HLpSolverITestIterativeDir(HLp, iterVec);
        
#ifdef HDSDP_CONJGRAD_DEBUG
        hdsdp_printf("CG Iteration: %d %6.2e \n", iter, resiNorm);
#endif
        
        if ( resiNorm != resiNorm ) {
            hdsdp_printf("Iterative solver failed \n");
            retcode = HDSDP_RETCODE_FAILED;
            goto exit_cleanup;
        }
        
        if ( resiNorm < dTol ) {
            break;
        }
    }
    
exit_cleanup:
    return retcode;
}

static hdsdp_retcode HLpSolverIIterativeSolve( hdsdp_lpsolver *HLp, double *dRhs, double *dSol ) {
    
    /* Apply preconditioned iterative method to solve normal system */
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    int nKrylovIter = 0;
    double dFactorSolveRatio = HLpKKTGetFactorSolveTimeRatio(HLp->Hkkt);
    
    /* Determine maximum iterations based on solution statistics */
    if ( dFactorSolveRatio < 50 ) {
        nKrylovIter = 2;
    } else if ( dFactorSolveRatio < 100 ) {
        nKrylovIter = 3;
    } else {
        nKrylovIter = 5;
    }
    
    double dIterativeTol = HLp->params.dIterativeTol;
    HDSDP_CALL(HLpSolverIConjGrad(HLp, nKrylovIter, dIterativeTol, dRhs, dSol));
    
exit_cleanup:
    return retcode;
}

/* Primal interior point method main algorithm */
static hdsdp_retcode HLpSolverITakePrimalStep( hdsdp_lpsolver *HLp ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    double dPotentialRho = HLpSolverIExtrapolatePotentialParam(HLp);
    double dShiftedThreshTol = HLp->params.dScalingThreshTol;
    
    /* Scaling matrix of the factorized normal matrix */
    double *dScalingMatrix = HLp->dScalingMatrix;
    /* Shifted scaling matrix */
    double *dShiftedScalingMatrix = HLp->dAuxiColVector1;
    /* Shifted scaling error matrix */
    double *dErrScalMatrix = HLp->dAuxiColVector2;
    double *dColTmp = HLp->dAuxiColVector3;
    double *dShiftedScalingMatrixSqr = HLp->dAuxiColVector4;
    double *dNormalEqnRhsTmp = HLp->dAuxiRowVector1;
    
    /* Compute shifted scaling matrix */
    for ( int iElem = 0; iElem < HLp->nCol; ++iElem ) {
        if ( HLp->dColVal[iElem] < dShiftedThreshTol ) {
            /* Small element. Let scaling matrix synchronize this part */
            dShiftedScalingMatrix[iElem] = HLp->dColVal[iElem];
            dErrScalMatrix[iElem] = 1.0;
        } else {
            /* Large element. Let scaling matrix keep unchanged */
            dShiftedScalingMatrix[iElem] = dScalingMatrix[iElem];
            dErrScalMatrix[iElem] = dShiftedScalingMatrix[iElem] / HLp->dColVal[iElem];
        }
        dShiftedScalingMatrixSqr[iElem] = dShiftedScalingMatrix[iElem] * dShiftedScalingMatrix[iElem];
    }
    
    /*  Set up the RHS of the shifted nromal system
     rtmp = -rp + A * (v .* ((v .* s) / mu - v./x)) - A * (vsqr .* (rd / mu)); */
    HDSDP_ZERO(dNormalEqnRhsTmp, double, HLp->nRow);
    for ( int iElem = 0; iElem < HLp->nCol; ++iElem ) {
        dColTmp[iElem] = (dShiftedScalingMatrix[iElem] * HLp->dColDual[iElem]) / HLp->dBarrierMu \
                          - dErrScalMatrix[iElem];
        dColTmp[iElem] = dColTmp[iElem] * dShiftedScalingMatrix[iElem];
    }
    csp_Axpy(HLp->nCol, HLp->colMatBeg, HLp->colMatIdx, HLp->colMatElem, 1.0, dColTmp, dNormalEqnRhsTmp);
    
    /* Only consider dual infeasibility when it is larger than the tightened tolerance */
    if ( HLp->dDualInfeasRel > HLp->params.dRelFeasTol * 1e-02 ) {
        for ( int iElem = 0; iElem < HLp->nCol; ++iElem ) {
            dColTmp[iElem] = dShiftedScalingMatrixSqr[iElem] * HLp->dDualInfeasVec[iElem] / HLp->dBarrierMu;
        }
        csp_Axpy(HLp->nCol, HLp->colMatBeg, HLp->colMatIdx, HLp->colMatElem, -1.0, dColTmp, dNormalEqnRhsTmp);
    }
    
    axpy(&HLp->nRow, &HDblConstantMinusOne, HLp->dPrimalInfeasVec, &HIntConstantOne, 
         dNormalEqnRhsTmp, &HIntConstantOne);
    
    /* Solve for dy / mu. Note that iterative solver and LP solver are integrated in the
       primal IPM solver: ratio test is done on each iterative of the inner iterative solver
       Hence when iterative solver returns, all the step directions and ratio test are available  */
    HDSDP_CALL(HLpSolverIIterativeSolve(HLp, dNormalEqnRhsTmp, HLp->dRowDualDirection));
    
    /* Recover other directions */
    scal(&HLp->nRow, &HLp->dBarrierMu, HLp->dRowDualDirection, &HIntConstantOne);
    HDSDP_ZERO(HLp->dColDualDirection, double, HLp->nCol);
    csp_ATxpy(HLp->nCol, HLp->colMatBeg, HLp->colMatIdx, HLp->colMatElem,
              -1.0, HLp->dRowDualDirection, HLp->dColDualDirection);
    axpy(&HLp->nCol, &HDblConstantMinusOne, HLp->dDualInfeasVec, &HIntConstantOne,
         HLp->dColDualDirection, &HIntConstantOne);
    for ( int iElem = 0; iElem < HLp->nCol; ++iElem ) {
        HLp->dColValDirection[iElem] = dErrScalMatrix[iElem] * dShiftedScalingMatrix[iElem] - \
        (dShiftedScalingMatrixSqr[iElem] * (HLp->dColDual[iElem] + HLp->dColDualDirection[iElem])) / HLp->dBarrierMu;
    }
    
    /* Final ratio text */
    HLpSolverIRatioTest(HLp, &HLp->pStep, &HLp->dStep);
    HLp->pStep = HDSDP_MIN(HLp->params.dPrimalUpdateStep * HLp->pStep, 1.0);
    HLp->dStep = HDSDP_MIN(HLp->params.dDualUpdateStep * HLp->dStep, 1.0);
    
    /* Go to next iteration */
    axpy(&HLp->nCol, &HLp->pStep, HLp->dColValDirection, &HIntConstantOne, HLp->dColVal, &HIntConstantOne);
    axpy(&HLp->nCol, &HLp->dStep, HLp->dColDualDirection, &HIntConstantOne, HLp->dColDual, &HIntConstantOne);
    axpy(&HLp->nRow, &HLp->dStep, HLp->dRowDualDirection, &HIntConstantOne, HLp->dRowDual, &HIntConstantOne);
    
    /* Additional ratio test on the dual problem */
    int isDualFeas = 1;
    HDSDP_MEMCPY(dColTmp, HLp->dColObj, double, HLp->nCol);
    csp_ATxpy(HLp->nCol, HLp->colMatBeg, HLp->colMatIdx, HLp->colMatElem, -1.0, HLp->dRowDual, dColTmp);
    
    for ( int iElem = 0; iElem < HLp->nCol; ++iElem ) {
        if ( dColTmp[iElem] < 0.0 ) {
            isDualFeas = 0;
            break;
        }
    }
    
    double dBarrierTarget = 0.0;
    double dBarrierLowerBnd = 0.0;
    double dComplGap = 0.0;
    double dXiSi = 0.0;
    
    /* Compute next barrier parameter */
    if ( isDualFeas ) {
        HDSDP_MEMCPY(HLp->dColDual, dColTmp, double, HLp->nCol);
        dBarrierTarget = dot(&HLp->nCol, HLp->dColVal, &HIntConstantOne, HLp->dColDual, &HIntConstantOne);
        dBarrierTarget = dBarrierTarget / dPotentialRho;
        dBarrierTarget = HDSDP_MIN(dBarrierTarget, HLp->dBarrierMu);
    } else {
        double dBarrierStep = HDSDP_MIN(HLp->pStep, HLp->dStep);
        dBarrierStep = HDSDP_MIN(dBarrierStep, 0.6);
        dBarrierTarget = HLp->dBarrierMu * (1.0 - dBarrierStep);
    }
    
    /* Compute the complementarity measure */
    for ( int iElem = 0; iElem < HLp->nCol; ++iElem ) {
        dXiSi = HLp->dColVal[iElem] * HLp->dColDual[iElem];
        HLp->dComplVec[iElem] = dXiSi;
        dComplGap += dXiSi;
    }
    
    dComplGap = dComplGap / HLp->nCol;
    
    /* Threshold barrier parameter */
    dBarrierLowerBnd = dComplGap / 10.0;
    dBarrierTarget = HDSDP_MAX(dBarrierTarget, dBarrierLowerBnd);
    
    /* Compute the proximity measure */
    HLp->dProxNorm = 0.0;
    
    for ( int iElem = 0; iElem < HLp->nCol; ++iElem ) {
        double dComplDeviate = fabs(HLp->dComplVec[iElem] / dComplGap - 1.0);
        HLp->dProxNorm = HDSDP_MAX(dComplDeviate, HLp->dProxNorm);
    }
    
    /* Dynamically adjust barrier parameter according to proximity norm */
    if ( HLp->dProxNorm < 1.0 ) {
        dBarrierTarget = dBarrierTarget * 0.3;
    }
    
    if ( (HLp->dProxNorm > 100) && (HLp->dPrimalInfeasRel > HLp->params.dRelFeasTol) ) {
        dBarrierTarget = HDSDP_MIN(HLp->dBarrierMu, dComplGap);
    }
    
    dBarrierLowerBnd = HLp->params.dRelFeasTol * HLp->params.dBarrierLowerBndCoeff;
    HLp->dBarrierMu = HDSDP_MAX(dBarrierTarget, dBarrierLowerBnd);
    
exit_cleanup:
    return retcode;
}

/* Homogeneous self-dual model with Mehrotra's predictor-corrector step */
static hdsdp_retcode HLpSolverITakeHSDStep( hdsdp_lpsolver *HLp ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    assert( 0 );
    
    
exit_cleanup:
    return retcode;
}

static void HDSDPIPrintSolutionStats( hdsdp_lpsolver *HLp ) {
    
    if ( HLp->LpStatus == HDSDP_UNKNOWN ) {
        hdsdp_printf("\nLP Status: %s \n", "Unknown");
    } else if ( HLp->LpStatus == HDSDP_DUAL_FEASIBLE ) {
        hdsdp_printf("\nLP Status: %s \n", "Dual feasible");
    } else if ( HLp->LpStatus == HDSDP_DUAL_OPTIMAL ) {
        hdsdp_printf("\nLP Status: %s \n", "Dual optimal ");
    } else if ( HLp->LpStatus == HDSDP_PRIMAL_DUAL_OPTIMAL ) {
        hdsdp_printf("\nLP Status: %s \n", "Primal dual optimal");
    } else if ( HLp->LpStatus == HDSDP_MAXITER ) {
        hdsdp_printf("\nLP Status: %s \n", "Maximum iteration");
    } else if ( HLp->LpStatus == HDSDP_SUSPECT_INFEAS_OR_UNBOUNDED ) {
        hdsdp_printf("\nLP Status: %s \n", "Suspected infeasible or unbounded");
    } else if ( HLp->LpStatus == HDSDP_INFEAS_OR_UNBOUNDED ) {
        hdsdp_printf("\nLP Status: %s \n", "Infeasible or unbounded");
    } else if ( HLp->LpStatus == HDSDP_TIMELIMIT ) {
        hdsdp_printf("\nLP Status: %s \n", "Time limit");
    } else if ( HLp->LpStatus == HDSDP_USER_INTERRUPT ) {
        hdsdp_printf("\nLP Status: %s \n", "User interrupt");
    } else if ( HLp->LpStatus == HDSDP_INTERNAL_ERROR ) {
        hdsdp_printf("\nLP Status: %s \n", "Internal error");
    } else if ( HLp->LpStatus == HDSDP_NUMERICAL ) {
        hdsdp_printf("\nLP Status: %s \n", "Numerical error");
    } else {
        assert( 0 );
    }
    
    hdsdp_printf("\nLP Solution statistic \n");
    hdsdp_printf("pObj: %+8.3e   dObj: %+8.3e \n", HLp->pObjVal, HLp->dObjVal);
    hdsdp_printf("Abs. pInf:  %8.3e    Rel. pInf:  %8.3e \n", HLp->dPrimalInfeas, HLp->dPrimalInfeasRel);
    hdsdp_printf("Abs. dInf:  %8.3e    Rel. dInf:  %8.3e \n", HLp->dDualInfeas, HLp->dDualInfeasRel);
    hdsdp_printf("Abs. Gap :  %8.3e    Rel. gap :  %8.3e \n", HLp->dPrimalDualGap, HLp->dPrimalDualGapRel);
    
    hdsdp_printf("\nLinear solver statistic \n");
    
    /* Remove solve in starting point */
    hdsdp_printf("Num. Factor: %3d (+1)    Factor Time: %5.2f \n", HLp->Hkkt->nFactor - 1, HLp->Hkkt->dFactorTime);
    hdsdp_printf("Num. Solves: %3d (+1)    Solve  Time: %5.2f \n", HLp->Hkkt->nSolve - 1, HLp->Hkkt->dSolveTime);
    
    hdsdp_printf("\nElapsed Time:  %5.3f seconds \n\n", HUtilGetTimeStamp() - HLp->dTStart);
    return;
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
    
    HDSDP_CALL(HPrimalStatsCreate(&HLp->pstats));
    HDSDP_CALL(HPrimalStatsInit(HLp->pstats, nCol, &HLp->params));
    
    HLp->LpStatus = HDSDP_UNKNOWN;
    
    HDSDP_INIT(HLp->dColVal, double, nCol);
    HDSDP_INIT(HLp->dColDual, double, nCol);
    HDSDP_INIT(HLp->dRowDual, double, nRow);
    
    HDSDP_MEMCHECK(HLp->dColVal);
    HDSDP_MEMCHECK(HLp->dColDual);
    HDSDP_MEMCHECK(HLp->dRowDual);
    
    HDSDP_INIT(HLp->dPrimalInfeasVec, double, nRow);
    HDSDP_INIT(HLp->dDualInfeasVec, double, nCol);
    HDSDP_INIT(HLp->dScalingMatrix, double, nCol);
    HDSDP_INIT(HLp->dComplVec, double, nCol);
    
    HDSDP_MEMCHECK(HLp->dPrimalInfeasVec);
    HDSDP_MEMCHECK(HLp->dDualInfeasVec);
    HDSDP_MEMCHECK(HLp->dScalingMatrix);
    HDSDP_MEMCHECK(HLp->dComplVec);
    
    HDSDP_INIT(HLp->dColValDirection, double, nCol);
    HDSDP_INIT(HLp->dColDualDirection, double, nCol);
    HDSDP_INIT(HLp->dRowDualDirection, double, nRow);
    
    HDSDP_MEMCHECK(HLp->dColValDirection);
    HDSDP_MEMCHECK(HLp->dColDualDirection);
    HDSDP_MEMCHECK(HLp->dRowDualDirection);
    
    HDSDP_INIT(HLp->dAuxiColVector1, double, nCol);
    HDSDP_INIT(HLp->dAuxiColVector2, double, nCol);
    HDSDP_INIT(HLp->dAuxiColVector3, double, nCol);
    HDSDP_INIT(HLp->dAuxiColVector4, double, nCol);
    
    HDSDP_MEMCHECK(HLp->dAuxiColVector1);
    HDSDP_MEMCHECK(HLp->dAuxiColVector2);
    HDSDP_MEMCHECK(HLp->dAuxiColVector3);
    HDSDP_MEMCHECK(HLp->dAuxiColVector4);
    
    HDSDP_INIT(HLp->dAuxiRowVector1, double, nRow);
    HDSDP_INIT(HLp->dAuxiRowVector2, double, nRow);
    
    HDSDP_MEMCHECK(HLp->dAuxiRowVector1);
    HDSDP_MEMCHECK(HLp->dAuxiRowVector2);
    
    HDSDP_INIT(HLp->dKrylovAuxVec1, double, nRow);
    HDSDP_INIT(HLp->dKrylovAuxVec2, double, nRow);
    HDSDP_INIT(HLp->dKrylovAuxVec3, double, nRow);
    HDSDP_INIT(HLp->dKrylovAuxVec4, double, nRow);
    HDSDP_INIT(HLp->dKrylovAuxVec5, double, nRow);
    HDSDP_INIT(HLp->dKrylovAuxVec6, double, nCol);
    
    HDSDP_MEMCHECK(HLp->dKrylovAuxVec1);
    HDSDP_MEMCHECK(HLp->dKrylovAuxVec2);
    HDSDP_MEMCHECK(HLp->dKrylovAuxVec3);
    HDSDP_MEMCHECK(HLp->dKrylovAuxVec4);
    HDSDP_MEMCHECK(HLp->dKrylovAuxVec5);
    HDSDP_MEMCHECK(HLp->dKrylovAuxVec6);
    
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
    
    HLp->dTStart = HUtilGetTimeStamp();
    
    /* Scaled problem data */
    HLpSolverICollectLpStats(HLp);
    retcode = HLpSolverIScaleData(HLp);
    
    if ( retcode != HDSDP_RETCODE_OK ) {
        HLp->LpStatus = HDSDP_INTERNAL_ERROR;
        goto exit_cleanup;
    }
    
    /* Compute starting point */
    retcode = HLpSolverIComputeStartingPoint(HLp, HLp->params.iStartMethod);
    
    if ( retcode != HDSDP_RETCODE_OK ) {
        HLp->LpStatus = HDSDP_INTERNAL_ERROR;
        goto exit_cleanup;
    }
    
    int nMaxIter = HLp->params.nMaxIter;
    int goOn = 0;
    int iPrimalStart = 0;
    
    hdsdp_printf("Opimizing an LP of %d variables and %d constraints\n", HLp->nCol, HLp->nRow);
    hdsdp_printf("Data statistics: |A| = %5.2e |b| = %5.2e |c| = %5.2e Nnz = %d \n",
                 HLp->lpstats.dAMatAbsNorm, HLp->lpstats.dRhsOneNorm, HLp->lpstats.dObjOneNorm, HLp->lpstats.nAMatNz);
    
    HLpSolverIComputeSolutionStats(HLp, 0);
    
    for ( int nIter = 1; nIter <= nMaxIter; ++nIter ) {
        
        if ( HLp->params.LpMethod == LP_ITER_PRIMAL ) {
            retcode = HLpSolverITakePrimalStep(HLp);
        } else {
            retcode = HLpSolverITakePrimalDualStep(HLp);
        }
        
        if ( retcode != HDSDP_RETCODE_OK ) {
            HLp->LpStatus = HDSDP_NUMERICAL;
            break;
        }
        
        goOn = HLpSolverIComputeSolutionStats(HLp, nIter);
        
        if ( HLp->dPrimalDualGap != HLp->dPrimalDualGap ) {
            HLp->LpStatus = HDSDP_NUMERICAL;
            break;
        }
        
        if ( HUtilGetTimeStamp() - HLp->dTStart > HLp->params.dTimeLimit ) {
            HLp->LpStatus = HDSDP_TIMELIMIT;
            break;
        }
        
        if ( !goOn ) {
            HLp->LpStatus = HDSDP_PRIMAL_DUAL_OPTIMAL;
            break;
        }
        
        iPrimalStart = HLpSolverICheckPrimalStats(HLp, nIter);
        
        if ( iPrimalStart ) {
            hdsdp_printf("Primal interior point method starts \n");
            HDSDP_CALL(HLpSolverIPreparePrimal(HLp));
            iPrimalStart = 0;
        }
    }
    
exit_cleanup:
    
    HDSDPIPrintSolutionStats(HLp);
    
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
    
    HPrimalStatsDestroy(&HLp->pstats);
    
    HDSDP_FREE(HLp->dColVal);
    HDSDP_FREE(HLp->dColDual);
    HDSDP_FREE(HLp->dRowDual);
    
    HDSDP_FREE(HLp->dPrimalInfeasVec);
    HDSDP_FREE(HLp->dDualInfeasVec);
    HDSDP_FREE(HLp->dScalingMatrix);
    HDSDP_FREE(HLp->dComplVec);
    
    HDSDP_FREE(HLp->dColValDirection);
    HDSDP_FREE(HLp->dColDualDirection);
    HDSDP_FREE(HLp->dRowDualDirection);
    
    HDSDP_FREE(HLp->dAuxiColVector1);
    HDSDP_FREE(HLp->dAuxiColVector2);
    HDSDP_FREE(HLp->dAuxiColVector3);
    HDSDP_FREE(HLp->dAuxiColVector4);
    
    HDSDP_FREE(HLp->dAuxiRowVector1);
    HDSDP_FREE(HLp->dAuxiRowVector2);
    
    HDSDP_FREE(HLp->dKrylovAuxVec1);
    HDSDP_FREE(HLp->dKrylovAuxVec2);
    HDSDP_FREE(HLp->dKrylovAuxVec3);
    HDSDP_FREE(HLp->dKrylovAuxVec4);
    HDSDP_FREE(HLp->dKrylovAuxVec5);
    HDSDP_FREE(HLp->dKrylovAuxVec6);
    
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
