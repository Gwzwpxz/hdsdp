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
    params.dScalingThreshTol = 1e-04;
    params.dBarrierLowerBnd = 1e-16;
    params.dTimeLimit = 3600.0;
    
    params.nThreads = 8;
    params.nScalIter = 1;
    params.nMaxIter = 1000;
    
    params.iStartMethod = MEHROTRA_START;
    params.iScalMethod = SCAL_GEOMETRIC;
    params.iPrimalMethod = 1;
    
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
    
    if ( HLp->lpstats.dObjInfNorm < 1e-08 ) {
        HLp->lpstats.isNoObj = 1;
    }
    
    if ( HLp->lpstats.dRhsInfNorm < 1e-08 ) {
        HLp->lpstats.isNoRhs = 1;
    }
    
    /* Adjust solver parameters */
    if ( HLp->lpstats.isNoObj ) {
        HLp->params.dScalingThreshTol = 1e-03 / sqrt(HLp->nCol);
    }
    
    return;
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
    
    /* Get complementarity XSe */
    for ( int iCol = 0; iCol < HLp->nCol; ++iCol ) {
        HLp->dComplVec[iCol] = HLp->dColVal[iCol] * HLp->dColDual[iCol];
    }

    return;
}

static int HLpSolverIComputeSolutionStats( hdsdp_lpsolver *HLp, int iIter ) {
    
    int goOn = 1;
    
    if ( iIter == 0 ) {
        hdsdp_printf("Using Hybrid Primal-Primal-Dual solver \n\n");
        hdsdp_printf("    %5s  %15s  %15s  %8s  %8s  %10s   %5s \n",
                     "nIter", "pObj", "dObj", "pInf", "Mu", "P/D Step", "T [PD]");
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
    HLp->dBarrierMu = dot(&HLp->nCol, HLp->dColVal, &HIntConstantOne, HLp->dColDual, &HIntConstantOne);
    HLp->dBarrierMu = HLp->dBarrierMu / HLp->nCol;
    
    if ( (HLp->dPrimalDualGapRel <= HLp->params.dRelOptTol)  &&
         (HLp->dPrimalInfeasRel  <= HLp->params.dRelFeasTol) &&
         (HLp->dDualInfeasRel    <= HLp->params.dRelFeasTol) &&
         (HLp->dPrimalDualGap    <= HLp->params.dAbsOptTol)  &&
         (HLp->dPrimalInfeas     <= HLp->params.dAbsFeasTol) &&
         (HLp->dDualInfeas       <= HLp->params.dAbsFeasTol) ) {
        
        goOn = 0;
    }
    
    /* Logging */
    hdsdp_printf("    %5d  %+15.8e  %+15.8e  %8.2e  %8.2e  %5.2f %5.2f  %4.1f \n",
                 iIter, HLp->pObjVal, HLp->dObjVal, HLp->dPrimalInfeasRel, HLp->dDualInfeasRel,
           HLp->pStep, HLp->dStep, HUtilGetTimeStamp() - HLp->dTStart);
    
    return goOn;
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

static int HLpSolverICollectIterationStats( hdsdp_lpsolver *HLp ) {
    
    
    
    return 1;
}

/* Primal interior point method main algorithm */
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
    HLp->pStep = HDSDP_MIN(0.995 * HLp->pStep, 1.0);
    HLp->dStep = HDSDP_MIN(0.995 * HLp->dStep, 1.0);
    
    /* Take step into next iteration */
    axpy(&HLp->nCol, &HLp->pStep, HLp->dColValDirection, &HIntConstantOne, HLp->dColVal, &HIntConstantOne);
    axpy(&HLp->nCol, &HLp->dStep, HLp->dColDualDirection, &HIntConstantOne, HLp->dColDual, &HIntConstantOne);
    axpy(&HLp->nRow, &HLp->dStep, HLp->dRowDualDirection, &HIntConstantOne, HLp->dRowDual, &HIntConstantOne);
    
exit_cleanup:
    return retcode;
}

/* Primal-dual interior point method main algorithm. Use Mehrotra's predictor-corrector step */
static hdsdp_retcode HLpSolverITakePrimalStep( hdsdp_lpsolver *HLp ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    
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
    
    hdsdp_printf("LP Solution statistic \n");
    hdsdp_printf("pObj: %+8.3e   dObj: %+8.3e \n", HLp->pObjVal, HLp->dObjVal);
    hdsdp_printf("Abs. pInf:  %8.3e    Rel. pInf:  %8.3e \n", HLp->dPrimalInfeas, HLp->dPrimalInfeasRel);
    hdsdp_printf("Abs. dInf:  %8.3e    Rel. dInf:  %8.3e \n", HLp->dDualInfeas, HLp->dDualInfeasRel);
    hdsdp_printf("Abs. Gap :  %8.3e    Rel. gap :  %8.3e \n", HLp->dPrimalDualGap, HLp->dPrimalDualGapRel);
    hdsdp_printf("Elapsed Time:  %5.3f seconds \n\n", HUtilGetTimeStamp() - HLp->dTStart);
    
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
    HDSDP_INIT(HLp->dAuxiColVector4, double, nCol);
    
    HDSDP_MEMCHECK(HLp->dAuxiColVector1);
    HDSDP_MEMCHECK(HLp->dAuxiColVector2);
    HDSDP_MEMCHECK(HLp->dAuxiColVector3);
    HDSDP_MEMCHECK(HLp->dAuxiColVector4);
    
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
        iPrimalStart = HLpSolverICollectIterationStats(HLp);
        
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
    HDSDP_FREE(HLp->dComplVec);
    HDSDP_FREE(HLp->dScalingMatrix);
    
    HDSDP_FREE(HLp->dColValDirection);
    HDSDP_FREE(HLp->dColDualDirection);
    HDSDP_FREE(HLp->dRowDualDirection);
    
    HDSDP_FREE(HLp->dAuxiColVector1);
    HDSDP_FREE(HLp->dAuxiColVector2);
    HDSDP_FREE(HLp->dAuxiColVector3);
    HDSDP_FREE(HLp->dAuxiColVector4);
    
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
