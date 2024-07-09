#ifndef def_hdsdp_lpsolve_h
#define def_hdsdp_lpsolve_h

#ifdef HEADERPATH
#include "interface/hdsdp_utils.h"
#include "linalg/hdsdp_linsolver.h"
#include "interface/def_hdsdp_lpkkt.h"
#else
#include "hdsdp_utils.h"
#include "hdsdp_linsolver.h"
#include "def_hdsdp_lpkkt.h"
#endif


#define LP_ITER_MEHROTRA (0)
#define LP_ITER_PRIMAL   (1)
#define LP_ITER_HSD      (2)
    
#define SCAL_NONE      (-1)
#define SCAL_RUIZ      (0)
#define SCAL_GEOMETRIC (1)
#define SCAL_L2NORM    (2)

#define MEHROTRA_START (0)
#define TRIVIAL_START  (1)

typedef struct {
    
    /* Tolerance of LP solver */
    double dAbsOptTol;
    double dAbsFeasTol;
    double dRelOptTol;
    double dRelFeasTol;
    
    /* KKT regularization */
    double dKKTPrimalReg;
    double dKKTDualReg;
    
    /* Potential parameter */
    double dPotentialRho;
    
    /* Primal and dual update */
    double dPrimalUpdateStep;
    double dDualUpdateStep;
    
    /* Tolerance of iterative solver */
    double dIterativeTol;
    /* Threshold of scaling matrix check */
    double dScalingThreshTol;
    /* Lower bound of barrier parameter*/
    double dBarrierLowerBnd;
    
    /* Solution time */
    double dTimeLimit;
    
    /* Starting point */
    int iStartMethod;
    
    /* Threads for linear system */
    int nThreads;
    
    /* Maximum iteration */
    int nMaxIter;
    
    /* Scaling */
    int iScalMethod;
    int nScalIter;
    
    /* Whether primal method is on */
    int iPrimalMethod;
    
    /* LP method */
    int LpMethod;
    
} hdsdp_lpsolver_params;

typedef struct {
    
    double dRhsOneNorm;
    double dRhsInfNorm;
    double dRhsTwoNorm;
    double dObjOneNorm;
    double dObjTwoNorm;
    double dObjInfNorm;
    double dAMatFNorm;
    double dAMatAbsNorm;
    
    int nAMatNz;
    
    int isNoObj;
    int isNoRhs;
    
} hdsdp_lpsolver_stats;

typedef struct {
    
    hdsdp_lpsolver_params *params;
    
    int nCol;
    int iIterLast;
    int isSuperLin;
    
    double *dPrimalIterHistory;
    double *dPrimalMuHistory;
    double dCondNumberEst;
    double dIterDiffMetric;
    double dIterDiffMetricScal;
    double dIterDiffMetricThresh;
    double dIterDiffMetricAggressive;
    
} hdsdp_primal_stats;

typedef struct {
    
    int nRow;
    int nCol;
    int nElem;
    
    /* LP Coefficient matrix and its transpose */
    int *colMatBeg;
    int *colMatIdx;
    double *colMatElem;
    
    int *colMatTransBeg;
    int *colMatTransIdx;
    double *colMatTransElem;
    
    /* Coefficient vectors */
    double *dRowRhs;
    double *dColObj;
    
    /* Scaling vector */
    double *dColScal;
    double *dRowScal;
    
    hdsdp_lpsolver_stats lpstats;
    hdsdp_lpsolver_params params;
    hdsdp_primal_stats *pstats;
    
    /* Solution status */
    hdsdp_status LpStatus;
    /* x */
    double *dColVal;
    /* s */
    double *dColDual;
    /* y */
    double *dRowDual;
    /* A * x - b * tau */
    double *dPrimalInfeasVec;
    /* A' * y + s - c * tau */
    double *dDualInfeasVec;
    /* Scaling matrix */
    double *dScalingMatrix;
    
    hdsdp_lp_kkt *Hkkt;
    
    double pStep;
    double dStep;
    
    double pObjVal;
    double dObjVal;
    double dPrimalDualGap;
    double dPrimalDualGapRel;
    double dProxNorm;
    
    double dPrimalInfeas;
    double dPrimalInfeasRel;
    double dDualInfeas;
    double dDualInfeasRel;
    
    double dBarrierMu;
    
    double *dColValDirection;
    double *dRowDualDirection;
    double *dColDualDirection;
    
    double *dAuxiColVector1;
    double *dAuxiColVector2;
    double *dAuxiColVector3;
    double *dAuxiColVector4;
    double *dAuxiRowVector1;
    double *dAuxiRowVector2;
    
    /* Embedded Krylov subspace solver */
    double *dKrylovAuxVec1;
    double *dKrylovAuxVec2;
    double *dKrylovAuxVec3;
    double *dKrylovAuxVec4;
    double *dKrylovAuxVec5;
    double *dKrylovAuxVec6;
    
    /* Self-dual embedding */
    double dKappa;
    double dTau;
    
    double dKappaDirection;
    double dTauDirection;
    
    double dTStart;

} hdsdp_lpsolver;

#endif /* def_hdsdp_lpsolve_h */
