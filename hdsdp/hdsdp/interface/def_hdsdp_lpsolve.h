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
    
    hdsdp_lpsolver_params *params;
    
    int nCol;
    double *dPrimalIterHistory;
    double *dPrimalMuHistory;
    double dCondNumberEst;
    double dIterDiffMetric;
    double dIterDiffMetricScal;
    double dIterDiffMetricThresh;
    
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
    
    hdsdp_lpsolver_params params;
    hdsdp_primal_stats *stats;
    
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
    
    /* XSe */
    double *dComplVec;
    
    /* Scaling matrix */
    double *dScalingMatrix;
    
    hdsdp_lp_kkt *Hkkt;
    
    double pStep;
    double dStep;
    
    double pObjVal;
    double dObjVal;
    double dPrimalDualGap;
    
    double dBarrierMu;
    
    double *dColValDirection;
    double *dRowDualDirection;
    double *dColDualDirection;
    
    double *dAuxiColVector1;
    double *dAuxiColVector2;
    double *dAuxiColVector3;
    double *dAuxiRowVector1;
    double *dAuxiRowVector2;
    
    /* Self-dual embedding */
    double dKappa;
    double dTau;
    
    double dKappaDirection;
    double dTauDirection;

} hdsdp_lpsolver;

#endif /* def_hdsdp_lpsolve_h */
