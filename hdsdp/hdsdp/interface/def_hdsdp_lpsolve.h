#ifndef def_hdsdp_lpsolve_h
#define def_hdsdp_lpsolve_h

#ifdef HEADERPATH
#include "interface/hdsdp_utils.h"
#include "linalg/hdsdp_linsolver.h"
#else
#include "hdsdp_utils.h"
#include "hdsdp_linsolver.h"
#endif

typedef struct {
    
    int nThreads;
    
    int nCol;
    int nRow;
    
    int *colMatBeg;
    int *colMatIdx;
    double *colMatElem;
    
    int *AugBeg;
    int *AugIdx;
    double *AugElem;
    
    double *colBackup;
    
    hdsdp_linsys *kkt; ///< Indefinite augmented system
    
    /* Algorithm parameters */
    double alpha;
    double beta;
    double gamma;
    double mu;
    
    /* Intermediate arrays */
    double *dd; ///< Size n
    double *xse; ///< Size n
    double *d1; ///< Size m + n
    double *d2; ///< Size m + n
    double *daux; ///< Size m + n
    
    /* Consecutive memory for [dx; dy; ds] */
    double *dx;
    double *dy;
    double *ds;
    double *dxcorr;
    double *dycorr;
    double *dscorr;
    
    double dkappa;
    double dkappacorr;
    double dtau;
    double dtaucorr;
    
    /* Primal-dual regularization of the augmented system */
    double pReg;
    double dReg;
    
    /* Signal for ill-conditioning Newton */
    int badNewton;
    
    /* Type of corrector used.
     0: no corrector
     1: Mehrotra's corrector
     2: Multiple centrality corrector
     */
    int iterType;

} hdsdp_lpsolver;

#endif /* def_hdsdp_lpsolve_h */
