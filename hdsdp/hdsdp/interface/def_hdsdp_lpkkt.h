#ifndef def_hdsdp_lpkkt_h
#define def_hdsdp_lpkkt_h

#ifdef HEADERPATH
#include "interface/def_hdsdp_lpsolve.h"
#include "interface/hdsdp_utils.h"
#include "linalg/hdsdp_linsolver.h"
#else
#include "def_hdsdp_lpsolve.h"
#include "hdsdp_utils.h"
#include "hdsdp_linsolver.h"
#endif

/* KKT solver for LP solver */
typedef struct {
    
    int nLpRow;
    int nLpCol;
    int nLpElem;
    
    int nKKTAugCol;
    int nKKTAugElem;
    
    const int *colMatBeg;
    const int *colMatIdx;
    const double *colMatElem;
    
    const int *colMatTransBeg;
    const int *colMatTransIdx;
    const double *colMatTransElem;
    
    int *KKTMatBeg;
    int *KKTMatIdx;
    double *KKTMatElem;
    
    lp_method LpMethod;
    
    /* Direct factorization of augmented system */
    hdsdp_linsys *KKTDirect;
    
    /* Iterative solver for normal equation */
    hdsdp_linsys *KKTIterative;
    
    int nSolve;
    int nFactor;
    
} hdsdp_lp_kkt;


#endif /* def_hdsdp_lpkkt_h */
