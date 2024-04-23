#ifndef lp_newton_h
#define lp_newton_h

#ifdef HEADERPATH
#include "interface/def_hdsdp_lpsolve.h"
#else
#include "def_hdsdp_lpsolve.h"
#endif

extern hdsdp_retcode LpNewtonCreate( hdsdp_lpsolver **pnewton, int nThreads );
extern hdsdp_retcode LpNewtonInit( hdsdp_lpsolver *newton, int nCol, int nRow, int *colMatBeg, int *colMatIdx, double *colMatElem );
extern hdsdp_retcode LpNewtonOneStep( hdsdp_lpsolver *newton, double *lpObj, double *lpRHS, int *colMatBeg, int *colMatIdx, double *colMatElem,
                                double *colVal, double *rowDual, double *colDual, double *kappa, double *tau,
                                double *pRes, double *dRes, double pObjVal, double dObjVal, double spxSize );
extern hdsdp_retcode LpNewtonInitRobust( hdsdp_lpsolver *newton, int nCol, int nRow, int *colMatBeg, int *colMatIdx, double *colMatElem );
extern hdsdp_retcode LpNewtonOneStepRobust( hdsdp_lpsolver *newton, double *lpObj, double *lpRHS, int *colMatBeg, int *colMatIdx, double *colMatElem,
                                double *colVal, double *rowDual, double *colDual, double *kappa, double *tau,
                                double *pRes, double *dRes, double pObjVal, double dObjVal, double spxSize );
extern void LpNewtonClear( hdsdp_lpsolver *newton );
extern void LpNewtonDestroy( hdsdp_lpsolver **pnewton );


#endif /* lp_newton_h */
