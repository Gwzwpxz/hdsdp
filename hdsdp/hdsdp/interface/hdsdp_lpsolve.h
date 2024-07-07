#ifndef lp_newton_h
#define lp_newton_h

#ifdef HEADERPATH
#include "interface/def_hdsdp_lpsolve.h"
#else
#include "def_hdsdp_lpsolve.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

extern hdsdp_retcode HLpSolverCreate( hdsdp_lpsolver **pHLp );
extern hdsdp_retcode HLpSolverInit( hdsdp_lpsolver *HLp, int nRow, int nCol );
extern hdsdp_retcode HLpSolverSetData( hdsdp_lpsolver *HLp, int *colMatBeg, int *colMatIdx, double *colMatElem,
                                      int *colMatTransBeg, int *colMatTransIdx, double *colMatTransElem, double *rowRHS, double *colObj );
extern hdsdp_retcode HLpSolverOptimize( hdsdp_lpsolver *HLp );
extern void HLpSolverClear( hdsdp_lpsolver *HLp );
extern void HLpSolverDestroy( hdsdp_lpsolver **pHLp );

#ifdef __cplusplus
}
#endif

#endif /* lp_newton_h */
