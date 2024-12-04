#ifndef hdsdp_lpkkt_h
#define hdsdp_lpkkt_h

#ifdef HEADERPATH
#include "interface/def_hdsdp_lpkkt.h"
#include "interface/hdsdp_lpsolve.h"
#else
#include "def_hdsdp_lpkkt.h"
#include "hdsdp_lpsolve.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

extern hdsdp_retcode HLpKKTCreate( hdsdp_lp_kkt **pkkt );
extern hdsdp_retcode HLpKKTInit( hdsdp_lp_kkt *kkt, int nRow, int nCol, int *colMatBeg, int *colMatIdx, double *colMatElem,
                                int *colMatTransBeg, int *colMatTransIdx, double *colMatTransElem );
extern hdsdp_retcode HLpKKTSetup( hdsdp_lp_kkt *kkt, int LpMethod, double *dScalingMatrix, double dPrimalReg, double dDualReg );
extern hdsdp_retcode HLpKKTSolveAugmented( hdsdp_lp_kkt *kkt, double *dLhsVec, double *dRhsVec );
extern hdsdp_retcode HLpKKTSolveNormalEqn( hdsdp_lp_kkt *kkt, int nRhs, double *dLhsVec, double *dRhsVec );
extern void HLpKKTMultiply( hdsdp_lp_kkt *kkt, double dXCoeff, double *dXVec, double *dYVec );
extern double HLpKKTGetFactorSolveTimeRatio( hdsdp_lp_kkt *kkt );
extern void HLpKKTClear( hdsdp_lp_kkt *kkt );
extern void HLpKKTDestroy( hdsdp_lp_kkt **pkkt );

#ifdef __cplusplus
}
#endif

#endif /* hdsdp_lpkkt_h */
