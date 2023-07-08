#ifndef hdsdp_lanczos_h
#define hdsdp_lanczos_h

#include "linalg/def_hdsdp_lanczos.h"

#ifdef __cplusplus
extern "C" {
#endif

extern hdsdp_retcode HLanczosCreate( hdsdp_lanczos **pHLanczos );
extern hdsdp_retcode HLanczosInit( hdsdp_lanczos *HLanczos, int nCol, int nSpaceDim );
extern void HLanczosSetData( hdsdp_lanczos *HLanczos, void *MMat, void (*Mvec) (void *, double *, double *) );
extern hdsdp_retcode HLanczosSolve( hdsdp_lanczos *HLanczos, double *dMinEVal );
extern void HLanczosClear( hdsdp_lanczos *HLanczos );
extern void HLanczosDestroy( hdsdp_lanczos **pHLanczos );

#ifdef __cplusplus
}
#endif

#endif /* hdsdp_lanczos_h */
