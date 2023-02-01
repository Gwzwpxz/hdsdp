#ifndef hdsdp_conic_sdp_h
#define hdsdp_conic_sdp_h

#include "interface/def_hdsdp_conic.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Dense SDP cone */
extern hdsdp_retcode sdpDenseConeCreateImpl( hdsdp_cone_sdp_dense **pCone );
extern hdsdp_retcode sdpDenseConeProcDataImpl( hdsdp_cone_sdp_dense *cone, int nRow, int nCol,
                                               int *coneMatBeg, int *coneMatIdx, double *coneMatElem );
extern void sdpDenseConeSetStartImpl( hdsdp_cone_sdp_dense *cone, double rResi );
extern void sdpDenseConeUpdateImpl( hdsdp_cone_sdp_dense *cone, double barHsdTau, double *rowDual );
extern double sdpDenseConeRatioTestImpl( hdsdp_cone_sdp_dense *cone, double barHsdTauStep, double *rowDualStep );
extern int64_t sdpDenseConeGetSymNnzImpl( hdsdp_cone_sdp_dense *cone );
extern void sdpDenseConeAddSymNnzImpl( hdsdp_cone_sdp_dense *cone, int *schurMatCol );
extern void sdpDenseConeClearImpl( hdsdp_cone_sdp_dense *cone );
extern void sdpDenseConeDestroyImpl( hdsdp_cone_sdp_dense **pCone );
extern void sdpDenseConeViewImpl( hdsdp_cone_sdp_dense *cone );

/* Sparse SDP cone */
extern hdsdp_retcode sdpSparseConeCreateImpl( hdsdp_cone_sdp_sparse **pCone );
extern hdsdp_retcode sdpSparseConeProcDataImpl( hdsdp_cone_sdp_sparse *cone, int nRow, int nCol,
                                               int *coneMatBeg, int *coneMatIdx, double *coneMatElem );
extern void sdpSparseConeSetStartImpl( hdsdp_cone_sdp_sparse *cone, double rResi );
extern void sdpSparseConeUpdateImpl( hdsdp_cone_sdp_sparse *cone, double barHsdTau, double *rowDual );
extern double sdpSparseConeRatioTestImpl( hdsdp_cone_sdp_sparse *cone, double barHsdTauStep, double *rowDualStep );
extern int64_t sdpSparseConeGetSymNnzImpl( hdsdp_cone_sdp_sparse *cone );
extern void sdpSparseConeAddSymNnzImpl( hdsdp_cone_sdp_sparse *cone, int *schurMatCol );
extern void sdpSparseConeClearImpl( hdsdp_cone_sdp_sparse *cone );
extern void sdpSparseConeDestroyImpl( hdsdp_cone_sdp_sparse **pCone );
extern void sdpSparseConeViewImpl( hdsdp_cone_sdp_sparse *cone );

#ifdef __cplusplus
}
#endif

#endif /* hdsdp_conic_sdp_h */
