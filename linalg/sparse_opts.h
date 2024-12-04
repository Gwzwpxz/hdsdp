#ifndef sparse_opts_h
#define sparse_opts_h

#ifdef HEADERPATH
#include "interface/hdsdp.h"
#else
#include "hdsdp.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

extern void csp_Axpy( int n, int *Ap, int *Ai, double *Ax, double a, double *x, double *y );
extern void csp_ATxpy( int n, int *Ap, int *Ai, double *Ax, double a, double *x, double *y );
extern double csp_sum_abs( int n, int *Ap, int *Ai, double *Ax );
extern double csp_fro_norm( int n, int *Ap, int *Ai, double *Ax );
extern void csp_aApB( int n, int nnz, double a, int *Al, double *Ax, double *Bx );
extern int csp_nnz_cols( int n, int *Ap );
extern void csp_trimultiply( int n, int *Sp, int *Si, double *Sx, double *X, double *aux, double *XSX );
extern double csp_dot_fds( int n, int *Ap, int *Ai, double *Ax, double *B );
extern void csp_max_rowabs( int n, int *Ap, int *Ai, double *Ax, double *row );
extern void csp_min_rownzabs( int m, int n, int *Ap, int *Ai, double *Ax, double *row );
extern void csp_max_colabs( int n, int *Ap, int *Ai, double *Ax, double *col );
extern void csp_min_colnzabs( int n, int *Ap, int *Ai, double *Ax, double *col );
extern void csp_rowscal( int n, int *Ap, int *Ai, double *Ax, double *row );
extern void csp_colscal( int n, int *Ap, int *Ai, double *Ax, double *col );
extern hdsdp_retcode csp_geoscal( int m, int n, int *Ap, int *Ai, double *Ax, double *D, double *E );
extern hdsdp_retcode csp_ruizscal( int m, int n, int *Ap, int *Ai, double *Ax, double *D, double *E, int maxIter );
extern void csp_l2scal( int m, int n, int *Ap, int *Ai, double *Ax, double *D, double *E );
extern void csp_dump( int n, int *Ap, int *Ai, double *Ax, double *v );

extern void tsp_decompress( int n, int nnz, int *Ci, double *Cx, int *Ai, int *Aj, double *Ax );
extern int tsp_r1_extract( int n, int nnz, int *Ai, int *Aj, double *Ax, double *sgn, double *a );
extern void tsp_scal( double a, int nnz, double *Ax );
extern double tsp_sum_abs( int nnz, int *Ai, int *Aj, double *Ax );
extern double tsp_fro_norm( int nnz, int *Ai, int *Aj, double *Ax );
extern void tsp_dump( int n, int nnz, int *Ai, int *Aj, double *Ax, double *v );
extern double tsp_quadform( int n, int nnz, int *Ai, int *Aj, double *Ax, double *v );

#ifdef __cplusplus
}
#endif

#endif /* sparse_opts_h */
