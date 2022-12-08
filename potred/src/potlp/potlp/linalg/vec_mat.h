#ifndef vec_mat_h
#define vec_mat_h

#include "pot_def.h"

static double potDblConstantOne = 1.0;
static double potDblConstantMinusOne = -1.0;
static double potDblConstantZero = 0.0;
static pot_int potIntConstantOne = 1;
static pot_int potIntConstantZero = 0;
static char potCharConstantLow = 'L';
static char potCharConstantTrans = 'T';
static char potCharConstantNoTrans = 'N';

extern double nrm2( pot_int *n, double *x, pot_int *incx );
extern void axpy( pot_int *n, double *a, double *x, pot_int *incx, double *y, pot_int *incy );
extern double dot( pot_int *n, double *x, pot_int *incx, double *y, pot_int *incy );
extern void scal( pot_int *n, double *sa, double *sx, pot_int *incx );
extern void rscl( pot_int *n, double *sa, double *sx, pot_int *incx );
extern pot_int idamax( pot_int *n, double *x, pot_int *incx );
extern pot_int idamin( pot_int *n, double *x, pot_int *incx );

extern pot_int psyev( pot_int n, double *U, double *d, double *Y,
                     double *work, pot_int *iwork, pot_int lwork, pot_int liwork );
extern void pgemv( pot_int m, pot_int n, double *M, double *v, double *y );
extern void symv( char *uplo, pot_int *n, double *alpha, double *a, pot_int *lda,
                  double *x, pot_int *incx, double *beta, double *y, pot_int *incy );
extern void gemv( char *trans, pot_int *m, pot_int *n, double *alpha,
                  double *a, pot_int *lda, double *x, pot_int *incx,
                  double *beta, double *y, pot_int *incy );
extern void syrk( char *uplo, char *trans, pot_int *n, pot_int *k, double *alpha,
                  double *a, pot_int *lda, double *beta, double *c, pot_int *ldc );
extern void potrf( char *uplo, pot_int *n, double *a, pot_int *lda, pot_int *info );
extern void potrs( char *uplo, pot_int *n, pot_int *nrhs, double *a, pot_int *lda,
                  double *b, pot_int *ldb, pot_int *info );

extern double quadform( pot_int *n, double *Q, double *x );
extern double sumlogdet( pot_int *n, double *x );
extern void vvscl( pot_int *n, double *s, double *x );
extern void vvrscl( pot_int *n, double *s, double *x );

extern void spMatAxpy( pot_int n, pot_int *Ap, pot_int *Ai, double *Ax, double a, double *x, double *y );
extern void spMatATxpy( pot_int n, pot_int *Ap, pot_int *Ai, double *Ax, double a, double *x, double *y );
extern void spMatMaxRowAbs( pot_int n, pot_int *Ap, pot_int *Ai, double *Ax, double *row );
extern void spMatMaxColAbs( pot_int n, pot_int *Ap, pot_int *Ai, double *Ax, double *col );
extern void spMatRowScal( pot_int n, pot_int *Ap, pot_int *Ai, double *Ax, double *row );
extern void spMatColScal( pot_int n, pot_int *Ap, pot_int *Ai, double *Ax, double *col );
extern int spMatBuildQMat( pot_int qm, pot_int qn, pot_int *Qp, pot_int *Qi, double *Qx,
                           pot_int am, pot_int an, pot_int *Ap, pot_int *Ai, double *Ax,
                           double *b, double *c );
extern int spMatRuizScal( pot_int m, pot_int n, pot_int *Ap, pot_int *Ai, double *Ax, double *D, double *E, int maxIter );
extern int spMatPCScal( pot_int m, pot_int n, pot_int *Ap, pot_int *Ai, double *Ax, double *D, double *E, int maxIter );
extern int spMatL2Scal( pot_int m, pot_int n, pot_int *Ap, pot_int *Ai, double *Ax, double *D, double *E );

#endif /* vec_mat_h */
