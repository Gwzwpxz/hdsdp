/** @file def\_hdsdp\_linsolver.h
 *  @brief HDSDP linear system solver
 */
#ifndef def_hdsdp_linsolver_h
#define def_hdsdp_linsolver_h

#ifdef HEADERPATH
#include "interface/hdsdp.h"
#include "external/qdldl.h"
#include "external/lapack_names.h"
#else
#include "hdsdp.h"
#include "qdldl.h"
#include "lapack_names.h"
#endif

/* In HDSDP, there are two cases where positve definite matrices need to
   be factorized: when inverting the dual matrix and when solving the Schur complement system.
   According to the sparsity of the matrix, we choose differnt data structures to do the job
 
   Dual matrix: sparse/dense Cholesky
   Schur matrix: sparse/dense Cholesky + Pre-conditioned conjugate gradient/residual
 
 */

typedef enum {
    /* Direct solver supports both Schur complement and matrix variable */
    HDSDP_LINSYS_DENSE_DIRECT,
    HDSDP_LINSYS_SMALL_DIRECT, /* Tailored for extremely small cones of size <= 3 */
    HDSDP_LINSYS_SPARSE_DIRECT,
    HDSDP_LINSYS_SPARSE_INDEFINITE,
    HDSDP_LINSYS_SPARSE_ITERATIVE,
    
    /* Iterative solver is only used for Schur complement */
    HDSDP_LINSYS_DENSE_ITERATIVE,
    HDSDP_LINSYS_DENSE_INDEFINITE
    
} linsys_type;

typedef struct {
    
    int nCol;
    void *chol;
    linsys_type LinType;
    
    hdsdp_retcode (*cholCreate) ( void **, int );
    void (*cholSetParam) ( void *, void * );
    hdsdp_retcode (*cholSymbolic) ( void *, int *, int * );
    hdsdp_retcode (*cholNumeric) ( void *, int *, int *, double * );
    hdsdp_retcode (*cholPsdCheck) ( void *, int *, int *, double *, int * );
    
    void (*cholFSolve) ( void *, int, double *, double * );
    void (*cholBSolve) ( void *, int, double *, double * );
    hdsdp_retcode (*cholSolve) ( void *, int, double *, double * );
    hdsdp_retcode (*cholGetDiag) ( void *, double * );
    void (*cholInvert) ( void *, double *, double * );
    
    void (*cholDestroy) ( void ** );
    
    int nSolves;
    int nFactorizes;
    
} hdsdp_linsys_fp;

/* Sparse direct */
typedef struct {
    
    int nCol;
    int mType;
    int *colMatBeg;
    int *colMatIdx;
    double *colMatElem;
    
    double *dWork;
    
    void *pt[64];
    int iparm[64];
    
} pardiso_linsys;


/* Both QDLDL and LDL accept upper-triangular part of the matrix as input.
   To directly use these packages, HDSDP constructs an internal map that converts
   lower-triangular components to upper-triangular */
typedef struct {
    
    int nCol;
    int *colMatUpperBeg;
    int *colMatUpperIdx;
    double *colMatUpperElem;
    int *iTransMap; /* Map of entries from lower-triangular to upper-triangular */
    
    int *iWork;
    
    double *dWork;
    int *P;
    
    int *Lnz;
    int *Lp;
    int *Li;
    double *Lx;

} qdldl_linsys;

typedef struct {
    
    int nCol;
    int *colMatUpperBeg;
    int *colMatUpperIdx;
    double *colMatUpperElem;
    int *iTransMap; /* Map of entries from lower-triangular to upper-triangular */
    
    int *Parent;
    int *Flag;
    double *D;
    double *Y;
    int *Pattern;
    
    int *Lnz;
    int *Lp;
    int *Li;
    double *Lx;
    
} ldl_linsys;

/* Dense direct */
typedef struct {
    
    int nCol;
    double *dFullMatElem;
    
} lapack_flinsys;

typedef struct {
    
    int nCol;
    double *dFullMatElem;
    int *iAuxiIPIV;
    int iLwork;
    double *dWork;
    
} lapack_indef_flinsys;

typedef struct {
    
    int nCol;
    double dSmallMatElem[9];
    
} small_linsys;

typedef struct {
    
    int useJacobi;
    int maxIter;
    int nRestartFreq;
    double absTol;
    double relTol;
    
} iterative_params;

typedef enum {
    
    ITERATIVE_STATUS_OK,
    ITERATIVE_STATUS_NUMERICAL,
    ITERATIVE_STATUS_MAXITER,
    ITERATIVE_STATUS_FAILED
    
} iter_status;

/* Dense iterative */
typedef struct {
    
    int nCol;
    
    double *fullMatElem;
    double *iterResi;
    double *iterResiNew;
    double *iterDirection;
    double *preInvResi;
    double *MTimesDirection;
    double *iterVec;
    double *rhsBuffer;
    
    /* Pre-conditioner */
    int useJacobi;
    double *JacobiPrecond;
    lapack_flinsys *lap;
    
    /* Statistics */
    double iterResiNorm;
    double solveTime;
    
    int nIters;
    iter_status solStatus;
    int nSolves;
    
    iterative_params params;
    
} iterative_linsys;

/* Sparse direct */
#define PARDISO_RET_OK          ( 0)
#define PARDISO_RET_INDEFINITE  (-4)
#define PARDISO_SYM_POSDEFINITE ( 2)
#define PARDISO_SYM_INDEFINITE  (-2)
#define PARDISO_PHASE_SYM       (11)      // Pardiso symbolic analysis
#define PARDISO_PHASE_SYM_FAC   (12)      // Pardiso symbolic analysis
#define PARDISO_PHASE_FAC       (22)      // Pardiso numerical factorization
#define PARDISO_PHASE_FORWARD  (331)
#define PARDISO_PHASE_BACKWARD (333)
#define PARDISO_PHASE_SOLVE     (33)      // Solve linear system
#define PARDISO_PHASE_FREE      (-1)      // Free internal data structure

#define PARDISO_PARAM_NONDEFAULT    (0)
#define PARDISO_PARAM_SYMBOLIC      (1)
#define PARDISO_PARAM_SYMBOLIC_MMD  (0)
#define PARDISO_PARAM_SYMBOLIC_ND   (2)
#define PARDISO_PARAM_REFINEMENT    (7)
#define PARDISO_PARAM_INPLACE       (5)
#define PARDISO_PARAM_PERTURBATION  (9)
#define PARDISO_PARAM_SCALING      (10)
#define PARDISO_PARAM_MATCHING     (12)
#define PARDISO_PARAM_FACNNZ       (17)
#define PARDISO_PARAM_FACFLOP      (18)
#define PARDISO_PARAM_THREADS      (33)
#define PARDISO_PARAM_INDEX        (34)
#define PARDISO_PARAM_INDEX_C       (1)
#define PARDISO_PARAM_DIAGONAL     (55)
#define PARDISO_PARAM_DIAGONAL_ON   (1)

#define set_pardiso_param(iparm, param, val) iparm[param] = val
#define get_pardiso_output(iparm, param) iparm[param]

#ifdef LINSYS_PARDISO
extern void pardisoinit ( void *, int *, int * );
extern void pardiso     ( void     *, int    *, int *, int *, int *, int *,
                          double   *, int    *, int *, int *, int *, int *,
                          int *, double      *, double   *, int * );
extern void pardiso_getdiag ( const void *, void *, void *, const int *, int * );
#endif

/* Dense direct */
#define LAPACK_RET_OK    ( 0 )
#define LAPACK_UPLOW_LOW ('L')
#define LAPACK_NOTRANS   ('N')
#define LAPACK_TRANS     ('T')
#define LAPACK_SIDE_LEFT ('L')
#define LAPACK_DIAG_NONUNIT ('N')

/* For different LAPACK routines */
extern void dtrsv( const char *uplo, const char *trans, const char *diag,
                   const int *n, const double *a, const int *lda, double *x,
                   const int *incx );
void dtrsm( const char *side, const char *uplo, const char *transa,
            const char *diag, const int *m, const int *n, const double *alpha,
            const double *a, const int *lda, double *b, const int *ldb );
void dpotrf ( const char *uplo, const int *n, double *a, const int *lda, int *info );
void dsytrf( const char *uplo, const int  *n, double *a, const int  *lda,
              int *ipiv, double *work, const int *lwork, int *info );
void dpotri( const char *uplo, const int *n, double *a, const int *lda, int *info );
void dpotrs( const char *uplo, const int *n, const int *nrhs, const double *a,
             const int *lda, double *b, const int *ldb, int *info );
void dsytrs( const char *uplo, const int *n, const int *nrhs, const double *a,
               const int *lda, const int *ipiv, double *b, const int *ldb, int *info );

#endif /* def_hdsdp_linsolver_h */
