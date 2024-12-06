/* ========================================================================== */
/* CXSparse/Include/cs.h file */
/* ========================================================================== */

/* This is the CXSparse/Include/cs.h file.  It has the same name (cs.h) as
   the CSparse/Include/cs.h file.  The 'make install' for SuiteSparse installs
   CXSparse, and this file, instead of CSparse.  The two packages have the same
   cs.h include filename, because CXSparse is a superset of CSparse.  Any user
   program that uses CSparse can rely on CXSparse instead, with no change to the
   user code.  The #include "cs.h" line will work for both versions, in user
   code, and the function names and user-visible typedefs from CSparse all
   appear in CXSparse.  For experimenting and changing the package itself, I
   recommend using CSparse since it's simpler and easier to modify.  For
   using the package in production codes, I recommend CXSparse since it has
   more features (support for complex matrices, and both int and long
   versions).
 */

/* ========================================================================== */

#ifndef _HDSDP_CS_
#define _HDSDP_CS_

#include "stddef.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct hdsdp_cs_sparse {
    int nzmax;
    int m;         /* number of rows */
    int n;         /* number of columns */
    int *p;        /* column pointers (size n+1) or col indices (size nzmax) */
    int *i;        /* row indices, size nzmax */
    double *x;     /* numerical values, size nzmax */
    int nz;        /* # of entries in triplet matrix, -1 for compressed-col */
} dcs ;

int dcs_entry (dcs *T, int i, int j, double x) ;
dcs *dcs_compress (const dcs *T) ;
double dcs_norm (const dcs *A) ;
int dcs_print (const dcs *A, int brief) ;

/* utilities */
void *dcs_calloc (int n, size_t size) ;
void *dcs_free (void *p) ;
void *dcs_realloc (void *p, int n, size_t size, int *ok) ;
dcs *dcs_spalloc (int m, int n, int nzmax, int values, int t) ;
dcs *dcs_spfree (dcs *A) ;
int dcs_sprealloc (dcs *A, int nzmax) ;
void *dcs_malloc (int n, size_t size) ;

/* utilities */
double dcs_cumsum (int *p, int *c, int n) ;
dcs *dcs_done (dcs *C, void *w, void *x, int ok) ;
int *dcs_idone (int *p, dcs *C, void *w, int ok) ;

/* Cholesky relevant */
typedef struct hdsdp_cs_symbolic  /* symbolic Cholesky, LU, or QR analysis */
{
    int *pinv ;     /* inverse row perm. for QR, fill red. perm for Chol */
    int *q ;        /* fill-reducing column permutation for LU and QR */
    int *parent ;   /* elimination tree for Cholesky and QR */
    int *cp ;       /* column pointers for Cholesky, row counts for QR */
    int *leftmost ; /* leftmost[i] = min(find(A(i,:))), for QR */
    int m2 ;        /* # of rows for QR, after adding fictitious rows */
    double lnz ;    /* # entries in L for LU or Cholesky; in V for QR */
    double unz ;    /* # entries in U for LU; in R for QR */
} dcss ;

typedef struct cs_numeric   /* numeric Cholesky, LU, or QR factorization */
{
    dcs *L ;         /* L for LU and Cholesky, V for QR */
    dcs *U ;         /* U for LU, R for QR, not used for Cholesky */
    int *pinv ;     /* partial pivoting for LU */
    double *B ;     /* beta [0..n-1] for QR */
} dcsn ;

dcsn *dcs_nfree (dcsn *N);
dcss *dcs_sfree (dcss *S);

int dcs_ereach (const dcs *A, int k, const int *parent, int *s, int *w);
int *dcs_etree (const dcs *A, int ata);
int dcs_fkeep (dcs *A, int (*fkeep) (int, int, double, void *), void *other);
int dcs_leaf (int i, int j, const int *first, int *maxfirst, int *prevleaf,
              int *ancestor, int *jleaf);
int dcs_lsolve (const dcs *L, double *x);
int dcs_ltsolve (const dcs *L, double *x);
dcs *dcs_multiply (const dcs *A, const dcs *B);
int *dcs_pinv (int const *p, int n);
int *dcs_pinv (int const *p, int n);
int dcs_pvec (const int *p, const double *b, double *x, int n);
int *dcs_post (const int *parent, int n);
int dcs_scatter (const dcs *A, int j, double beta, int *w, double *x, int mark,
             dcs *C, int nz);
dcs *dcs_symperm (const dcs *A, const int *pinv, int values);
int dcs_tdfs (int j, int k, int *head, const int *next, int *post, int *stack) ;
dcs *dcs_transpose (const dcs *A, int values) ;

/* Cholesky routines */
dcss *dcs_schol (int order, const dcs *A);
dcsn *dcs_chol (const dcs *A, const dcss *S);
int dcs_cholsol (int order, const dcs *A, double *b);

/* LDL routines */
void ldl_symbolic (int n, int Ap [ ], int Ai [ ], int Lp [ ],
    int Parent [ ], int Lnz [ ], int Flag [ ], int P [ ],
    int Pinv [ ]) ;

int ldl_numeric (int n, int Ap [ ], int Ai [ ], double Ax [ ],
    int Lp [ ], int Parent [ ], int Lnz [ ], int Li [ ],
    double Lx [ ], double D [ ], double Y [ ], int Pattern [ ],
    int Flag [ ], int P [ ], int Pinv [ ]) ;

void ldl_lsolve (int n, double X [ ], int Lp [ ], int Li [ ],
    double Lx [ ]) ;

void ldl_dsolve (int n, double X [ ], double D [ ]) ;

void ldl_ltsolve (int n, double X [ ], int Lp [ ], int Li [ ],
    double Lx [ ]) ;

void ldl_perm  (int n, double X [ ], double B [ ], int P [ ]) ;
void ldl_permt (int n, double X [ ], double B [ ], int P [ ]) ;

int ldl_valid_perm (int n, int P [ ], int Flag [ ]) ;
int ldl_valid_matrix ( int n, int Ap [ ], int Ai [ ]) ;

#define LDL_symbolic ldl_symbolic
#define LDL_numeric ldl_numeric
#define LDL_lsolve ldl_lsolve
#define LDL_dsolve ldl_dsolve
#define LDL_ltsolve ldl_ltsolve

#define IS_CSC(A) (A && (A->nz == -1))
#define IS_TRIPLET(A) (A && (A->nz >= 0))
#define CS_MAX(a,b) (((a) > (b)) ? (a) : (b))
#define CS_MIN(a,b) (((a) < (b)) ? (a) : (b))
#define CS_FLIP(i) (-(i)-2)
#define CS_UNFLIP(i) (((i) < 0) ? CS_FLIP(i) : (i))
#define CS_MARK(w,j) { w [j] = CS_FLIP (w [j]) ; }
#define CS_MARKED(w,j) (w [j] < 0)

#ifdef __cplusplus
}
#endif

#endif
