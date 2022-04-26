#include "schurmat.h"
#include "densemat.h"
#include "sparsemat.h"
#include "rankkmat.h"
#include "vec.h"

// Error type
static char etype[] = "Schur Matrix operation Error";

static DSDP_INT schurMatselectTypeDense( schurMat *sMat ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    sMat->denseM = (dsMat *) calloc(1, sizeof(dsMat));
    denseMatInit(sMat->denseM);
    retcode = denseMatAlloc(sMat->denseM, sMat->m, TRUE);
    // Get the diagonal entries
    for (DSDP_INT i = 0; i < sMat->m; ++i) {
        sMat->diag[i] = &fullIdx(sMat->denseM->array, sMat->m, i, i);
    }
    return retcode;
}

static DSDP_INT schurMatselectTypeSparse( schurMat *sMat, DSDP_INT *nnzIdx ) {
    // Transform the index array into CSC format
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    sMat->spsM = (spsMat *) calloc(1, sizeof(spsMat));
    spsMatInit(sMat->spsM);
    // Count total number of nonzeros
    DSDP_INT nnzs = 0, m = sMat->m, i, j;
    for (i = 0; i < nsym(m); ++i) {
        if (nnzIdx[i] > 0) { nnzs += 1; }
    }
    spsMatAllocData(sMat->spsM, sMat->m, nnzs);
    DSDP_INT *Mp = sMat->spsM->p, *Mi = sMat->spsM->i;
    double *Mx = sMat->spsM->x; Mp[0] = 0;
    for (i = 0, nnzs = 0; i < m; ++i) {
        assert(nnzIdx[nnzs] > 0);
        sMat->diag[i] = &Mx[nnzs];
        for (j = 1; j < m - i; ++j) {
            if (nnzIdx[nnzs] > 0) { nnzs += 1; Mi[nnzs] = j; }
        }
        Mp[i] = nnzs;
    }
    assert(nnzs == sMat->spsM->nnz);
    // Symbolic factorization for the Schur matrix
    spsMatSymbolic(sMat->spsM);
    return retcode;
}

extern void schurMatInit( schurMat *sMat ) {
    // Initialize Schur matrix
    sMat->m = 0; sMat->stype = SCHUR_TYPE_UNKNOWN;
    sMat->spsM = NULL; sMat->denseM = NULL;
    sMat->diag = NULL; sMat->isFactorized = FALSE;
    sMat->isillCond = FALSE;
}

extern DSDP_INT schurMatAlloc( schurMat *sMat, DSDP_INT dim ) {
    // Allocate preliminary arrays
    DSDP_INT retcode = DSDP_RETCODE_OK;
    sMat->m = dim;
    sMat->diag = (double **) calloc(dim, sizeof(double *));
    return retcode;
}

extern DSDP_INT schurMatselectType( schurMat *sMat, DSDP_INT stype, DSDP_INT *nnzIdx ) {
    // Allocate data dependent type for the schur matrix
    // Also filling in the diag vector by the addresses
    DSDP_INT retcode = DSDP_RETCODE_OK;
    sMat->stype = stype;
    if (stype == SCHUR_TYPE_DENSE) {
        retcode = schurMatselectTypeDense(sMat);
    } else {
        retcode = schurMatselectTypeSparse(sMat, nnzIdx);
    }
    return retcode;
}

extern void schurMatFree( schurMat *sMat ) {
    // Free the Schur matrix
    sMat->m = 0; DSDP_FREE(sMat->diag);
    if (sMat->stype == SCHUR_TYPE_DENSE) {
        denseMatFree(sMat->denseM); DSDP_FREE(sMat->denseM);
    } else if (sMat->stype == SCHUR_TYPE_SPARSE) {
        spsMatFree(sMat->spsM); DSDP_FREE(sMat->spsM);
    } else {
        assert( FALSE );
    }
    
    DSDP_FREE(sMat->diag);
}

extern void schurMatGetdiag( schurMat *sMat, vec *diag ) {
    // Extract diagonal
    for (DSDP_INT i = 0; i < sMat->m; ++i) {
        diag->x[i] = *sMat->diag[i];
    }
}

extern void schurMatAdddiag( schurMat *sMat, double d ) {
    // Add element to the diagonal
    for (DSDP_INT i = 0; i < sMat->m; ++i) {
        *sMat->diag[i] += d;
    }
}

extern void schurMatFactorize( schurMat *sMat ) {
    // Factorize the schur matrix
    if (sMat->stype == SCHUR_TYPE_DENSE) {
        denseMatFactorize(sMat->denseM);
    } else if (sMat->stype == SCHUR_TYPE_SPARSE) {
        DSDP_INT ispsd = FALSE;
        
        if (!sMat->isillCond) {
            spsMatIspd(sMat->spsM, &ispsd);
            if (!ispsd) {
                // TODO: Add an alternative to factorize the matrix when it is not positive definite
            }
        } else {
            // TODO: Add an alternative to factorize the matrix when it is not positive definite
        }
    } else {
        assert( FALSE );
    }
    
    sMat->isFactorized = TRUE;
}

extern void schurMatSolve( schurMat *sMat, DSDP_INT nrhs, double *B, double *aux ) {
    // Solve M * x = b
    if (sMat->stype == SCHUR_TYPE_DENSE) {
        denseArrSolveInp(sMat->denseM, nrhs, B);
    } else if (sMat->stype == SCHUR_TYPE_SPARSE ){
        spsArrSolveInp(sMat->spsM, nrhs, B, aux);
    } else {
        assert( FALSE );
    }
}

extern void schurMatMx( schurMat *sMat, vec *x, vec *Ax ) {
    if (sMat->stype == SCHUR_TYPE_DENSE) {
        denseMataAxpby(sMat->denseM, 1.0, x, 0.0, Ax);
    } else if (sMat->stype == SCHUR_TYPE_SPARSE ){
        spsMatAx2(sMat->spsM, x, Ax);
    } else {
        assert( FALSE );
    }
}

extern void schurMatReset( schurMat *sMat, DSDP_INT cleanM ) {
    if (sMat->stype == SCHUR_TYPE_DENSE) {
        if (cleanM) {
            denseMatReset(sMat->denseM);
        } else {
            denseMatResetFactor(sMat->denseM);
            sMat->isFactorized = FALSE;
        }
    } else if (sMat->stype == SCHUR_TYPE_SPARSE ){
        spsMatReset(sMat->spsM);
        sMat->isFactorized = FALSE;
    } else {
        assert( FALSE );
    }
}