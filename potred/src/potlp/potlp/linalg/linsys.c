#include "linsys.h"
#include "qdldl.h"
#include "pardiso.h"

typedef struct {
    
    pot_int    n; ///< Dimension of the linear system
    
    pot_int *iWork; ///< Working array
    double *dWork; ///< Working array
    pot_int    *P; ///< Symbolic tree
    
    pot_int    *Lnz;
    pot_int    *Lp; ///< CSC representation of L
    pot_int    *Li;
    double *Lx;
    
} qdldl_linsys;

typedef struct {
    
    pot_int  n;
    
    pot_int *Ap;
    pot_int *Ai;
    double  *Ax;
    
    double *dWork;
    void   *pt[64];
    int    iparm[64];
    
} pds_linsys;

#ifdef LINSYSDUMMY
static int dummyCreate( void **pdummy, int n ) {
    
    return RETCODE_OK;
}

static int dummySymbolic( void *dummy, int *Ap, int *Ai ) {
    
    return RETCODE_OK;
}

static int dummyNumeric( void *dummy, int *Ap, int *Ai, double *Ax ) {
    
    return RETCODE_OK;
}

static int dummySolve( void *dummy, double *bx ) {
    
    return RETCODE_OK;
}

static void dummyDestroy( void **pdummy ) {
    
    if ( !pdummy ) {
        return;
    }
    
    POTLP_FREE(*pdummy);
    return;
}
#endif

#ifdef QDLDL
static int ldlCreate( void **pldl, int n ) {
    
    int retcode = RETCODE_OK;
    qdldl_linsys *qdldl = NULL;
    
    POTLP_INIT(qdldl, qdldl_linsys, 1);
    
    if ( !qdldl ) {
        retcode = RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    qdldl->n = n;
    *pldl = qdldl;

exit_cleanup:
    
    return retcode;
}

static int ldlSymbolic( void *ldl, int *Ap, int *Ai ) {
    
    int retcode = RETCODE_OK;
    qdldl_linsys *qdldl = (qdldl_linsys *) ldl;
    
    POTLP_INIT(qdldl->P, int, qdldl->n);
    POTLP_INIT(qdldl->iWork, int, qdldl->n * 4);
    POTLP_INIT(qdldl->dWork, double, qdldl->n * 3);
    POTLP_INIT(qdldl->Lp, int, qdldl->n + 1);
    POTLP_INIT(qdldl->Lnz, int, qdldl->n);
    
    if ( !qdldl->P || !qdldl->Lp || !qdldl->iWork ||
         !qdldl->dWork || !qdldl->Lnz ) {
        retcode = RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    int ldlret = QDLDL_etree(qdldl->n, Ap, Ai, qdldl->iWork, qdldl->Lnz, qdldl->P);
    
    if ( ldlret == -1 || ldlret == -2 ) {
        retcode = RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    POTLP_INIT(qdldl->Li, int, ldlret);
    POTLP_INIT(qdldl->Lx, double, ldlret);
    
    if ( !qdldl->Li || !qdldl->Lx ) {
        retcode = RETCODE_FAILED;
        goto exit_cleanup;
    }
    
exit_cleanup:
    return retcode;
}

static int ldlNumeric( void *ldl, int *Ap, int *Ai, double *Ax ) {
    
    int retcode = RETCODE_OK;
    qdldl_linsys *qdldl = (qdldl_linsys *) ldl;
    
    int n = qdldl->n;
    
    double *D = qdldl->dWork;
    double *Dinv = qdldl->dWork + n;
    double *fWork = Dinv + n;
    
    int *bWork = qdldl->iWork;
    int *iWork = bWork + n;
    
    int nPos = QDLDL_factor(n, Ap, Ai, Ax, qdldl->Lp, qdldl->Li, qdldl->Lx,
                            D, Dinv, qdldl->Lnz, qdldl->P, bWork, iWork, fWork);
    
    if ( nPos == -1 ) {
        retcode = RETCODE_FAILED;
        goto exit_cleanup;
    }
    
exit_cleanup:
    return retcode;
}

static int ldlSolve( void *ldl, double *bx ) {
    
    qdldl_linsys *qdldl = (qdldl_linsys *) ldl;
    QDLDL_solve(qdldl->n, qdldl->Lp, qdldl->Li, qdldl->Lx,
                qdldl->dWork + qdldl->n, bx);
    
    return RETCODE_OK;
}

static void ldlDestroy( void **pldl ) {
    
    if ( !pldl ) {
        return;
    }
    
    qdldl_linsys *qdldl = (qdldl_linsys *) (*pldl);
    
    if ( qdldl ) {
        
        POTLP_FREE(qdldl->P);
        POTLP_FREE(qdldl->iWork);
        POTLP_FREE(qdldl->dWork);
        POTLP_FREE(qdldl->Lnz);
        POTLP_FREE(qdldl->Lp);
        POTLP_FREE(qdldl->Li);
        POTLP_FREE(qdldl->Lx);
        
        POTLP_ZERO(qdldl, qdldl_linsys, 1);
    }
    
    POTLP_FREE(*pldl);
    
    return;
}
#else
static int pdsCreate( void **pldl, int n ) {
    
    int retcode = RETCODE_OK;
    pds_linsys *pds = NULL;
    
    POTLP_INIT(pds, pds_linsys, 1);
    
    if ( !pds ) {
        retcode = RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    pds->n = n;
    *pldl = pds;
    
    /* Initialize pardiso */
    POTLP_ZERO(pds->pt, void *, 64);
    POTLP_ZERO(pds->iparm, int, 64);
    
    int mtype = PARDISO_SYM_INDEFINITE;
    pardisoinit(pds->pt, &mtype, pds->iparm);
    
    set_pardiso_param(pds->iparm, PARDISO_PARAM_NONDEFAULT, 1);
    set_pardiso_param(pds->iparm, PARDISO_PARAM_SYMBOLIC, PARDISO_PARAM_SYMBOLIC_MMD);
    set_pardiso_param(pds->iparm, PARDISO_PARAM_PERTURBATION, 3);
    set_pardiso_param(pds->iparm, PARDISO_PARAM_INPLACE, 1);
    set_pardiso_param(pds->iparm, PARDISO_PARAM_INDEX, PARDISO_PARAM_INDEX_C);
    
exit_cleanup:
    return retcode;
}

#define SHOW_ORDERING(ord) // printf("Using ordering %s. \n", ord)
static int pdsSymbolic( void *ldl, pot_int *Ap, pot_int *Ai ) {
    
    int retcode = RETCODE_OK;
    pds_linsys *pds = (pds_linsys *) ldl;
    
    pds->Ap = Ap; pds->Ai = Ai;
    
    POTLP_INIT(pds->dWork, double, pds->n);
    
    if ( !pds->dWork ) {
        retcode = RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    pot_int maxfct = 1, mnum = 1, mtype = PARDISO_SYM_INDEFINITE, phase = PARDISO_PHASE_SYM;
    pot_int idummy = 0, msg = 0, pdsret = PARDISO_RET_OK;
    
    /* Use amd first */
    int amdFactorNnz = 0;
    set_pardiso_param(pds->iparm, PARDISO_PARAM_SYMBOLIC, PARDISO_PARAM_SYMBOLIC_MMD);
    set_pardiso_param(pds->iparm, PARDISO_PARAM_FACNNZ, -1);
    
    pardiso(pds->pt, &maxfct, &mnum, &mtype, &phase,
            &pds->n, NULL, Ap, Ai, &idummy, &idummy,
            pds->iparm, &msg, NULL, NULL, &pdsret);
    
    if ( pdsret != PARDISO_RET_OK ) {
        retcode = RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    amdFactorNnz = get_pardiso_output(pds->iparm, PARDISO_PARAM_FACNNZ);
    
    /* See if nested dissection is better */
    pot_int ndFactorNnz = 0;
    set_pardiso_param(pds->iparm, PARDISO_PARAM_SYMBOLIC, PARDISO_PARAM_SYMBOLIC_ND);
    set_pardiso_param(pds->iparm, PARDISO_PARAM_FACNNZ, -1);
    pardiso(pds->pt, &maxfct, &mnum, &mtype, &phase,
            &pds->n, NULL, Ap, Ai, &idummy, &idummy,
            pds->iparm, &msg, NULL, NULL, &pdsret);
    
    if ( pdsret != PARDISO_RET_OK ) {
        retcode = RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    ndFactorNnz = get_pardiso_output(pds->iparm, PARDISO_PARAM_FACNNZ);
    if ( ndFactorNnz > amdFactorNnz ) {
        set_pardiso_param(pds->iparm, PARDISO_PARAM_SYMBOLIC, PARDISO_PARAM_SYMBOLIC_MMD);
        SHOW_ORDERING("Minimum degree");
    } else {
        SHOW_ORDERING("Nested dissection");
    }
    
exit_cleanup:
    return retcode;
}

static int pdsScalNumeric( void *ldl, pot_int *Ap, pot_int *Ai, double *Ax ) {
    
    /* Rescue the IPM if factorization fails due to the highly indefinite system */
    int retcode = RETCODE_OK;
    pds_linsys *pds = (pds_linsys *) ldl;
    pds->Ax = Ax;
    
    pot_int maxfct = 1, mnum = 1, mtype = PARDISO_SYM_INDEFINITE, phase = PARDISO_PHASE_SYM_FAC;
    pot_int idummy = 0, msg = 0, pdsret = PARDISO_RET_OK;
    
    set_pardiso_param(pds->iparm, PARDISO_PARAM_SCALING, 1);
    set_pardiso_param(pds->iparm, PARDISO_PARAM_MATCHING, 1);
    
    pardiso(pds->pt, &maxfct, &mnum, &mtype, &phase,
            &pds->n, Ax, Ap, Ai, &idummy, &idummy,
            pds->iparm, &msg, NULL, NULL, &pdsret);
    
    if ( pdsret != PARDISO_RET_OK ) {
        retcode = RETCODE_FAILED;
        goto exit_cleanup;
    }
    
exit_cleanup:
    return retcode;
}

static int pdsNumeric( void *ldl, pot_int *Ap, pot_int *Ai, double *Ax ) {
    
    int retcode = RETCODE_OK;
    pds_linsys *pds = (pds_linsys *) ldl;
    pds->Ax = Ax;
    
    pot_int maxfct = 1, mnum = 1, mtype = PARDISO_SYM_INDEFINITE, phase = PARDISO_PHASE_FAC;
    pot_int idummy = 0, msg = 0, pdsret = PARDISO_RET_OK;
    
    pardiso(pds->pt, &maxfct, &mnum, &mtype, &phase,
            &pds->n, Ax, Ap, Ai, &idummy, &idummy,
            pds->iparm, &msg, NULL, NULL, &pdsret);
    
    if ( pdsret != PARDISO_RET_OK ) {
        retcode = RETCODE_FAILED;
        goto exit_cleanup;
    }
    
exit_cleanup:
    return retcode;
}

static int pdsSolve( void *ldl, double *bx ) {
    
    pot_int retcode = RETCODE_OK;
    pds_linsys *pds = (pds_linsys *) ldl;
    
    pot_int maxfct = 1, mnum = 1, mtype = PARDISO_SYM_INDEFINITE, phase = PARDISO_PHASE_SOLVE;
    pot_int idummy = 0, nrhs = 1, msg = 0, pdsret = PARDISO_RET_OK;
    
    pardiso(pds->pt, &maxfct, &mnum, &mtype, &phase,
            &pds->n, pds->Ax, pds->Ap, pds->Ai, &idummy, &nrhs,
            pds->iparm, &msg, bx, pds->dWork, &pdsret);
    
    if ( pdsret != PARDISO_RET_OK ) {
        retcode = RETCODE_FAILED;
        return retcode;
    }
    
    return retcode;
}

static void pdsDestroy( void **pldl ) {
    
    if ( !pldl ) {
        return;
    }
    
    pds_linsys *pds = (pds_linsys *) (*pldl);
    
    if ( pds ) {
        
        pot_int maxfct = 1, mnum = 1, mtype = PARDISO_SYM_INDEFINITE, phase = PARDISO_PHASE_FREE;
        pot_int idummy = 0, msg = 0, pdsret = PARDISO_RET_OK;
        pardiso(pds->pt, &maxfct, &mnum, &mtype, &phase,
                &pds->n, pds->Ax, pds->Ap, pds->Ai, &idummy, &idummy,
                pds->iparm, &msg, NULL, NULL, &pdsret);
        
        POTLP_FREE(pds->dWork);
        POTLP_ZERO(pds, pds_linsys, 1);
    }
    
    
    POTLP_FREE(*pldl);

    return;
}
#endif

extern pot_int potLinsysCreate( pot_linsys **ppotLinsys ) {
    
    pot_int retcode = RETCODE_OK;
    
    if ( !ppotLinsys ) {
        retcode = RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    pot_linsys *potLinsys = NULL;
    POTLP_INIT(potLinsys, pot_linsys, 1);
    
    if ( !potLinsys ) {
        retcode = RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    POTLP_ZERO(potLinsys, pot_linsys, 1);
    
    potLinsys->backUpLin = 0;
    
#ifdef LINSYSDUMMY
    potLinsys->LCreate = dummyCreate;
    potLinsys->LDestroy = dummyDestroy;
    potLinsys->LSFac = dummySymbolic;
    potLinsys->LNFac = dummyNumeric;
    potLinsys->LNFacBackup = dummyNumeric;
    potLinsys->LSolve = dummySolve;
#else
#ifdef QDLDL
    potLinsys->LCreate = ldlCreate;
    potLinsys->LDestroy = ldlDestroy;
    potLinsys->LSFac = ldlSymbolic;
    potLinsys->LNFac = ldlNumeric;
    potLinsys->LNFacBackup = ldlNumeric;
    potLinsys->LSolve = ldlSolve;
#else
    potLinsys->LCreate = pdsCreate;
    potLinsys->LDestroy = pdsDestroy;
    potLinsys->LSFac = pdsSymbolic;
    potLinsys->LNFac = pdsNumeric;
    potLinsys->LNFacBackup = pdsScalNumeric;
    potLinsys->LSolve = pdsSolve;
#endif
#endif
    
    *ppotLinsys = potLinsys;
    
exit_cleanup:
    return retcode;
}

extern pot_int potLinsysInit( pot_linsys *potLinsys, pot_int nCol ) {
    
    pot_int retcode = RETCODE_OK;
    
    potLinsys->nCol = nCol;
    retcode = potLinsys->LCreate(&potLinsys->solver, nCol);
    
    if ( retcode != RETCODE_OK ) {
        retcode = RETCODE_FAILED;
        goto exit_cleanup;
    }
    
exit_cleanup:
    return retcode;
}

extern pot_int potLinsysSymFactorize( pot_linsys *potLinsys, pot_int *colMatBeg, pot_int *colMatIdx ) {
    
    pot_int retcode = RETCODE_OK;
    retcode = potLinsys->LSFac(potLinsys->solver, colMatBeg, colMatIdx);
    
    if ( retcode != RETCODE_OK ) {
        retcode = RETCODE_FAILED;
        goto exit_cleanup;
    }
    
exit_cleanup:
    return retcode;
}

extern pot_int potLinsysNumFactorize( pot_linsys *potLinsys, pot_int *colMatBeg, pot_int *colMatIdx, double *colMatElem ) {
    
    pot_int retcode = RETCODE_OK;
    
    if ( potLinsys->backUpLin ) {
        retcode = potLinsys->LNFacBackup(potLinsys->solver, colMatBeg, colMatIdx, colMatElem);
        if ( retcode != RETCODE_OK ) {
            retcode = RETCODE_FAILED;
            goto exit_cleanup;
        }
    } else {
        retcode = potLinsys->LNFac(potLinsys->solver, colMatBeg, colMatIdx, colMatElem);
        if ( retcode != RETCODE_OK ) {
            /* Use a backup factorization */
            retcode = potLinsys->LNFacBackup(potLinsys, colMatBeg, colMatIdx, colMatElem);
            if ( retcode != RETCODE_OK ) {
                retcode = RETCODE_FAILED;
                goto exit_cleanup;
            }
        }
    }
    
exit_cleanup:
    potLinsys->backUpLin = 0;
    return retcode;
}

extern pot_int potLinsysSolve( pot_linsys *potLinsys, double *rhsVec, double *solVec ) {
    
    pot_int retcode = RETCODE_OK;
    if ( solVec ) {
        POTLP_MEMCPY(solVec, rhsVec, double, potLinsys->nCol);
        retcode = potLinsys->LSolve(potLinsys->solver, solVec);
    } else {
        retcode = potLinsys->LSolve(potLinsys->solver, rhsVec);
    }
    
    return retcode;
}

extern void potLinsysSwitchToBackup( pot_linsys *potLinsys ) {
    
    potLinsys->backUpLin = 1;
}

extern void potLinsysClear( pot_linsys *potLinsys ) {
    
    if ( !potLinsys ) {
        return;
    }
    
    potLinsys->LDestroy(&potLinsys->solver);
    POTLP_ZERO(potLinsys, pot_linsys, 1);
    return;
}

extern void potLinsysDestroy( pot_linsys **ppotLinsys ) {
    
    if ( !ppotLinsys ) {
        return;
    }
    
    potLinsysClear(*ppotLinsys);
    POTLP_FREE(*ppotLinsys);
    
    return;
}
