#ifdef HEADERPATH
#include "interface/hdsdp_utils.h"
#include "linalg/def_hdsdp_linsolver.h"
#include "linalg/hdsdp_linsolver.h"
#include "linalg/dense_opts.h"
#include "linalg/vec_opts.h"
#include "external/hdsdp_cs.h"
#else
#include "hdsdp_utils.h"
#include "def_hdsdp_linsolver.h"
#include "hdsdp_linsolver.h"
#include "dense_opts.h"
#include "vec_opts.h"
#include "hdsdp_cs.h"
#endif

#include <math.h>

#ifdef HDSDP_LINSYS_PROFILE
static int iPrintPardiso = 1;
static double dPardisoNumeric = 0.0;
static double dPardisoPsdCheck = 0.0;
static double dPardisoForward = 0.0;
static double dPardisoBackward = 0.0;
static double dPardisoSolve = 0.0;
static double dPardisoInvert = 0.0;
static int nPardisoNumeric = 0;
static int nPardisoPsdCheck = 0;
static int nPardisoForward = 0;
static int nPardisoBackward = 0;
static int nPardisoSolve = 0;
static int nPardisoInvert = 0;

static int iPrintQdldl = 1;
static double dQdldlNumeric = 0.0;
static double dQdldlPsdCheck = 0.0;
static double dQdldlForward = 0.0;
static double dQdldlBackward = 0.0;
static double dQdldlSolve = 0.0;
static double dQdldlInvert = 0.0;
static int nQdldlNumeric = 0;
static int nQdldlPsdCheck = 0;
static int nQdldlForward = 0;
static int nQdldlBackward = 0;
static int nQdldlSolve = 0;
static int nQdldlInvert = 0;

static int iPrintLapack = 1;
static double dLapackNumeric = 0.0;
static double dLapackPsdCheck = 0.0;
static double dLapackForward = 0.0;
static double dLapackBackward = 0.0;
static double dLapackSolve = 0.0;
static double dLapackInvert = 0.0;
static int nLapackNumeric = 0;
static int nLapackPsdCheck = 0;
static int nLapackForward = 0;
static int nLapackBackward = 0;
static int nLapackSolve = 0;
static int nLapackInvert = 0;
#endif

/* Sparse direct solver interface */
#ifdef LINSYS_PARDISO
static hdsdp_retcode pardisoLinSolverCreate( void **pchol, int nCol ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    HDSDP_NULLCHECK(pchol);
    pardiso_linsys *pds = NULL;
    HDSDP_INIT(pds, pardiso_linsys, 1);
    HDSDP_MEMCHECK(pds);
    
    pds->nCol = nCol;
    HDSDP_ZERO(pds->pt, void *, 64);
    HDSDP_ZERO(pds->iparm, int, 64);
    
    int mType = PARDISO_SYM_POSDEFINITE;
    pds->mType = mType;
    pardisoinit(pds->pt, &mType, pds->iparm);
    
    /* Non-default setting */
    set_pardiso_param(pds->iparm, PARDISO_PARAM_NONDEFAULT, 1);
    /* Symbolic ordering */
    set_pardiso_param(pds->iparm, PARDISO_PARAM_SYMBOLIC, PARDISO_PARAM_SYMBOLIC_MMD);
    /* Pivoting perturbation */
    set_pardiso_param(pds->iparm, PARDISO_PARAM_PERTURBATION, 13);
    /* Solve in-place*/
    set_pardiso_param(pds->iparm, PARDISO_PARAM_INPLACE, 1);
    /* Refinement */
    set_pardiso_param(pds->iparm, PARDISO_PARAM_REFINEMENT, 0);
    /* Use 0-based index*/
    set_pardiso_param(pds->iparm, PARDISO_PARAM_INDEX, PARDISO_PARAM_INDEX_C);
    /* Get diagonal elements */
    set_pardiso_param(pds->iparm, PARDISO_PARAM_DIAGONAL, PARDISO_PARAM_DIAGONAL_ON);
    
    *pchol = (void *) pds;
    
exit_cleanup:
    return retcode;
}

static hdsdp_retcode pardisoLinSolverCreateIndefinite( void **pchol, int nCol ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    HDSDP_NULLCHECK(pchol);
    pardiso_linsys *pds = NULL;
    HDSDP_INIT(pds, pardiso_linsys, 1);
    HDSDP_MEMCHECK(pds);
    
    pds->nCol = nCol;
    HDSDP_ZERO(pds->pt, void *, 64);
    HDSDP_ZERO(pds->iparm, int, 64);
    
    int mType = PARDISO_SYM_INDEFINITE;
    pds->mType = mType;
    pardisoinit(pds->pt, &mType, pds->iparm);
    
    set_pardiso_param(pds->iparm, PARDISO_PARAM_NONDEFAULT, 1);
    set_pardiso_param(pds->iparm, PARDISO_PARAM_SYMBOLIC, PARDISO_PARAM_SYMBOLIC_MMD);
    
    /* More perturbation for indefinite matrices */
    set_pardiso_param(pds->iparm, PARDISO_PARAM_PERTURBATION, 13);
    set_pardiso_param(pds->iparm, PARDISO_PARAM_INPLACE, 1);
    set_pardiso_param(pds->iparm, PARDISO_PARAM_REFINEMENT, 0);
    set_pardiso_param(pds->iparm, PARDISO_PARAM_INDEX, PARDISO_PARAM_INDEX_C);
    set_pardiso_param(pds->iparm, PARDISO_PARAM_DIAGONAL, PARDISO_PARAM_DIAGONAL_ON);
    
    *pchol = pds;
    
exit_cleanup:
    return retcode;
}

static void pardisoLinSolverSetThreads( void *chol, void *pThreads ) {
    
    pardiso_linsys *pds = (pardiso_linsys *) chol;
    
    int nThreads = *((int *) pThreads);
    set_pardiso_param(pds->iparm, PARDISO_PARAM_THREADS, nThreads);
    
    return;
}

#define SHOW_ORDERING(ord) printf("Using ordering %s. \n", ord)
static hdsdp_retcode pardisoLinSolverSymbolic( void *chol, int *colMatBeg, int *colMatIdx ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    pardiso_linsys *pds = (pardiso_linsys *) chol;
    
    pds->colMatBeg = colMatBeg;
    pds->colMatIdx = colMatIdx;
    
    HDSDP_INIT(pds->dWork, double, pds->nCol);
    HDSDP_MEMCHECK(pds->dWork);
    
    int maxfct = 1, mnum = 1, mtype = pds->mType, phase = PARDISO_PHASE_SYM;
    int idummy = 0, msg = 0, pdsret = PARDISO_RET_OK;
    
    /* Use amd first */
    int amdFactorNnz = 0;
    set_pardiso_param(pds->iparm, PARDISO_PARAM_SYMBOLIC, PARDISO_PARAM_SYMBOLIC_MMD);
    set_pardiso_param(pds->iparm, PARDISO_PARAM_FACNNZ, -1);
    
    pardiso(pds->pt, &maxfct, &mnum, &mtype, &phase,
            &pds->nCol, NULL, colMatBeg, colMatIdx, &idummy, &idummy,
            pds->iparm, &msg, NULL, NULL, &pdsret);
    
    if ( pdsret != PARDISO_RET_OK ) {
        retcode = HDSDP_RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    amdFactorNnz = get_pardiso_output(pds->iparm, PARDISO_PARAM_FACNNZ);
    
    /* Try nested dissection */
    int ndFactorNnz = 0;
    set_pardiso_param(pds->iparm, PARDISO_PARAM_SYMBOLIC, PARDISO_PARAM_SYMBOLIC_ND);
    set_pardiso_param(pds->iparm, PARDISO_PARAM_FACNNZ, -1);
    
    pardiso(pds->pt, &maxfct, &mnum, &mtype, &phase,
            &pds->nCol, NULL, colMatBeg, colMatIdx, &idummy, &idummy,
            pds->iparm, &msg, NULL, NULL, &pdsret);
    
    if ( pdsret != PARDISO_RET_OK ) {
        retcode = HDSDP_RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    ndFactorNnz = get_pardiso_output(pds->iparm, PARDISO_PARAM_FACNNZ);
    
    if ( ndFactorNnz > amdFactorNnz ) {
        /* Finally decide to use amd */
        set_pardiso_param(pds->iparm, PARDISO_PARAM_SYMBOLIC, PARDISO_PARAM_SYMBOLIC_MMD);
        set_pardiso_param(pds->iparm, PARDISO_PARAM_FACNNZ, -1);
        pardiso(pds->pt, &maxfct, &mnum, &mtype, &phase,
                &pds->nCol, NULL, colMatBeg, colMatIdx, &idummy, &idummy,
                pds->iparm, &msg, NULL, NULL, &pdsret);
#ifdef HDSDP_LINSYS_DEBUG
        printf("Using AMD ordering. nFactorNnz = %d \n", amdFactorNnz);
#endif
    } else {
#ifdef HDSDP_LINSYS_DEBUG
        printf("Using ND ordering. nFactorNnz = %d \n", ndFactorNnz);
#endif
    }
    
exit_cleanup:
    
    return retcode;
}

/* Stable Cholesky solve */
static hdsdp_retcode pardisoLinSolverStableNumeric( void *chol, int *colMatBeg, int *colMatIdx, double *colMatElem ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    pardiso_linsys *pds = (pardiso_linsys *) chol;
    
    pds->colMatBeg = colMatBeg;
    pds->colMatIdx = colMatIdx;
    pds->colMatElem = colMatElem;
    
    int maxfct = 1, mnum = 1, mtype = pds->mType, phase = PARDISO_PHASE_FAC;
    int idummy = 0, msg = 0, pdsret = PARDISO_RET_OK;
    
    set_pardiso_param(pds->iparm, PARDISO_PARAM_SCALING, 1);
    set_pardiso_param(pds->iparm, PARDISO_PARAM_MATCHING, 1);
    
    pardiso(pds->pt, &maxfct, &mnum, &mtype, &phase,
            &pds->nCol, colMatElem, colMatBeg, colMatIdx, &idummy, &idummy,
            pds->iparm, &msg, NULL, NULL, &pdsret);
    
    if ( pdsret != PARDISO_RET_OK ) {
        retcode = HDSDP_RETCODE_FAILED;
        goto exit_cleanup;
    }

exit_cleanup:
    return retcode;
}

static hdsdp_retcode pardisoLinSolverNumeric( void *chol, int *colMatBeg, int *colMatIdx, double *colMatElem ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    pardiso_linsys *pds = (pardiso_linsys *) chol;
    
    pds->colMatBeg = colMatBeg;
    pds->colMatIdx = colMatIdx;
    pds->colMatElem = colMatElem;
    
    int maxfct = 1, mnum = 1, mtype = pds->mType, phase = PARDISO_PHASE_FAC;
    int idummy = 0, msg = 0, pdsret = PARDISO_RET_OK;
    
#ifdef HDSDP_LINSYS_PROFILE
    double dTBeg = HUtilGetTimeStamp();
#endif
    
    pardiso(pds->pt, &maxfct, &mnum, &mtype, &phase,
            &pds->nCol, colMatElem, colMatBeg, colMatIdx, &idummy, &idummy,
            pds->iparm, &msg, NULL, NULL, &pdsret);
    
#ifdef HDSDP_LINSYS_PROFILE
    dPardisoNumeric += HUtilGetTimeStamp() - dTBeg;
    nPardisoNumeric += 1;
#endif
    
    if ( pdsret != PARDISO_RET_OK ) {
        retcode = HDSDP_RETCODE_FAILED;
        goto exit_cleanup;
    }

exit_cleanup:
    return retcode;
}

static hdsdp_retcode pardisoLinSolverPsdCheck( void *chol, int *colMatBeg, int *colMatIdx, double *colMatElem, int *isPsd ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    pardiso_linsys *pds = (pardiso_linsys *) chol;
    
    pds->colMatBeg = colMatBeg;
    pds->colMatIdx = colMatIdx;
    pds->colMatElem = colMatElem;
    
    /* Indefinite system does not have PD check */
    if ( pds->mType != PARDISO_SYM_POSDEFINITE ) {
        return HDSDP_RETCODE_FAILED;
    }
    
    int maxfct = 1, mnum = 1, mtype = pds->mType, phase = PARDISO_PHASE_FAC;
    int idummy = 0, msg = 0, pdsret = PARDISO_RET_OK;
    
#ifdef HDSDP_LINSYS_PROFILE
    double dTBeg = HUtilGetTimeStamp();
#endif
    
    pardiso(pds->pt, &maxfct, &mnum, &mtype, &phase,
            &pds->nCol, colMatElem, colMatBeg, colMatIdx, &idummy, &idummy,
            pds->iparm, &msg, NULL, NULL, &pdsret);
    
#ifdef HDSDP_LINSYS_PROFILE
    dPardisoPsdCheck += HUtilGetTimeStamp() - dTBeg;
    nPardisoPsdCheck += 1;
#endif
    
    if ( pdsret == PARDISO_RET_OK ) {
        *isPsd = 1;
    } else if ( pdsret == PARDISO_RET_INDEFINITE ) {
        *isPsd = 0;
    } else {
        retcode = HDSDP_RETCODE_FAILED;
        goto exit_cleanup;
    }

exit_cleanup:
    return retcode;
}

/* If nRhs > 1, then the solution must be overwritten */
static void pardisoLinSolverForwardN( void *chol, int nRhs, double *rhsVec, double *solVec ) {
    
    pardiso_linsys *pds = (pardiso_linsys *) chol;
    
    int maxfct = 1, mnum = 1, mtype = pds->mType, phase = PARDISO_PHASE_FORWARD;
    int idummy = 0, msg = 0, pdsret = PARDISO_RET_OK;
    
#ifdef HDSDP_LINSYS_PROFILE
    double dTBeg = HUtilGetTimeStamp();
#endif
    
    /* Do not overwrite */
    if ( solVec ) {
        set_pardiso_param(pds->iparm, PARDISO_PARAM_INPLACE, 0);
        pardiso(pds->pt, &maxfct, &mnum, &mtype, &phase,
                &pds->nCol, pds->colMatElem, pds->colMatBeg, pds->colMatIdx, &idummy, &nRhs,
                pds->iparm, &msg, rhsVec, solVec, &pdsret);
    } else {
        
        assert( nRhs == 1 );
        set_pardiso_param(pds->iparm, PARDISO_PARAM_INPLACE, 1);
        pardiso(pds->pt, &maxfct, &mnum, &mtype, &phase,
                &pds->nCol, pds->colMatElem, pds->colMatBeg, pds->colMatIdx, &idummy, &nRhs,
                pds->iparm, &msg, rhsVec, pds->dWork, &pdsret);
    }
    
#ifdef HDSDP_LINSYS_PROFILE
    dPardisoForward += HUtilGetTimeStamp() - dTBeg;
    nPardisoForward += 1;
#endif
    
    return;
}

static void pardisoLinSolverBackwardN( void *chol, int nRhs, double *rhsVec, double *solVec ) {
    
    pardiso_linsys *pds = (pardiso_linsys *) chol;
    
    int maxfct = 1, mnum = 1, mtype = pds->mType, phase = PARDISO_PHASE_BACKWARD;
    int idummy = 0, msg = 0, pdsret = PARDISO_RET_OK;
    
#ifdef HDSDP_LINSYS_PROFILE
    double dTBeg = HUtilGetTimeStamp();
#endif
    
    /* Do not overwrite */
    if ( solVec ) {
        set_pardiso_param(pds->iparm, PARDISO_PARAM_INPLACE, 0);
        pardiso(pds->pt, &maxfct, &mnum, &mtype, &phase,
                &pds->nCol, pds->colMatElem, pds->colMatBeg, pds->colMatIdx, &idummy, &nRhs,
                pds->iparm, &msg, rhsVec, solVec, &pdsret);
    } else {
        
        assert( nRhs == 1 );
        set_pardiso_param(pds->iparm, PARDISO_PARAM_INPLACE, 1);
        pardiso(pds->pt, &maxfct, &mnum, &mtype, &phase,
                &pds->nCol, pds->colMatElem, pds->colMatBeg, pds->colMatIdx, &idummy, &nRhs,
                pds->iparm, &msg, rhsVec, pds->dWork, &pdsret);
    }
    
#ifdef HDSDP_LINSYS_PROFILE
    dPardisoBackward += HUtilGetTimeStamp() - dTBeg;
    nPardisoBackward += 1;
#endif
    
    return;
}

static hdsdp_retcode pardisoLinSolverSolveN( void *chol, int nRhs, double *rhsVec, double *solVec ) {
    
    pardiso_linsys *pds = (pardiso_linsys *) chol;
    
    int maxfct = 1, mnum = 1, mtype = pds->mType, phase = PARDISO_PHASE_SOLVE;
    int idummy = 0, msg = 0, pdsret = PARDISO_RET_OK;
    
#ifdef HDSDP_LINSYS_PROFILE
    double dTBeg = HUtilGetTimeStamp();
#endif
    
    /* Do not overwrite */
    if ( solVec ) {
        set_pardiso_param(pds->iparm, PARDISO_PARAM_INPLACE, 0);
        pardiso(pds->pt, &maxfct, &mnum, &mtype, &phase,
                &pds->nCol, pds->colMatElem, pds->colMatBeg, pds->colMatIdx, &idummy, &nRhs,
                pds->iparm, &msg, rhsVec, solVec, &pdsret);
    } else {
        
        assert( nRhs == 1 );
        set_pardiso_param(pds->iparm, PARDISO_PARAM_INPLACE, 1);
        pardiso(pds->pt, &maxfct, &mnum, &mtype, &phase,
                &pds->nCol, pds->colMatElem, pds->colMatBeg, pds->colMatIdx, &idummy, &nRhs,
                pds->iparm, &msg, rhsVec, pds->dWork, &pdsret);
    }
    
#ifdef HDSDP_LINSYS_PROFILE
    dPardisoSolve += HUtilGetTimeStamp() - dTBeg;
    nPardisoSolve += 1;
#endif
    
    return HDSDP_RETCODE_OK;
}

static hdsdp_retcode pardisoLinSolverGetDiag( void *chol, double *diagVec ) {
    
    /* Get diagonal of Cholesky factor */
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    pardiso_linsys *pds = (pardiso_linsys *) chol;
    
    int one = 1, pdsret = PARDISO_RET_OK;
    pardiso_getdiag((const void **) pds->pt, (void *) diagVec, (void *) pds->dWork, &one, &pdsret);
    
    if ( pdsret != PARDISO_RET_OK ) {
        retcode = HDSDP_RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    /* For pardiso we need to take squareroot */
    for ( int i = 0; i < pds->nCol; ++i ) {
        diagVec[i] = sqrt(diagVec[i]);
    }
    
exit_cleanup:
    return retcode;
}

static void pardisoLinSolverInvert( void *chol, double *dFullMatrix, double *dAuxiMatrix ) {
    
    pardiso_linsys *pds = (pardiso_linsys *) chol;
    HDSDP_ZERO(dAuxiMatrix, double, pds->nCol * pds->nCol);
    
    double *pElem = dAuxiMatrix;
    for ( int iRow = 0; iRow < pds->nCol; ++iRow ) {
        *pElem = 1.0;
        pElem += pds->nCol + 1;
    }
    
#ifdef HDSDP_LINSYS_PROFILE
    double dTBeg = HUtilGetTimeStamp();
#endif
    
    pardisoLinSolverSolveN(chol, pds->nCol, dAuxiMatrix, dFullMatrix);
    
#ifdef HDSDP_LINSYS_PROFILE
    dPardisoInvert += HUtilGetTimeStamp() - dTBeg;
    nPardisoInvert += 1;
#endif
    
    return;
}

static void pardisoLinSolverClear( void *chol ) {
    
    if ( !chol ) {
        return;
    }
    
    pardiso_linsys *pds = (pardiso_linsys *) chol;
    
    int maxfct = 1, mnum = 1, mtype = pds->mType, phase = PARDISO_PHASE_FREE;
    int idummy = 0, msg = 0, pdsret = PARDISO_RET_OK;
    
    pardiso(pds->pt, &maxfct, &mnum, &mtype, &phase,
            &pds->nCol, pds->colMatElem, pds->colMatBeg, pds->colMatIdx, &idummy, NULL,
            pds->iparm, &msg, NULL, pds->dWork, &pdsret);
    
    HDSDP_FREE(pds->dWork);
    HDSDP_ZERO(pds, pardiso_linsys, 1);
    
    return;
}

static void pardisoLinSolverDestroy( void **pchol ) {
    
    if ( !pchol ) {
        return;
    }
    
    pardisoLinSolverClear(*pchol);
    HDSDP_FREE(*pchol);
    
    return;
}
#else

/* QDLDL solver for Quasi-definite systems */
static hdsdp_retcode qdldlLinSolverCreate( void **pchol, int nCol ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    HDSDP_NULLCHECK(pchol);
    qdldl_linsys *qdldl = NULL;
    HDSDP_INIT(qdldl, qdldl_linsys, 1);
    HDSDP_MEMCHECK(qdldl);
    
    qdldl->nCol = nCol;
    *pchol = qdldl;
    
exit_cleanup:
    return retcode;
}

static void qdldlLinSolverSetThreads( void *chol, void *dummy ) {
    
    (void) dummy;
    
    return;
}

static hdsdp_retcode qdldlLinSolverSymbolic( void *chol, int *colMatBeg, int *colMatIdx ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    qdldl_linsys *qdldl = (qdldl_linsys *) chol;
    /* Allocate memory for factorization */
    HDSDP_INIT(qdldl->P, int, qdldl->nCol);
    HDSDP_INIT(qdldl->iWork, int, qdldl->nCol * 4);
    HDSDP_INIT(qdldl->dWork, double, qdldl->nCol * 4);
    HDSDP_INIT(qdldl->Lp, int, qdldl->nCol + 1);
    HDSDP_INIT(qdldl->Lnz, int, qdldl->nCol);
    
    HDSDP_MEMCHECK(qdldl->P);
    HDSDP_MEMCHECK(qdldl->iWork);
    HDSDP_MEMCHECK(qdldl->dWork);
    HDSDP_MEMCHECK(qdldl->Lp);
    HDSDP_MEMCHECK(qdldl->Lnz);
        
    /* Reconstruct the upper-triangular part of the matrix */
    int nNz = colMatBeg[qdldl->nCol];
    HDSDP_INIT(qdldl->colMatUpperBeg, int, qdldl->nCol + 1);
    HDSDP_INIT(qdldl->colMatUpperIdx, int, nNz);
    HDSDP_INIT(qdldl->colMatUpperElem, double, nNz);
    HDSDP_INIT(qdldl->iTransMap, int, nNz);
    
    HDSDP_MEMCHECK(qdldl->colMatUpperBeg);
    HDSDP_MEMCHECK(qdldl->colMatUpperIdx);
    HDSDP_MEMCHECK(qdldl->colMatUpperElem);
    HDSDP_MEMCHECK(qdldl->iTransMap);
    
    int *iTransBuffer = qdldl->iWork;
    /* Transpose the matrix to obtain the upper-triangular part */
    for ( int iElem = 0; iElem < colMatBeg[qdldl->nCol]; ++iElem ) {
        iTransBuffer[colMatIdx[iElem]] += 1;
    }
    
    dcs_cumsum(qdldl->colMatUpperBeg, iTransBuffer, qdldl->nCol);
    
    int iPos = 0;
    for ( int iCol = 0; iCol < qdldl->nCol; ++iCol ) {
        for ( int iElem = colMatBeg[iCol]; iElem < colMatBeg[iCol + 1]; ++iElem ) {
            qdldl->colMatUpperIdx[iPos = iTransBuffer[colMatIdx[iElem]]++] = iCol;
            qdldl->iTransMap[iElem] = iPos;
        }
    }
    
    /* Symbolic factorization */
    int ldlret = QDLDL_etree(qdldl->nCol, qdldl->colMatUpperBeg, qdldl->colMatUpperIdx,
                             qdldl->iWork, qdldl->Lnz, qdldl->P);
    
    if ( ldlret < 0 ) {
        retcode = HDSDP_RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    /* Allocate memory for numerical factorization */
    HDSDP_INIT(qdldl->Li, int, ldlret);
    HDSDP_INIT(qdldl->Lx, double, ldlret);
    HDSDP_MEMCHECK(qdldl->Li);
    HDSDP_MEMCHECK(qdldl->Lx);
    
exit_cleanup:
    return retcode;
}

static hdsdp_retcode qdldlLinSolverNumeric( void *chol, int *colMatBeg, int *colMatIdx, double *colMatElem ) {
    /* Given that the internal mapping from lower to upper triangular is fixed after the symbolic factorization
       phase. colMatBeg and colMatIdx are assumed to be unchanged */
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    qdldl_linsys *qdldl = (qdldl_linsys *) chol;
    
    int nCol = qdldl->nCol;
    double *D = qdldl->dWork;
    double *Dinv = D + nCol;
    double *fWork = Dinv + nCol;
    
    int *bWork = qdldl->iWork;
    int *iWork = bWork + nCol;
    
    /* Map upper to lower-triangular */
    for ( int iElem = 0; iElem < qdldl->colMatUpperBeg[nCol]; ++iElem ) {
        qdldl->colMatUpperElem[qdldl->iTransMap[iElem]] = colMatElem[iElem];
    }
    
    /* Numerical factorization */
    int nPos = QDLDL_factor(qdldl->nCol, qdldl->colMatUpperBeg, qdldl->colMatUpperIdx,
                            qdldl->colMatUpperElem, qdldl->Lp, qdldl->Li, qdldl->Lx,
                            D, Dinv, qdldl->Lnz, qdldl->P, bWork, iWork, fWork);
    
    if ( nPos == -1 ) {
        retcode = HDSDP_RETCODE_FAILED;
        goto exit_cleanup;
    }
    
exit_cleanup:
    return retcode;
}

static hdsdp_retcode qdldlLinSolverPsdCheck( void *chol, int *colMatBeg, int *colMatIdx, double *colMatElem, int *isPsd ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    qdldl_linsys *qdldl = (qdldl_linsys *) chol;
    
    int nCol = qdldl->nCol;
    double *D = qdldl->dWork;
    double *Dinv = D + nCol;
    double *Dinvsqrt = Dinv + nCol;
    double *fWork = Dinvsqrt + nCol;
    
    int *bWork = qdldl->iWork;
    int *iWork = bWork + nCol;
    
    /* Map upper to lower-triangular */
    for ( int iElem = 0; iElem < qdldl->colMatUpperBeg[nCol]; ++iElem ) {
        qdldl->colMatUpperElem[qdldl->iTransMap[iElem]] = colMatElem[iElem];
    }
    
    /* Numerical factorization */
    int nPos = QDLDL_factor(qdldl->nCol, qdldl->colMatUpperBeg, qdldl->colMatUpperIdx,
                            qdldl->colMatUpperElem, qdldl->Lp, qdldl->Li, qdldl->Lx,
                            D, Dinv, qdldl->Lnz, qdldl->P, bWork, iWork, fWork);
    
    if ( nPos < nCol ) {
        *isPsd = 0;
    } else {
        *isPsd = 1;
    }
    
    /* Create Dinvsqrt */
    for ( int iCol = 0; iCol < nCol; ++iCol ) {
        Dinvsqrt[iCol] = sqrt(Dinv[iCol]);
    }
    
exit_cleanup:
    return retcode;
}

static void qdldlLinSolverForwardN( void *chol, int nRhs, double *rhsVec, double *solVec ) {

    qdldl_linsys *qdldl = (qdldl_linsys *) chol;
    
    double *dSol = NULL;
    double *Dinvsqrt = qdldl->dWork + 2 * qdldl->nCol;
    
    if ( solVec ) {
        dSol = solVec;
        HDSDP_MEMCPY(dSol, rhsVec, double, qdldl->nCol * nRhs);
    } else {
        assert( nRhs == 1 );
        dSol = rhsVec;
    }
    
    /* Solve and scale by 1./sqrt(D) */
    for ( int iRhs = 0; iRhs < nRhs; ++iRhs ) {
        QDLDL_Lsolve(qdldl->nCol, qdldl->Lp, qdldl->Li, qdldl->Lx, dSol);
        for ( int iCol = 0; iCol < qdldl->nCol; ++iCol ) {
            dSol[iCol] = dSol[iCol] * Dinvsqrt[iCol];
        }
        dSol += qdldl->nCol;
    }
    
    return;
}

static void qdldlLinSolverBackwardN( void *chol, int nRhs, double *rhsVec, double *solVec ) {

    qdldl_linsys *qdldl = (qdldl_linsys *) chol;
    
    double *dSol = NULL;
    double *Dinvsqrt = qdldl->dWork + 2 * qdldl->nCol;
    
    if ( solVec ) {
        dSol = solVec;
        HDSDP_MEMCPY(dSol, rhsVec, double, qdldl->nCol * nRhs);
    } else {
        assert( nRhs == 1 );
        dSol = rhsVec;
    }
    
    /* Scale by 1./sqrt(D) and solve */
    for ( int iRhs = 0; iRhs < nRhs; ++iRhs ) {
        for ( int iCol = 0; iCol < qdldl->nCol; ++iCol ) {
            dSol[iCol] = dSol[iCol] * Dinvsqrt[iCol];
        }
        QDLDL_Ltsolve(qdldl->nCol, qdldl->Lp, qdldl->Li, qdldl->Lx, dSol);
        dSol += qdldl->nCol;
    }
    
    return;
}

static hdsdp_retcode qdldlLinSolverSolveN( void *chol, int nRhs, double *rhsVec, double *solVec ) {

    qdldl_linsys *qdldl = (qdldl_linsys *) chol;
    double *Dinv = qdldl->dWork + qdldl->nCol;
    
    double *dSol = NULL;
    if ( solVec ) {
        dSol = solVec;
        HDSDP_MEMCPY(dSol, rhsVec, double, qdldl->nCol * nRhs);
    } else {
        assert( nRhs == 1 );
        dSol = rhsVec;
    }
    
    /* Solve */
    for ( int iRhs = 0; iRhs < nRhs; ++iRhs ) {
        QDLDL_solve(qdldl->nCol, qdldl->Lp, qdldl->Li, qdldl->Lx, Dinv, dSol);
        dSol += qdldl->nCol;
    }
    
    return HDSDP_RETCODE_OK;
}

static hdsdp_retcode qdldlLinSolverGetDiag( void *chol, double *diagVec ) {

    qdldl_linsys *qdldl = (qdldl_linsys *) chol;
    double *D = qdldl->dWork;
    
    for ( int iCol = 0; iCol < qdldl->nCol; ++iCol ) {
        diagVec[iCol] = sqrt(D[iCol]);
    }
    
    return HDSDP_RETCODE_OK;
}

static void qdldlLinSolverInvert( void *chol, double *dFullMatrix, double *dAuxiMatrix ) {
    
    qdldl_linsys *qdldl = (qdldl_linsys *) chol;
    HDSDP_ZERO(dAuxiMatrix, double, qdldl->nCol * qdldl->nCol);
    
    double *pElem = dAuxiMatrix;
    for ( int iRow = 0; iRow < qdldl->nCol; ++iRow ) {
        *pElem = 1.0;
        pElem += qdldl->nCol + 1;
    }
    
    qdldlLinSolverSolveN(chol, qdldl->nCol, dAuxiMatrix, dFullMatrix);
    
    return;
}

static void qdldlLinSolverClear( void *chol ) {
    
    if ( !chol ) {
        return;
    }
    
    qdldl_linsys *qdldl = (qdldl_linsys *) chol;
    
    HDSDP_FREE(qdldl->colMatUpperBeg);
    HDSDP_FREE(qdldl->colMatUpperIdx);
    HDSDP_FREE(qdldl->colMatUpperElem);
    HDSDP_FREE(qdldl->iTransMap);
    HDSDP_FREE(qdldl->iWork);
    HDSDP_FREE(qdldl->dWork);
    HDSDP_FREE(qdldl->P);
    HDSDP_FREE(qdldl->Lnz);
    HDSDP_FREE(qdldl->Lp);
    HDSDP_FREE(qdldl->Li);
    HDSDP_FREE(qdldl->Lx);
    
    HDSDP_ZERO(qdldl, qdldl_linsys, 1);
    
    return;
}

static void qdldlLinSolverDestroy( void **pchol ) {
    
    if ( !pchol ) {
        return;
    }
    
    qdldlLinSolverClear(*pchol);
    HDSDP_FREE(*pchol);
    
    return;
}

/* LDL for indefinite systems */
static hdsdp_retcode ldlLinSolverCreate( void **pchol, int nCol ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    HDSDP_NULLCHECK(pchol);
    ldl_linsys *ldl = NULL;
    HDSDP_INIT(ldl, ldl_linsys, 1);
    HDSDP_MEMCHECK(ldl);
    
    ldl->nCol = nCol;
    *pchol = ldl;
    
exit_cleanup:
    return retcode;
}

static void ldlLinSolverSetThreads( void *chol, void *dummy ) {
    
    (void) dummy;
    
    return;
}

static hdsdp_retcode ldlLinSolverSymbolic( void *chol, int *colMatBeg, int *colMatIdx ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    ldl_linsys *ldl = (ldl_linsys *) chol;
    
    HDSDP_INIT(ldl->Parent, int, ldl->nCol);
    HDSDP_INIT(ldl->Flag, int, ldl->nCol);
    HDSDP_INIT(ldl->Lp, int, ldl->nCol + 1);
    HDSDP_INIT(ldl->D, double, ldl->nCol);
    HDSDP_INIT(ldl->Y, double, ldl->nCol);
    HDSDP_INIT(ldl->Pattern, int, ldl->nCol);
    HDSDP_INIT(ldl->Lnz, int, ldl->nCol);
    
    HDSDP_MEMCHECK(ldl->Parent);
    HDSDP_MEMCHECK(ldl->Flag);
    HDSDP_MEMCHECK(ldl->Y);
    HDSDP_MEMCHECK(ldl->Pattern);
    HDSDP_MEMCHECK(ldl->D);
    HDSDP_MEMCHECK(ldl->Lnz);
    HDSDP_MEMCHECK(ldl->Lp);
    
    int nNz = colMatBeg[ldl->nCol];
    
    HDSDP_INIT(ldl->colMatUpperBeg, int, ldl->nCol + 1);
    HDSDP_INIT(ldl->colMatUpperIdx, int, nNz);
    HDSDP_INIT(ldl->colMatUpperElem, double, nNz);
    HDSDP_INIT(ldl->iTransMap, int, nNz);
    
    HDSDP_MEMCHECK(ldl->colMatUpperBeg);
    HDSDP_MEMCHECK(ldl->colMatUpperIdx);
    HDSDP_MEMCHECK(ldl->colMatUpperElem);
    HDSDP_MEMCHECK(ldl->iTransMap);
    
    int *iTransBuffer = ldl->Parent;
        
    /* Transpose the matrix to obtain the upper-triangular part */
    for ( int iElem = 0; iElem < colMatBeg[ldl->nCol]; ++iElem ) {
        iTransBuffer[colMatIdx[iElem]] += 1;
    }
    
    dcs_cumsum(ldl->colMatUpperBeg, iTransBuffer, ldl->nCol);
    
    int iPos = 0;
    for ( int iCol = 0; iCol < ldl->nCol; ++iCol ) {
        for ( int iElem = colMatBeg[iCol]; iElem < colMatBeg[iCol + 1]; ++iElem ) {
            ldl->colMatUpperIdx[iPos = iTransBuffer[colMatIdx[iElem]]++] = iCol;
            ldl->iTransMap[iElem] = iPos;
        }
    }
    
    /* Symbolic factorization without reordering */
    ldl_symbolic(ldl->nCol, ldl->colMatUpperBeg, ldl->colMatUpperIdx,
                 ldl->Lp, ldl->Parent, ldl->Lnz, ldl->Flag, NULL, NULL);
    
    HDSDP_INIT(ldl->Li, int, ldl->Lp[ldl->nCol]);
    HDSDP_INIT(ldl->Lx, double, ldl->Lp[ldl->nCol]);
    HDSDP_MEMCHECK(ldl->Li);
    HDSDP_MEMCHECK(ldl->Lx);
    
exit_cleanup:
    return retcode;
}

static hdsdp_retcode ldlLinSolverNumeric( void *chol, int *colMatBeg, int *colMatIdx, double *colMatElem ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    ldl_linsys *ldl = (ldl_linsys *) chol;
    
    int nCol = ldl->nCol;
    
    for ( int iElem = 0; iElem < ldl->colMatUpperBeg[nCol]; ++iElem ) {
        ldl->colMatUpperElem[ldl->iTransMap[iElem]] = colMatElem[iElem];
    }
    
    int nPos = ldl_numeric(nCol, ldl->colMatUpperBeg, ldl->colMatUpperIdx,
                           ldl->colMatUpperElem, ldl->Lp, ldl->Parent, ldl->Lnz,
                           ldl->Li, ldl->Lx, ldl->D, ldl->Y, ldl->Pattern, ldl->Flag, NULL, NULL);
    
    if ( nPos != nCol ) {
        retcode = HDSDP_RETCODE_FAILED;
        goto exit_cleanup;
    }
    
exit_cleanup:
    return retcode;
}

static hdsdp_retcode ldlLinSolverPsdCheck( void *chol, int *dummy1, int *dummy2, double *dummy3, int *dummy4 ) {
    
    (void) chol;
    (void) dummy1;
    (void) dummy2;
    (void) dummy3;
    (void) dummy4;
    
    return HDSDP_RETCODE_FAILED;
}

static void ldlLinSolverForwardN( void *chol, int dummy1, double *dummy2, double *dummy3 ) {

    (void) chol;
    (void) dummy1;
    (void) dummy2;
    (void) dummy3;
    
    return;
}

static void ldlLinSolverBackwardN( void *chol, int dummy1, double *dummy2, double *dummy3 ) {

    (void) chol;
    (void) dummy1;
    (void) dummy2;
    (void) dummy3;
    
    return;
}

static hdsdp_retcode ldlLinSolverSolveN( void *chol, int nRhs, double *rhsVec, double *solVec ) {
    
    ldl_linsys *ldl = (ldl_linsys *) chol;
    
    double *dSol = NULL;
    if ( solVec ) {
        dSol = solVec;
        HDSDP_MEMCPY(dSol, rhsVec, double, ldl->nCol * nRhs);
    } else {
        assert( nRhs == 1 );
        dSol = rhsVec;
    }
    
    /* Solve */
    for ( int iRhs = 0; iRhs < nRhs; ++iRhs ) {
        ldl_lsolve(ldl->nCol, dSol, ldl->Lp, ldl->Li, ldl->Lx);
        ldl_dsolve(ldl->nCol, dSol, ldl->D);
        ldl_ltsolve(ldl->nCol, dSol, ldl->Lp, ldl->Li, ldl->Lx);
        dSol += ldl->nCol;
    }
    
    return HDSDP_RETCODE_OK;
}

static hdsdp_retcode ldlLinSolverGetDiag( void *chol, double *dummy ) {
    
    (void) chol;
    (void) dummy;
    
    return HDSDP_RETCODE_FAILED;
}

static void ldlLinSolverInvert( void *chol, double *dFullMatrix, double *dAuxiMatrix ) {
    
    ldl_linsys *ldl = (ldl_linsys *) chol;
    HDSDP_ZERO(dAuxiMatrix, double, ldl->nCol * ldl->nCol);
    
    double *pElem = dAuxiMatrix;
    for ( int iRow = 0; iRow < ldl->nCol; ++iRow ) {
        *pElem = 1.0;
        pElem += ldl->nCol + 1;
    }
    
    ldlLinSolverSolveN(chol, ldl->nCol, dAuxiMatrix, dFullMatrix);
    
    return;
}

static void ldlLinSolverClear( void *chol ) {
    
    if ( !chol ) {
        return;
    }
    
    ldl_linsys *ldl = (ldl_linsys *) chol;
    
    HDSDP_FREE(ldl->colMatUpperBeg);
    HDSDP_FREE(ldl->colMatUpperIdx);
    HDSDP_FREE(ldl->colMatUpperElem);
    
    HDSDP_FREE(ldl->iTransMap);
    HDSDP_FREE(ldl->Parent);
    HDSDP_FREE(ldl->Flag);
    HDSDP_FREE(ldl->D);
    HDSDP_FREE(ldl->Y);
    HDSDP_FREE(ldl->Pattern);
    
    HDSDP_FREE(ldl->Lnz);
    HDSDP_FREE(ldl->Lp);
    HDSDP_FREE(ldl->Li);
    HDSDP_FREE(ldl->Lx);
    
    HDSDP_ZERO(ldl, ldl_linsys, 1);
    
    return;
}

static void ldlLinSolverDestroy( void **pchol ) {
    
    if ( !pchol ) {
        return;
    }
    
    ldlLinSolverClear(*pchol);
    HDSDP_FREE(*pchol);
    
    return;
}
#endif

/* Dense direct solver interface */
static hdsdp_retcode lapackLinSolverCreate( void **pchol, int nCol ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    HDSDP_NULLCHECK(pchol);
    lapack_flinsys *lap = NULL;
    HDSDP_INIT(lap, lapack_flinsys, 1);
    HDSDP_MEMCHECK(lap);
    
    HDSDP_ZERO(lap, lapack_flinsys, 1);
    lap->nCol = nCol;
    
    HDSDP_INIT(lap->dFullMatElem, double, nCol * nCol);
    HDSDP_MEMCHECK(lap->dFullMatElem);
    
    *pchol = lap;
    
exit_cleanup:
    return retcode;
}

static void lapackLinSolverSetParamsDummy( void *chol, void *param ) {
    
    (void) chol;
    (void) param;
    
    return;
}

static hdsdp_retcode lapackLinSolverSymbolic( void *chol, int *dummy1, int *dummy2 ) {
    
    (void) chol;
    (void) dummy1;
    (void) dummy2;
    
    return HDSDP_RETCODE_OK;
}

static hdsdp_retcode lapackLinSolverNumeric( void *chol, int *dummy1, int *dummy2, double *dFullMatrix ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    lapack_flinsys *lap = (lapack_flinsys *) chol;
    
    HDSDP_MEMCPY(lap->dFullMatElem, dFullMatrix, double, lap->nCol * lap->nCol);
    
    int info = LAPACK_RET_OK;
    char uplolow = LAPACK_UPLOW_LOW;
    
#ifdef HDSDP_LINSYS_PROFILE
    double dTBeg = HUtilGetTimeStamp();
#endif
    
    dpotrf(&uplolow, &lap->nCol, lap->dFullMatElem, &lap->nCol, &info);
    
#ifdef HDSDP_LINSYS_PROFILE
    dLapackNumeric += HUtilGetTimeStamp() - dTBeg;
    nLapackNumeric += 1;
#endif
    
    if ( info != LAPACK_RET_OK ) {
        retcode = HDSDP_RETCODE_FAILED;
        goto exit_cleanup;
    }
    
exit_cleanup:
    return retcode;
}

static hdsdp_retcode lapackLinSolverPsdCheck( void *chol, int *dummy1, int *dummy2, double *dFullMatrix, int *isPsd ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    lapack_flinsys *lap = (lapack_flinsys *) chol;
    
    HDSDP_MEMCPY(lap->dFullMatElem, dFullMatrix, double, lap->nCol * lap->nCol);
    
    int info = LAPACK_RET_OK;
    char uplolow = LAPACK_UPLOW_LOW;
    
#ifdef HDSDP_LINSYS_PROFILE
    double dTBeg = HUtilGetTimeStamp();
#endif
    
    dpotrf(&uplolow, &lap->nCol, lap->dFullMatElem, &lap->nCol, &info);
    
#ifdef HDSDP_LINSYS_PROFILE
    dLapackPsdCheck += HUtilGetTimeStamp() - dTBeg;
    nLapackPsdCheck += 1;
#endif
    
    if ( info == LAPACK_RET_OK ) {
        *isPsd = 1;
    } else if ( info > 0 ) {
        *isPsd = 0;
    } else {
        retcode = HDSDP_RETCODE_FAILED;
        goto exit_cleanup;
    }
 
exit_cleanup:
    return retcode;
}

static void lapackLinSolverForwardN( void *chol, int nRhs, double *rhsVec, double *solVec ) {
    
    lapack_flinsys *lap = (lapack_flinsys *) chol;
    
    char uplolow = LAPACK_UPLOW_LOW, sideleft = LAPACK_SIDE_LEFT, notrans = LAPACK_NOTRANS, nonunit = LAPACK_DIAG_NONUNIT;
    double done = 1.0;
    
#ifdef HDSDP_LINSYS_PROFILE
    double dTBeg = HUtilGetTimeStamp();
#endif
    if ( solVec ) {
        HDSDP_MEMCPY(solVec, rhsVec, double, nRhs * lap->nCol);
        dtrsm(&sideleft, &uplolow, &notrans, &nonunit, &lap->nCol, &nRhs,
              &done, lap->dFullMatElem, &lap->nCol, solVec, &lap->nCol);
    } else {
        dtrsm(&sideleft, &uplolow, &notrans, &nonunit, &lap->nCol, &nRhs,
              &done, lap->dFullMatElem, &lap->nCol, rhsVec, &lap->nCol);
    }
#ifdef HDSDP_LINSYS_PROFILE
    dLapackForward += HUtilGetTimeStamp() - dTBeg;
    nLapackForward += 1;
#endif
 
    return;
}

static void lapackLinSolverBackwardN( void *chol, int nRhs, double *rhsVec, double *solVec ) {
    
    lapack_flinsys *lap = (lapack_flinsys *) chol;
    
    char uplolow = LAPACK_UPLOW_LOW, sideleft = LAPACK_SIDE_LEFT, trans = LAPACK_TRANS, nonunit = LAPACK_DIAG_NONUNIT;
    double done = 1.0;
    
#ifdef HDSDP_LINSYS_PROFILE
    double dTBeg = HUtilGetTimeStamp();
#endif
    if ( solVec ) {
        HDSDP_MEMCPY(solVec, rhsVec, double, nRhs * lap->nCol);
        dtrsm(&sideleft, &uplolow, &trans, &nonunit, &lap->nCol, &nRhs,
              &done, lap->dFullMatElem, &lap->nCol, solVec, &lap->nCol);
    } else {
        dtrsm(&sideleft, &uplolow, &trans, &nonunit, &lap->nCol, &nRhs,
              &done, lap->dFullMatElem, &lap->nCol, rhsVec, &lap->nCol);
    }
#ifdef HDSDP_LINSYS_PROFILE
    dLapackBackward += HUtilGetTimeStamp() - dTBeg;
    nLapackBackward += 1;
#endif
 
    return;
}

static hdsdp_retcode lapackLinSolverSolveN( void *chol, int nRhs, double *rhsVec, double *solVec ) {
    
    lapack_flinsys *lap = (lapack_flinsys *) chol;
    
    int info = LAPACK_RET_OK;
    char uplolow = LAPACK_UPLOW_LOW;
    
#ifdef HDSDP_LINSYS_PROFILE
    double dTBeg = HUtilGetTimeStamp();
#endif
    if ( solVec ) {
        HDSDP_MEMCPY(solVec, rhsVec, double, nRhs * lap->nCol);
        dpotrs(&uplolow, &lap->nCol, &nRhs, lap->dFullMatElem,
               &lap->nCol, solVec, &lap->nCol, &info);
    } else {
        dpotrs(&uplolow, &lap->nCol, &nRhs, lap->dFullMatElem,
               &lap->nCol, rhsVec, &lap->nCol, &info);
    }
    
#ifdef HDSDP_LINSYS_PROFILE
    dLapackSolve += HUtilGetTimeStamp() - dTBeg;
    nLapackSolve += 1;
#endif
    
    assert( info == LAPACK_RET_OK );
    
    return HDSDP_RETCODE_OK;
}

static hdsdp_retcode lapackLinSolverGetDiag( void *chol, double *diagVec ) {
    
    lapack_flinsys *lap = (lapack_flinsys *) chol;
    
    for ( int iRow = 0; iRow < lap->nCol; ++iRow ) {
        diagVec[iRow] = lap->dFullMatElem[iRow * lap->nCol + iRow];
    }
    
    return HDSDP_RETCODE_OK;
}

static void lapackLinSolverInvert( void *chol, double *dFullMatrix, double *dAuxiMatrix ) {
    
    (void) dAuxiMatrix;
    lapack_flinsys *lap = (lapack_flinsys *) chol;
    HDSDP_MEMCPY(dFullMatrix, lap->dFullMatElem, double, lap->nCol * lap->nCol);
    
    int info = LAPACK_RET_OK;
    char uplolow = LAPACK_UPLOW_LOW;
    
#ifdef HDSDP_LINSYS_PROFILE
    double dTBeg = HUtilGetTimeStamp();
#endif
    dpotri(&uplolow, &lap->nCol, dFullMatrix, &lap->nCol, &info);
    HUtilMatSymmetrize(lap->nCol, dFullMatrix);
#ifdef HDSDP_LINSYS_PROFILE
    dLapackInvert += HUtilGetTimeStamp() - dTBeg;
    nLapackInvert += 1;
#endif
    
    assert( info == LAPACK_RET_OK );
    
    return;
}

static void lapackLinSolverClear( void *chol ) {
    
    if ( !chol ) {
        return;
    }
    
    lapack_flinsys *lap = (lapack_flinsys *) chol;
    HDSDP_FREE(lap->dFullMatElem);
    
    HDSDP_ZERO(lap, lapack_flinsys, 1);
    
    return;
}

static void lapackLinSolverDestroy( void **pchol ) {
    
    if ( !pchol ) {
        return;
    }
    
    lapackLinSolverClear(*pchol);
    HDSDP_FREE(*pchol);
    
    return;
}

/* Dense iterative solver interface */
static hdsdp_retcode conjGradLinSolverCreate( void **pchol, int nCol ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    HDSDP_NULLCHECK(pchol);
    
    iterative_linsys *cg = NULL;
    HDSDP_INIT(cg, iterative_linsys, 1);
    HDSDP_MEMCHECK(cg);
    HDSDP_ZERO(cg, iterative_linsys, 1);
    
    cg->nCol = nCol;
    
    HDSDP_INIT(cg->iterResi, double, nCol);
    HDSDP_INIT(cg->iterResiNew, double, nCol);
    HDSDP_INIT(cg->iterDirection, double, nCol);
    HDSDP_INIT(cg->preInvResi, double, nCol);
    HDSDP_INIT(cg->MTimesDirection, double, nCol);
    HDSDP_INIT(cg->iterVec, double, nCol);
    HDSDP_INIT(cg->rhsBuffer, double, nCol);
    
    cg->useJacobi = 1;
    HDSDP_INIT(cg->JacobiPrecond, double, nCol);
    HDSDP_CALL(lapackLinSolverCreate((void **)&cg->lap, nCol));
    
    /* Set parameters */
    cg->params.absTol = 1e-06;
    cg->params.relTol = 1e-06;
    cg->params.useJacobi = 1;
    cg->params.maxIter = HDSDP_MAX(50, nCol / 20);
    cg->params.nRestartFreq = 20;
    
    *pchol = cg;
        
exit_cleanup:
    return retcode;
}

static void conjGradLinSolverSetParam( void *chol, void *iterParams ) {
    
    iterative_linsys *cg = (iterative_linsys *) chol;
    iterative_params *params = (iterative_params *) iterParams;
    
    if ( params->absTol > 0 ) {
        cg->params.absTol = params->absTol;
    }
    
    if ( params->relTol > 0 ) {
        cg->params.relTol = params->relTol;
    }
    
    if ( params->maxIter > 0 ) {
        cg->params.maxIter = params->maxIter;
    }
    
    if ( params->nRestartFreq > 0 )
    cg->params.nRestartFreq = params->nRestartFreq;
    
    return;
}

static hdsdp_retcode conjGradLinSolverSymbolic( void *chol, int *dummy1, int *dummy2 ) {
    
    (void) chol;
    (void) dummy1;
    (void) dummy2;
    
    return HDSDP_RETCODE_OK;
}

static hdsdp_retcode conjGradBuildPreconditioner( iterative_linsys *cg );
static hdsdp_retcode conjGradLinSolverNumeric( void *chol, int *dummy1, int *dummy2, double *dFullMatrix ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    (void) chol;
    (void) dummy1;
    (void) dummy2;
    
    iterative_linsys *cg = (iterative_linsys *) chol;
    cg->fullMatElem = dFullMatrix;
    HDSDP_CALL(conjGradBuildPreconditioner(cg));
    
exit_cleanup:
    return retcode;
}

static hdsdp_retcode conjGradLinSolverPsdCheck( void *chol, int *dummy1, int *dummy2, double *dummy3, int *dummy4 ) {
    
    (void) chol;
    (void) dummy1;
    (void) dummy2;
    (void) dummy3;
    (void) dummy4;
    
    return HDSDP_RETCODE_FAILED;
}

static void conjGradLinSolverForwardN( void *chol, int dummy1, double *dummy2, double *dummy3 ) {
    
    (void) chol;
    (void) dummy1;
    (void) dummy2;
    (void) dummy3;
    
    return;
}

static void conjGradLinSolverBackwardN( void *chol, int dummy1, double *dummy2, double *dummy3 ) {
    
    (void) chol;
    (void) dummy1;
    (void) dummy2;
    (void) dummy3;
    
    return;
}

static hdsdp_retcode conjGradBuildPreconditioner( iterative_linsys *cg ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    /*
     Jacobi preconditioner just uses the diagonal elements of a matrix
     https://netlib.org/linalg/html_templates/node55.html
     TODO: consider sparsity-approximate inverse method mentioned in
     https://math.stackexchange.com/questions/978051/find-diagonal-of-inverse-matrix
     */
    if ( cg->useJacobi ) {
        for ( int iCol = 0; iCol < cg->nCol; ++iCol ) {
            cg->JacobiPrecond[iCol] = cg->fullMatElem[iCol + iCol * cg->nCol];
        }
#ifdef HDSDP_CONJGRAD_DEBUG
        hdsdp_printf("Using Jacobi pre-conditioner. \n");
#endif
    } else {
        /* Use Cholesky as the preconditioner */
#ifdef HDSDP_CONJGRAD_DEBUG
        hdsdp_printf("Using Cholesky pre-conditioner. \n");
#endif
        retcode = lapackLinSolverNumeric((void *) cg->lap, NULL, NULL, cg->fullMatElem);
    }
    
    return retcode;
}

static void conjGradApplyPreconditioner( iterative_linsys *cg, double *rhsVec ) {
    
    /* Apply preconditioner rhsVec <- P \ rhsVec */
    
    if ( cg->useJacobi ) {
        vvrscl(&cg->nCol, cg->JacobiPrecond, rhsVec);
    } else {
        lapackLinSolverSolveN((void *) cg->lap, 1, rhsVec, NULL);
    }
    
    return;
}

static hdsdp_retcode conjGradSolve( iterative_linsys *cg, double *rhsVec, double *solVec ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    int nCGMaxIter = cg->params.maxIter;
    cg->useJacobi = cg->params.useJacobi;
    
    /* CG working scalers */
    double dDotMTimesd = 0.0;
    double resiDotPreInvResi = 0.0;
    double alpha = 0.0;
    double beta = 0.0;
    double resiNorm = 0.0;
    double rhsNorm = 0.0;
    
    double *iterResi = cg->iterResi;
    double *iterResiNew = cg->iterResiNew;
    double *iterDirection = cg->iterDirection;
    double *MTimesDirection = cg->MTimesDirection;
    double *preInvResi = cg->preInvResi;
    double *iterVec = cg->iterVec;
    
    double cgStartTime = HUtilGetTimeStamp();
    double cgDuration = 0.0;
    
    int iter = 0;
    
    /* Get restart frequency */
    int nRestartFreq = HDSDP_MAX(cg->params.nRestartFreq, 20);
    
    int incx = 1;
    double one = 1.0;
    double minusOne = -1.0;
    
    HDSDP_ZERO(iterVec, double, cg->nCol);
    HDSDP_MEMCPY(iterResi, rhsVec, double, cg->nCol);
    
    resiNorm = nrm2(&cg->nCol, iterResi, &incx);
    rhsNorm = nrm2(&cg->nCol, rhsVec, &incx);
    
    double cgTol = HDSDP_MIN(cg->params.absTol, rhsNorm * cg->params.relTol);
    cgTol = HDSDP_MAX(cgTol, 0.1 * cg->params.absTol);
    
    /* Initial value is already a qualified solution */
    if ( resiNorm < cgTol ) {
        cg->solStatus = ITERATIVE_STATUS_OK;
        goto exit_cleanup;
    }
    
    HDSDP_MEMCPY(iterDirection, iterResi, double, cg->nCol);
    
    /* Initial iteration preparation */
    conjGradApplyPreconditioner(cg, iterDirection);
    HDSDP_MEMCPY(preInvResi, iterDirection, double, cg->nCol);
    fds_symv(cg->nCol, 1.0, cg->fullMatElem, iterDirection, 0.0, MTimesDirection);
    
    for ( iter = 0; iter < nCGMaxIter; ++iter ) {
        
        resiDotPreInvResi = dot(&cg->nCol, preInvResi, &incx, iterResi, &incx);
        dDotMTimesd = dot(&cg->nCol, iterDirection, &incx, MTimesDirection, &incx);
        alpha = resiDotPreInvResi / dDotMTimesd;
        axpy(&cg->nCol, &alpha, iterDirection, &incx, iterVec, &incx);
        
        if ( iter % nRestartFreq == 5 && cg->useJacobi ) {
            /* Restart CG solver */
#ifdef HDSDP_CONJGRAD_DEBUG
            printf("Conjugate gradient restarts at iteration %d. \n", iter);
#endif
            HDSDP_MEMCPY(iterResi, rhsVec, double, cg->nCol);
            fds_symv(cg->nCol, 1.0, cg->fullMatElem, iterVec, 0.0, iterDirection);
            axpy(&cg->nCol, &minusOne, iterDirection, &incx, iterResi, &incx);
            HDSDP_MEMCPY(iterDirection, iterResi, double, cg->nCol);
            conjGradApplyPreconditioner(cg, iterDirection);
            fds_symv(cg->nCol, 1.0, cg->fullMatElem, iterDirection, 0.0, MTimesDirection);
            HDSDP_MEMCPY(preInvResi, iterResi, double, cg->nCol);
            conjGradApplyPreconditioner(cg, preInvResi);
            continue;
        }
        
        HDSDP_MEMCPY(iterResiNew, iterResi, double, cg->nCol);
        double minusAlpha = -alpha;
        axpy(&cg->nCol, &minusAlpha, MTimesDirection, &incx, iterResiNew, &incx);
        HDSDP_MEMCPY(preInvResi, iterResiNew, double, cg->nCol);
        conjGradApplyPreconditioner(cg, preInvResi);
        beta = dot(&cg->nCol, iterResiNew, &incx, preInvResi, &incx);
        beta = beta / resiDotPreInvResi;
        axpby(&cg->nCol, &one, preInvResi, &incx, &beta, iterDirection, &incx);
        fds_symv(cg->nCol, 1.0, cg->fullMatElem, iterDirection, 0.0, MTimesDirection);
        HDSDP_MEMCPY(iterResi, iterResiNew, double, cg->nCol);
        resiNorm = nrm2(&cg->nCol, iterResi, &incx);
        
        if ( resiNorm != resiNorm ) {
            cg->solStatus = ITERATIVE_STATUS_NUMERICAL;
            retcode = HDSDP_RETCODE_FAILED;
            goto exit_cleanup;
        }
        
        if ( iter > 20 && resiNorm > 0.01 * rhsNorm ) {
            cg->solStatus = ITERATIVE_STATUS_MAXITER;
            break;
        }
        
#ifdef HDSDP_CONJGRAD_DEBUG
         printf("Conjugate iteration %d. Resi: %10.3e \n", iter, resiNorm);
#endif
        
        if ( resiNorm < cgTol ) {
            cg->solStatus = ITERATIVE_STATUS_OK;
            break;
        }
    }
    
    if ( iter >= nCGMaxIter || cg->solStatus == ITERATIVE_STATUS_MAXITER ) {
        cg->solStatus = ITERATIVE_STATUS_MAXITER;
        if ( cg->useJacobi ) {
#ifdef HDSDP_CONJGRAD_DEBUG
        printf("Conjugate gradient switched to Cholesky pre-conditioner.\n");
#endif
            cg->useJacobi = 0;
            cg->params.useJacobi = 0;
            HDSDP_CALL(conjGradBuildPreconditioner(cg));
            HDSDP_CALL(conjGradSolve(cg, rhsVec, solVec));
        } else {
#ifdef HDSDP_CONJGRAD_DEBUG
        printf("Conjugate gradient failed. \n");
#endif
            cg->solStatus = ITERATIVE_STATUS_FAILED;
            retcode = HDSDP_RETCODE_FAILED;
        }
    }
    
    HDSDP_MEMCPY(solVec, iterVec, double, cg->nCol);

exit_cleanup:
    
    /* Collect solution time statistics */
    cgDuration = HUtilGetTimeStamp() - cgStartTime;
    cg->nSolves += 1;
    cg->solveTime += cgDuration;
    cg->nIters = iter;
    
    return retcode;
}

static hdsdp_retcode conjGradLinSolverSolveN( void *chol, int nRhs, double *rhsVec, double *solVec ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    iterative_linsys *cg = (iterative_linsys *) chol;
    
    if ( solVec ) {
        for ( int iCol = 0; iCol < nRhs; ++iCol ) {
            HDSDP_CALL(conjGradSolve(chol, rhsVec + iCol * cg->nCol, solVec + iCol * cg->nCol))
        }
    } else {
        /* Overwrite solution */
        for ( int iCol = 0; iCol < nRhs; ++iCol ) {
            HDSDP_MEMCPY(cg->rhsBuffer, rhsVec + iCol * cg->nCol, double, cg->nCol);
            HDSDP_CALL(conjGradSolve(chol, cg->rhsBuffer, rhsVec + iCol * cg->nCol));
        }
    }
    
exit_cleanup:
    return retcode;
}

static hdsdp_retcode conjGradLinSolverGetDiag( void *chol, double *dummy ) {
    
    (void) chol;
    (void) dummy;
    
    return HDSDP_RETCODE_FAILED;
}

static void conjGradLinSolverInvert( void *chol, double *dummy1, double *dummy2 ) {
    
    (void) chol;
    (void) dummy1;
    (void) dummy2;
    
    return;
}

static void conjGradLinSolverClear( void *chol ) {
    
    if ( !chol ) {
        return;
    }
    
    iterative_linsys *cg = (iterative_linsys *) chol;
    
    HDSDP_FREE(cg->iterResi);
    HDSDP_FREE(cg->iterResiNew);
    HDSDP_FREE(cg->iterDirection);
    HDSDP_FREE(cg->preInvResi);
    HDSDP_FREE(cg->MTimesDirection);
    HDSDP_FREE(cg->iterVec);
    HDSDP_FREE(cg->rhsBuffer);
    HDSDP_FREE(cg->JacobiPrecond);
    
    lapackLinSolverDestroy((void **) &cg->lap);
    
    return;
}

static void conjGradLinSolverDestroy( void **pchol ) {
    
    if ( !pchol ) {
        return;
    }
    
    conjGradLinSolverClear(*pchol);
    HDSDP_FREE(*pchol);
    
    return;
}

static hdsdp_retcode lapackIndefiniteLinSolverCreate( void **pchol, int nCol ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    HDSDP_NULLCHECK(pchol);
    lapack_indef_flinsys *indlap = NULL;
    HDSDP_INIT(indlap, lapack_indef_flinsys, 1);
    HDSDP_MEMCHECK(indlap);
    
    HDSDP_ZERO(indlap, lapack_indef_flinsys, 1);
    indlap->nCol = nCol;
    
    HDSDP_INIT(indlap->dFullMatElem, double, nCol * nCol);
    
    /* Allocate workspace */
    indlap->iLwork = nCol * 8;
    HDSDP_INIT(indlap->dWork, double, indlap->iLwork);
    HDSDP_MEMCHECK(indlap->dWork);
    HDSDP_INIT(indlap->iAuxiIPIV, int, nCol);
    HDSDP_MEMCHECK(indlap->iAuxiIPIV);
    
    *pchol = indlap;
    
exit_cleanup:
    return retcode;
}

static void lapackIndefiniteLinSolverSetParamsDummy( void *chol, void *param ) {
    
    (void) chol;
    (void) param;
    
    return;
}

static hdsdp_retcode lapackIndefiniteLinSolverSymbolic( void *chol, int *dummy1, int *dummy2 ) {
    
    (void) chol;
    (void) dummy1;
    (void) dummy2;
    
    return HDSDP_RETCODE_OK;
}

static hdsdp_retcode lapackIndefiniteLinSolverNumeric( void *chol, int *dummy1, int *dummy2, double *dFullMatrix ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    lapack_indef_flinsys *indlap = (lapack_indef_flinsys *) chol;
    
    HDSDP_MEMCPY(indlap->dFullMatElem, dFullMatrix, double, indlap->nCol * indlap->nCol);
    
    int info = LAPACK_RET_OK;
    char uplolow = LAPACK_UPLOW_LOW;
    
    /* Do blocked LDL */
    dsytrf(&uplolow, &indlap->nCol, indlap->dFullMatElem, &indlap->nCol,
           indlap->iAuxiIPIV, indlap->dWork, &indlap->iLwork, &info);
    
    if ( info != LAPACK_RET_OK ) {
        retcode = HDSDP_RETCODE_FAILED;
        goto exit_cleanup;
    }
    
exit_cleanup:
    return retcode;
}

static hdsdp_retcode lapackIndefiniteLinSolverPsdCheck( void *chol, int *dummy1, int *dummy2, double *dFullMatrix, int *isPsd ) {
    
    (void) chol;
    (void) dummy1;
    (void) dummy2;
    (void) dFullMatrix;
    (void) isPsd;
    
    assert( 0 );
    return HDSDP_RETCODE_FAILED;
}

static void lapackIndefiniteLinSolverForwardN( void *chol, int dummy1, double *dummy2, double *dummy3 ) {
    
    (void) chol;
    (void) dummy1;
    (void) dummy2;
    (void) dummy3;
    
    return;
}

static void lapackIndefiniteLinSolverBackwardN( void *chol, int dummy1, double *dummy2, double *dummy3 ) {
    
    (void) chol;
    (void) dummy1;
    (void) dummy2;
    (void) dummy3;
    
    return;
}

static hdsdp_retcode lapackIndefiniteLinSolverSolveN( void *chol, int nRhs, double *rhsVec, double *solVec ) {
    
    lapack_indef_flinsys *indlap = (lapack_indef_flinsys *) chol;
    
    int info = LAPACK_RET_OK;
    char uplolow = LAPACK_UPLOW_LOW;
    
    if ( solVec ) {
        HDSDP_MEMCPY(solVec, rhsVec, double, nRhs * indlap->nCol);
        dsytrs(&uplolow, &indlap->nCol, &nRhs, indlap->dFullMatElem,
               &indlap->nCol, indlap->iAuxiIPIV, solVec, &indlap->nCol, &info);
    } else {
        dsytrs(&uplolow, &indlap->nCol, &nRhs, indlap->dFullMatElem,
               &indlap->nCol, indlap->iAuxiIPIV, rhsVec, &indlap->nCol, &info);
    }
    
    assert( info == LAPACK_RET_OK );
    
    return HDSDP_RETCODE_OK;
}

static hdsdp_retcode lapackIndefiniteLinSolverGetDiag( void *chol, double *dummy ) {
    
    (void) chol;
    (void) dummy;
    
    return HDSDP_RETCODE_FAILED;
}

static void lapackIndefiniteLinSolverInvert( void *chol, double *dummy1, double *dummy2 ) {
    
    (void) chol;
    (void) dummy1;
    (void) dummy2;
    
    return;
}

static void lapackIndefiniteLinSolverClear( void *chol ) {
    
    if ( !chol ) {
        return;
    }
    
    lapack_indef_flinsys *indlap = (lapack_indef_flinsys *) chol;
    HDSDP_FREE(indlap->dFullMatElem);
    HDSDP_FREE(indlap->dWork);
    HDSDP_FREE(indlap->iAuxiIPIV);
    
    HDSDP_ZERO(indlap, lapack_indef_flinsys, 1);

    return;
}

static void lapackIndefiniteLinSolverDestroy( void **pchol ) {
    
    if ( !pchol ) {
        return;
    }
    
    lapackIndefiniteLinSolverClear(*pchol);
    HDSDP_FREE(*pchol);
    
    return;
}

static hdsdp_retcode HFpLinsysSwitchToIndefinite( hdsdp_linsys_fp *HLin ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    /* Switch to indefinite solver */
    assert( HLin->LinType == HDSDP_LINSYS_DENSE_ITERATIVE );
    
    iterative_linsys *cg = (iterative_linsys *) HLin->chol;
    double *dataMatElem = cg->fullMatElem;
    
    HLin->LinType = HDSDP_LINSYS_DENSE_INDEFINITE;
    HLin->cholDestroy(&HLin->chol);
    
    HLin->cholCreate = lapackIndefiniteLinSolverCreate;
    HLin->cholSetParam = lapackIndefiniteLinSolverSetParamsDummy;
    HLin->cholSymbolic = lapackIndefiniteLinSolverSymbolic;
    HLin->cholNumeric = lapackIndefiniteLinSolverNumeric;
    HLin->cholPsdCheck = lapackIndefiniteLinSolverPsdCheck;
    HLin->cholFSolve = lapackIndefiniteLinSolverForwardN;
    HLin->cholBSolve = lapackIndefiniteLinSolverBackwardN;
    HLin->cholSolve = lapackIndefiniteLinSolverSolveN;
    HLin->cholGetDiag = lapackIndefiniteLinSolverGetDiag;
    HLin->cholInvert = lapackIndefiniteLinSolverInvert;
    HLin->cholDestroy = lapackIndefiniteLinSolverDestroy;
    
    HDSDP_CALL(HLin->cholCreate(&HLin->chol, HLin->nCol));
    HDSDP_CALL(HLin->cholNumeric(HLin->chol, NULL, NULL, dataMatElem));
    
exit_cleanup:
    return retcode;
}

extern hdsdp_retcode HFpLinsysCreate( hdsdp_linsys_fp **pHLin, int nCol, linsys_type Ltype ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    HDSDP_NULLCHECK(pHLin);
    
    hdsdp_linsys_fp *HLinsys = NULL;
    HDSDP_INIT(HLinsys, hdsdp_linsys_fp, 1);
    HDSDP_MEMCHECK(HLinsys);
    
    HLinsys->LinType = Ltype;
    HLinsys->nCol = nCol;
    
    switch ( Ltype ) {
        case HDSDP_LINSYS_DENSE_DIRECT:
            HLinsys->cholCreate = lapackLinSolverCreate;
            HLinsys->cholSetParam = lapackLinSolverSetParamsDummy;
            HLinsys->cholSymbolic = lapackLinSolverSymbolic;
            HLinsys->cholNumeric = lapackLinSolverNumeric;
            HLinsys->cholPsdCheck = lapackLinSolverPsdCheck;
            HLinsys->cholFSolve = lapackLinSolverForwardN;
            HLinsys->cholBSolve = lapackLinSolverBackwardN;
            HLinsys->cholSolve = lapackLinSolverSolveN;
            HLinsys->cholGetDiag = lapackLinSolverGetDiag;
            HLinsys->cholInvert = lapackLinSolverInvert;
            HLinsys->cholDestroy = lapackLinSolverDestroy;
            break;
        case HDSDP_LINSYS_SPARSE_DIRECT:
#ifdef LINSYS_PARDISO
            HLinsys->cholCreate = pardisoLinSolverCreate;
            HLinsys->cholSetParam = pardisoLinSolverSetThreads;
            HLinsys->cholSymbolic = pardisoLinSolverSymbolic;
            HLinsys->cholNumeric = pardisoLinSolverNumeric;
            HLinsys->cholPsdCheck = pardisoLinSolverPsdCheck;
            HLinsys->cholFSolve = pardisoLinSolverForwardN;
            HLinsys->cholBSolve = pardisoLinSolverBackwardN;
            HLinsys->cholSolve = pardisoLinSolverSolveN;
            HLinsys->cholGetDiag = pardisoLinSolverGetDiag;
            HLinsys->cholInvert = pardisoLinSolverInvert;
            HLinsys->cholDestroy = pardisoLinSolverDestroy;
            break;
#else
            HLinsys->cholCreate = qdldlLinSolverCreate;
            HLinsys->cholSetParam = qdldlLinSolverSetThreads;
            HLinsys->cholSymbolic = qdldlLinSolverSymbolic;
            HLinsys->cholNumeric = qdldlLinSolverNumeric;
            HLinsys->cholPsdCheck = qdldlLinSolverPsdCheck;
            HLinsys->cholFSolve = qdldlLinSolverForwardN;
            HLinsys->cholBSolve = qdldlLinSolverBackwardN;
            HLinsys->cholSolve = qdldlLinSolverSolveN;
            HLinsys->cholGetDiag = qdldlLinSolverGetDiag;
            HLinsys->cholInvert = qdldlLinSolverInvert;
            HLinsys->cholDestroy = qdldlLinSolverDestroy;
            break;
#endif
        case HDSDP_LINSYS_SPARSE_INDEFINITE:
#ifdef LINSYS_PARDISO
            HLinsys->cholCreate = pardisoLinSolverCreateIndefinite;
            HLinsys->cholSetParam = pardisoLinSolverSetThreads;
            HLinsys->cholSymbolic = pardisoLinSolverSymbolic;
            HLinsys->cholNumeric = pardisoLinSolverStableNumeric;
            HLinsys->cholPsdCheck = pardisoLinSolverPsdCheck;
            HLinsys->cholFSolve = pardisoLinSolverForwardN;
            HLinsys->cholBSolve = pardisoLinSolverBackwardN;
            HLinsys->cholSolve = pardisoLinSolverSolveN;
            HLinsys->cholGetDiag = pardisoLinSolverGetDiag;
            HLinsys->cholInvert = pardisoLinSolverInvert;
            HLinsys->cholDestroy = pardisoLinSolverDestroy;
            break;
#else       
            HLinsys->cholCreate = ldlLinSolverCreate;
            HLinsys->cholSetParam = ldlLinSolverSetThreads;
            HLinsys->cholSymbolic = ldlLinSolverSymbolic;
            HLinsys->cholNumeric = ldlLinSolverNumeric;
            HLinsys->cholPsdCheck = ldlLinSolverPsdCheck;
            HLinsys->cholFSolve = ldlLinSolverForwardN;
            HLinsys->cholBSolve = ldlLinSolverBackwardN;
            HLinsys->cholSolve = ldlLinSolverSolveN;
            HLinsys->cholGetDiag = ldlLinSolverGetDiag;
            HLinsys->cholInvert = ldlLinSolverInvert;
            HLinsys->cholDestroy = ldlLinSolverDestroy;
            break;
#endif
        case HDSDP_LINSYS_DENSE_ITERATIVE:
            HLinsys->cholCreate = conjGradLinSolverCreate;
            HLinsys->cholSetParam = conjGradLinSolverSetParam;
            HLinsys->cholSymbolic = conjGradLinSolverSymbolic;
            HLinsys->cholNumeric = conjGradLinSolverNumeric;
            HLinsys->cholPsdCheck = conjGradLinSolverPsdCheck;
            HLinsys->cholFSolve = conjGradLinSolverForwardN;
            HLinsys->cholBSolve = conjGradLinSolverBackwardN;
            HLinsys->cholSolve = conjGradLinSolverSolveN;
            HLinsys->cholGetDiag = conjGradLinSolverGetDiag;
            HLinsys->cholInvert = conjGradLinSolverInvert;
            HLinsys->cholDestroy = conjGradLinSolverDestroy;
            break;
        case HDSDP_LINSYS_SPARSE_ITERATIVE:
#ifdef LINSYS_PARDISO
            HLinsys->cholCreate = pardisoLinSolverCreate;
            HLinsys->cholSetParam = pardisoLinSolverSetThreads;
            HLinsys->cholSymbolic = pardisoLinSolverSymbolic;
            HLinsys->cholNumeric = pardisoLinSolverNumeric;
            HLinsys->cholPsdCheck = pardisoLinSolverPsdCheck;
            HLinsys->cholFSolve = pardisoLinSolverForwardN;
            HLinsys->cholBSolve = pardisoLinSolverBackwardN;
            HLinsys->cholSolve = pardisoLinSolverSolveN;
            HLinsys->cholGetDiag = pardisoLinSolverGetDiag;
            HLinsys->cholInvert = pardisoLinSolverInvert;
            HLinsys->cholDestroy = pardisoLinSolverDestroy;
            break;
#else
            HLinsys->cholCreate = qdldlLinSolverCreate;
            HLinsys->cholSetParam = qdldlLinSolverSetThreads;
            HLinsys->cholSymbolic = qdldlLinSolverSymbolic;
            HLinsys->cholNumeric = qdldlLinSolverNumeric;
            HLinsys->cholPsdCheck = qdldlLinSolverPsdCheck;
            HLinsys->cholFSolve = qdldlLinSolverForwardN;
            HLinsys->cholBSolve = qdldlLinSolverBackwardN;
            HLinsys->cholSolve = qdldlLinSolverSolveN;
            HLinsys->cholGetDiag = qdldlLinSolverGetDiag;
            HLinsys->cholInvert = qdldlLinSolverInvert;
            HLinsys->cholDestroy = qdldlLinSolverDestroy;
            break;
#endif
        /* Not implemented */
        case HDSDP_LINSYS_SMALL_DIRECT:
            retcode = HDSDP_RETCODE_FAILED;
            goto exit_cleanup;
        default:
            retcode = HDSDP_RETCODE_FAILED;
            goto exit_cleanup;
    }
    
    *pHLin = HLinsys;
    
    HDSDP_CALL(HLinsys->cholCreate(&HLinsys->chol, nCol));
    
exit_cleanup:
    return retcode;
}

extern void HFpLinsysSetParam( hdsdp_linsys_fp *HLin, double relTol, double absTol, int nThreads, int maxIter, int nRestartFreq ) {
    
    if ( HLin->LinType == HDSDP_LINSYS_DENSE_ITERATIVE ) {
        
        iterative_params params;
        params.absTol = absTol;
        params.relTol = relTol;
        params.maxIter = maxIter;
        params.nRestartFreq = nRestartFreq;
        HLin->cholSetParam(HLin->chol, &params);
        
    } else if ( HLin->LinType == HDSDP_LINSYS_SPARSE_DIRECT ) {
        
        void *pnThreads = (void *) (&nThreads);
        HLin->cholSetParam(HLin->chol, pnThreads);
        
    }
    
    return;
}

extern hdsdp_retcode HFpLinsysSymbolic( hdsdp_linsys_fp *HLin, int *colMatBeg, int *colMatIdx ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    HDSDP_CALL(HLin->cholSymbolic(HLin->chol, colMatBeg, colMatIdx));
    
exit_cleanup:
    return retcode;
}

extern hdsdp_retcode HFpLinsysNumeric( hdsdp_linsys_fp *HLin, int *colMatBeg, int *colMatIdx, double *colMatElem ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    retcode = HLin->cholNumeric(HLin->chol, colMatBeg, colMatIdx, colMatElem);
    
    /* Cholesky fails. Switch to indefinite solver */
    if ( retcode == HDSDP_RETCODE_FAILED && HLin->LinType == HDSDP_LINSYS_DENSE_ITERATIVE ) {
        hdsdp_printf("KKT system is almost indefinite. Switch to LDL. \n");
        HDSDP_CALL(HFpLinsysSwitchToIndefinite(HLin));
    }
    
    HLin->nFactorizes += 1;
    
exit_cleanup:
    return retcode;
}

#ifdef LINSYS_PARDISO
extern hdsdp_retcode HFpLinsysSwitchToBackUp( hdsdp_linsys_fp *HLin ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    assert( HLin->LinType == HDSDP_LINSYS_SPARSE_INDEFINITE );
    HLin->cholNumeric = pardisoLinSolverStableNumeric;
    
exit_cleanup:
    return retcode;
}
#endif

extern hdsdp_retcode HFpLinsysPsdCheck( hdsdp_linsys_fp *HLin, int *colMatBeg, int *colMatIdx, double *colMatElem, int *isPsd ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    HDSDP_CALL(HLin->cholPsdCheck(HLin->chol, colMatBeg, colMatIdx, colMatElem, isPsd));
    
    HLin->nFactorizes += 1;
    
exit_cleanup:
    return retcode;
}

extern void HFpLinsysFSolve( hdsdp_linsys_fp *HLin, int nRhs, double *rhsVec, double *solVec ) {
    
    HLin->cholFSolve(HLin->chol, nRhs, rhsVec, solVec);
    return;
}

extern void HFpLinsysBSolve( hdsdp_linsys_fp *HLin, int nRhs, double *rhsVec, double *solVec ) {
    
    HLin->cholBSolve(HLin->chol, nRhs, rhsVec, solVec);
    
    return;
}

extern hdsdp_retcode HFpLinsysSolve( hdsdp_linsys_fp *HLin, int nRhs, double *rhsVec, double *solVec ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    retcode = HLin->cholSolve(HLin->chol, nRhs, rhsVec, solVec);
    
    /* Check NAN */
    if ( solVec ) {
        if ( solVec[0] != solVec[0] ) {
            retcode = HDSDP_RETCODE_FAILED;
        }
    }
    
    if ( rhsVec[0] != rhsVec[0] ) {
        retcode = HDSDP_RETCODE_FAILED;
    }
    
    if ( retcode != HDSDP_RETCODE_OK && HLin->LinType == HDSDP_LINSYS_DENSE_ITERATIVE ) {
        hdsdp_printf("KKT system is unstable. Switch to LDL. \n");
        HDSDP_CALL(HFpLinsysSwitchToIndefinite(HLin));
        HDSDP_CALL(HFpLinsysSolve(HLin, nRhs, rhsVec, solVec));
    }
    
    HLin->nSolves += 1;
    
exit_cleanup:
    
    return retcode;
}

extern hdsdp_retcode HFpLinsysGetDiag( hdsdp_linsys_fp *HLin, double *diagElem ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    HLin->cholGetDiag(HLin->chol, diagElem);
       
    return retcode;
}

extern void HFpLinsysInvert( hdsdp_linsys_fp *HLin, double *dFullMatrix, double *dAuxiMatrix ) {
    
    HLin->cholInvert(HLin->chol, dFullMatrix, dAuxiMatrix);
    
    return;
}

extern void HFpLinsysClear( hdsdp_linsys_fp *HLin ) {
    
    if ( !HLin ) {
        return;
    }
    
    HLin->nCol = 0;
    HLin->cholDestroy(&HLin->chol);
    
    return;
}

extern void HFpLinsysDestroy( hdsdp_linsys_fp **HLin ) {
    
    if ( !HLin ) {
        return;
    }
    
    HFpLinsysClear(*HLin);
    HDSDP_FREE(*HLin);
    
#ifdef HDSDP_LINSYS_PROFILE
    
    if ( iPrintPardiso ) {
        hdsdp_printf("\n-------------------------------\n");
        hdsdp_printf("Sparse Linear system statistics \n");
        hdsdp_printf("-------------------------------\n");
        hdsdp_printf("\nPardiso time: \nNumeric: %f \nSolve: %f \nCheck: %f \nInvert %f \nForward: %f \nBackward: %f \n",
               dPardisoNumeric, dPardisoSolve, dPardisoPsdCheck, dPardisoInvert, dPardisoForward, dPardisoBackward);
        hdsdp_printf("\nPardiso call: \nNumeric: %d \nSolve: %d \nCheck: %d \nInvert %d \nForward: %d \nBackward: %d \n",
               nPardisoNumeric, nPardisoSolve, nPardisoPsdCheck, nPardisoInvert, nPardisoForward, nPardisoBackward);
        iPrintPardiso = 0;
    }
    
    dPardisoSolve = 0.0;
    dPardisoNumeric = 0.0;
    dPardisoPsdCheck = 0.0;
    dPardisoInvert = 0.0;
    dPardisoForward = 0.0;
    dPardisoBackward = 0.0;
    
    nPardisoSolve = -1;
    nPardisoNumeric = -1;
    nPardisoPsdCheck = -1;
    nPardisoInvert = -1;
    nPardisoForward = -1;
    nPardisoBackward = -1;
    
    if ( iPrintLapack ) {
        hdsdp_printf("\n------------------------------\n");
        hdsdp_printf("Dense Linear system statistics  \n");
        hdsdp_printf("------------------------------\n");
        hdsdp_printf("\nLapack time: \nNumeric: %f \nSolve: %f \nCheck: %f \nInvert %f \nForward: %f \nBackward: %f \n",
               dLapackNumeric, dLapackSolve, dLapackPsdCheck, dLapackInvert, dLapackForward, dLapackBackward);
        hdsdp_printf("\nLapack call: \nNumeric: %d \nSolve: %d \nCheck: %d \nInvert %d \nForward: %d \nBackward: %d \n",
               nLapackNumeric, nLapackSolve, nLapackPsdCheck, nLapackInvert, nLapackForward, nLapackBackward);
        iPrintLapack = 0;
    }
    
    dLapackSolve = 0.0;
    dLapackNumeric = 0.0;
    dLapackPsdCheck = 0.0;
    dLapackInvert = 0.0;
    dLapackForward = 0.0;
    dLapackBackward = 0.0;
    
    nLapackSolve = -1;
    nLapackNumeric = -1;
    nLapackPsdCheck = -1;
    nLapackInvert = -1;
    nLapackForward = -1;
    nLapackBackward = -1;
    
#endif
    
    return;
}
