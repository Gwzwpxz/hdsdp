/** @file qp\_avg.h
 * QP solver for averaging acceleration from potential reduction
 * Cache efficient and not aware of sparsity
 */

#include "pot_def.h"
#include "pot_utils.h"
#include "pot_solver.h"
#include "qp_avg.h"
#include "vec_mat.h"

#include <math.h>

static void potQPIInit( pot_qpsolver *potQP ) {
    
    /* Initialize QP solver with s = lbd = [1,...,1], alpha = [1/k,...,1/k], nu = 0 */
    for ( int i = 0; i < potQP->nRow; ++i ) {
        potQP->s[i] = 1.0;
        potQP->lbd[i] = 1.0;
    }
    
    for ( int i = 0; i < potQP->nQuadCol; ++i ) {
        potQP->alpha[i] = 1.0 / potQP->nQuadCol;
    }
        
    return;
}

static double potQPIEvalPObj( pot_qpsolver *potQP ) {
    
    return quadform(&potQP->nQuadCol, potQP->QMatElem, potQP->alpha);
}

static double potQPIGetMu( pot_qpsolver *potQP ) {
    
    double mu = dot(&potQP->nRow, potQP->s, &potIntConstantOne,
                    potQP->lbd, &potIntConstantOne);
    
    return (mu / potQP->nRow);
}

static double potQPIRatioTest( pot_int nQPRow, double *s, double *ds, double *lbd, double *dlbd ) {
    
    double stepTemp = POTLP_INFINITY;
    double ratio;
    
    for ( int i = 0; i < nQPRow; ++i ) {
        ratio = ds[i] / s[i];
        stepTemp = POTLP_MIN(ratio, stepTemp);
    }
    
    for ( int i = 0; i < nQPRow; ++i ) {
        ratio = dlbd[i] / lbd[i];
        stepTemp = POTLP_MIN(ratio, stepTemp);
    }
    
    return fabs(1.0 / stepTemp);
}

#ifndef POTQP_DEBUG
#define POTQP_DEBUG printf
#endif
/* Solve QP problem using infeasible start primal-dual interior point method */
static pot_int potQPISolve( pot_qpsolver *potQP, double pObjTarget ) {
    
    pot_int retcode = RETCODE_OK;
    
    pot_int nFreeCol = potQP->nQuadCol, nQPRow = potQP->nRow;
    
    potQPIInit(potQP);
    
    double *Q = potQP->QMatElem;
    double *V = potQP->colMatElem;
    double *s = potQP->s;
    double *ds = potQP->ds;
    double *lbd = potQP->lbd;
    double *dlbd = potQP->dlbd;
    double *alpha = potQP->alpha;
    double *dalpha = potQP->dalpha;
    double *alphabest = potQP->alphabest;
    double *M = potQP->M;
    double *SLV = potQP->SLV;
    double *d1 = potQP->d1;
    double *d2 = potQP->d2;
    double *Valpha = potQP->Valpha;
    double *buffer = potQP->buffer;
    
    double eTd1, eTd2, eTalpha, dnu, step, musigma;
    pot_int nQuadElem = nFreeCol * nFreeCol;
    pot_int nColMatElem = nFreeCol * nQPRow;
    double mu = 1.0, sigma = 0.1, nu = 0.0;
    double pObjVal = POTLP_INFINITY;
    double pObjBest = pObjVal;
    
    int iter, info;
    
#ifdef POTQP_DEBUG
    POTQP_DEBUG("-------------------------------------"
                "-------------------------------------\n");
#endif
    
    /* Start the interior point iteration */
    for ( iter = 0; iter < 100 && mu > 1e-10; ++iter ) {
        
        mu = potQPIGetMu(potQP);
        musigma = mu * sigma;
        
        POTLP_MEMCPY(SLV, potQP->colMatElem, double, nColMatElem);
        
        for ( int i = 0; i < nQPRow; ++i ) {
            buffer[i] = sqrtl(lbd[i]) / sqrtl(s[i]);
        }
        
        /* SLV = diag(sqrt(sinvl)) * V; */
        for ( int i = 0; i < nFreeCol; ++i ) {
            vvscl(&nQPRow, buffer, SLV + i * nQPRow);
        }
        
        /* VTSLV = SLV' * SLV; */
        syrk(&potCharConstantLow, &potCharConstantTrans, &nFreeCol, &nQPRow,
             &potDblConstantOne, SLV, &nQPRow, &potDblConstantZero, M, &nFreeCol);
        
        for ( int i = 0; i < nQPRow; ++i ) {
            ds[i] = musigma / s[i];
        }
        
        /* M = Q + VTSLV + eye(k) * 1e-15; */
        axpy(&nQuadElem, &potDblConstantOne, Q, &potIntConstantOne, M, &potIntConstantOne);
        
        /* r = M * alpha + e * nu - V' * (lbd + mu * sigma * sinv); */
        for ( int i = 0; i < nFreeCol; ++i ) {
            d1[i] = 1.0; d2[i] = nu;
        }
        
        for ( int i = 0; i < nQPRow; ++i ) {
            Valpha[i] = -lbd[i] - ds[i];
        }
        
        symv(&potCharConstantLow, &nFreeCol, &potDblConstantOne, M, &nFreeCol,
             alpha, &potIntConstantOne, &potDblConstantOne, d2, &potIntConstantOne);
        gemv(&potCharConstantTrans, &nQPRow, &nFreeCol, &potDblConstantOne,
             V, &nQPRow, Valpha, &potIntConstantOne, &potDblConstantOne, d2, &potIntConstantOne);
        
        /* Cholesky. Some perturbation is needed */
        for ( int i = 0; i < nFreeCol; ++i ) {
            M[i * nFreeCol + i] += 1e-15;
        }
        
        potrf(&potCharConstantLow, &nFreeCol, M, &nFreeCol, &info);
        if ( info ) {
            retcode = RETCODE_FAILED;
            goto exit_cleanup;
        }
        
        /* Solve */
        potrs(&potCharConstantLow, &nFreeCol, &potIntConstantOne,
              M, &nFreeCol, d1, &nFreeCol, &info);
        if ( info ) {
            retcode = RETCODE_FAILED;
            goto exit_cleanup;
        }
        
        potrs(&potCharConstantLow, &nFreeCol, &potIntConstantOne,
              M, &nFreeCol, d2, &nFreeCol, &info);
        if ( info ) {
            retcode = RETCODE_FAILED;
            goto exit_cleanup;
        }
        
        /* dnu = -(e' * d2) / (e' * d1); */
        eTd1 = eTd2 = eTalpha = 0.0;
        for ( int i = 0; i < nFreeCol; ++i ) {
            eTd1 += d1[i]; eTd2 += d2[i]; eTalpha += alpha[i];
        }
        dnu = (eTalpha - 1.0 - eTd2) / eTd1;
        
        /* dalpha = -d1 * dnu - d2; */
        for ( int i = 0; i < nFreeCol; ++i ) {
            dalpha[i] = -d1[i] * dnu - d2[i];
        }
        
        /* Vaplusda = V * (alpha + dalpha); */
        POTLP_MEMCPY(d2, dalpha, double, nFreeCol);
        axpy(&nFreeCol, &potDblConstantOne, alpha, &potIntConstantOne, d2, &potIntConstantOne);
        gemv(&potCharConstantNoTrans, &nQPRow, &nFreeCol, &potDblConstantOne, V, &nQPRow,
             d2, &potIntConstantOne, &potDblConstantZero, Valpha, &potIntConstantOne);
        
        /* dlbd = mu * sigma * sinv - sinvl .* Vaplusda; */
        for ( int i = 0; i < nQPRow; ++i ) {
            dlbd[i] = ds[i] - buffer[i] * buffer[i] * Valpha[i];
        }
        
        /* ds = Vaplusda - s; */
        for ( int i = 0; i < nQPRow; ++i ) {
            ds[i] = Valpha[i] - s[i];
        }
        
        /* Ratio test */
        step = potQPIRatioTest(nQPRow, s, ds, lbd, dlbd);
        step = step * 0.995;
        step = POTLP_MIN(step, 1.0);
        
        /* Update variables */
        axpy(&nQPRow, &step, ds, &potIntConstantOne, s, &potIntConstantOne);
        axpy(&nQPRow, &step, dlbd, &potIntConstantOne, lbd, &potIntConstantOne);
        axpy(&nFreeCol, &step, dalpha, &potIntConstantOne, alpha, &potIntConstantOne);
        nu = nu + step * dnu;
        
        pObjVal = potQPIEvalPObj(potQP);
        
        if ( pObjVal < pObjBest ) {
            pObjBest = pObjVal;
            POTLP_MEMCPY(alphabest, alpha, double, nFreeCol);
        }
        
#ifdef POTQP_DEBUG
        if ( iter % 2 == 0 ) {
            gemv(&potCharConstantNoTrans, &nQPRow, &nFreeCol, &potDblConstantOne, V, &nQPRow,
                 alpha, &potIntConstantOne, &potDblConstantZero, buffer, &potIntConstantOne);
            for ( int i = 0; i < nQPRow; ++i ) {
                buffer[i] -= s[i];
            }
            double pInf = nrm2(&nQPRow, buffer, &potIntConstantOne);
            gemv(&potCharConstantTrans, &nQPRow, &nFreeCol, &potDblConstantMinusOne, V, &nQPRow,
                 lbd, &potIntConstantOne, &potDblConstantZero, buffer, &potIntConstantOne);
            symv(&potCharConstantLow, &nFreeCol, &potDblConstantOne, potQP->QMatElem, &nFreeCol,
                 alpha, &potIntConstantOne, &potDblConstantOne, buffer, &potIntConstantOne);
            for ( int i = 0; i < nFreeCol; ++i ) {
                buffer[i] += nu;
            }
            double dInf = nrm2(&nFreeCol, buffer, &potIntConstantOne);
            printf("  %3d %10.6e %10.3e %10.3e %10.3e \n", iter + 1, pObjVal, pInf, dInf, mu);
        }
#endif
    }
    
#ifdef POTQP_DEBUG
    POTQP_DEBUG("   [QP] Obj. Reduction rate %f \n", pObjBest / pObjTarget);
    POTQP_DEBUG("-------------------------------------"
                "-------------------------------------\n");
#endif
    
    /* Subspace does not contain a better point */
    if ( pObjBest > pObjTarget * 0.999 ) {
        retcode = RETCODE_FAILED;
    }
    
exit_cleanup:
    
    return retcode;
}

extern pot_int potQPCreate( pot_qpsolver **ppotQP ) {
    
    pot_int retcode = RETCODE_OK;
    
    if ( !ppotQP ) {
        retcode = RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    pot_qpsolver *potQP = NULL;
    POTLP_INIT(potQP, pot_qpsolver, 1);
    
    if ( !potQP ) {
        retcode = RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    POTLP_ZERO(potQP, pot_qpsolver, 1);
    *ppotQP = potQP;
    
exit_cleanup:
    
    return retcode;
}

extern pot_int potQPInit( pot_qpsolver *potQP, pot_int nQuadCol, pot_int nRowAll, pot_int nRow ) {
    
    pot_int retcode = RETCODE_OK;
    
    if ( !potQP ) {
        retcode = RETCODE_OK;
        goto exit_cleanup;
    }
    
    if ( nRow < nQuadCol ) {
        retcode = RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    potQP->nQuadCol = nQuadCol;
    potQP->nRow = nRow;
    potQP->nRowAll = nRowAll;
    
    POTLP_INIT(potQP->QMatElem, double, nQuadCol * nQuadCol);
    POTLP_INIT(potQP->colMatElem, double, nQuadCol * nRow);
    POTLP_INIT(potQP->s, double, nRow);
    POTLP_INIT(potQP->ds, double, nRow);
    POTLP_INIT(potQP->lbd, double, nRow);
    POTLP_INIT(potQP->dlbd, double, nRow);
    POTLP_INIT(potQP->alpha, double, nQuadCol);
    POTLP_INIT(potQP->dalpha, double, nQuadCol);
    POTLP_INIT(potQP->alphabest, double, nQuadCol);
    POTLP_INIT(potQP->M,  double, nQuadCol * nQuadCol);
    POTLP_INIT(potQP->SLV, double, nRow * nQuadCol);
    POTLP_INIT(potQP->d1, double, nQuadCol);
    POTLP_INIT(potQP->d2, double, nQuadCol);
    POTLP_INIT(potQP->Valpha, double, nRow);
    POTLP_INIT(potQP->buffer, double, nRow);
    
    if ( !potQP->QMatElem || !potQP->colMatElem || !potQP->s || !potQP->ds || !potQP->lbd ||
         !potQP->dlbd || !potQP->alpha || !potQP->dalpha || !potQP->alphabest || !potQP->M ||
         !potQP->SLV || !potQP->d1 || !potQP->d2 || !potQP->Valpha || !potQP->buffer ) {
        retcode = RETCODE_FAILED;
        goto exit_cleanup;
    }
    
exit_cleanup:
    
    return retcode;
}

/** @brief Load QP into the IPM solver
 *
 */
extern void potQPLoadProb( pot_qpsolver *potQP, int nIter, pot_int nCol, pot_int nRow, pot_int nCone, double *gradWindow, double *xVarWindow ) {
    
    /* Set up Q = X' * A' * A * X */
    double *QMatElem = potQP->QMatElem;
    syrk(&potCharConstantLow, &potCharConstantTrans, &potQP->nQuadCol, &nRow, &potDblConstantOne,
         gradWindow, &nRow, &potDblConstantZero, QMatElem, &potQP->nQuadCol);
    
    /* Copy Xc = X(coneidx, :) */
    double *Xdst = potQP->colMatElem;
    double *Xsrc = xVarWindow + nCol - nCone;
        
    for ( int i = 0; i < potQP->nQuadCol; ++i ) {
        POTLP_MEMCPY(Xdst, Xsrc, double, nCone);
        Xdst += nCone; Xsrc += nCol;
    }
    
    potQP->xWindow = xVarWindow;
    
    return;
}

/** @brief Solve the averaging QP
 *
 */
extern pot_int potQPSolveProb( pot_qpsolver *potQP, double targetObj, double *xAvg ) {
    
    pot_int retcode = RETCODE_OK;
    POT_CALL(potQPISolve(potQP, targetObj));
    
    gemv(&potCharConstantNoTrans, &potQP->nRowAll, &potQP->nQuadCol, &potDblConstantOne,
         potQP->xWindow, &potQP->nRowAll, potQP->alphabest, &potIntConstantOne,
         &potDblConstantZero, xAvg, &potIntConstantOne);
    
exit_cleanup:
    
    return retcode;
}

extern void potQPClear( pot_qpsolver *potQP ) {
    
    if ( !potQP ) {
        return;
    }
    
    POTLP_FREE(potQP->QMatElem);
    POTLP_FREE(potQP->colMatElem);
    POTLP_FREE(potQP->s);
    POTLP_FREE(potQP->ds);
    POTLP_FREE(potQP->lbd);
    POTLP_FREE(potQP->dlbd);
    POTLP_FREE(potQP->alpha);
    POTLP_FREE(potQP->dalpha);
    POTLP_FREE(potQP->alphabest);
    POTLP_FREE(potQP->M);
    POTLP_FREE(potQP->SLV);
    POTLP_FREE(potQP->d1);
    POTLP_FREE(potQP->d2);
    POTLP_FREE(potQP->Valpha);
    POTLP_FREE(potQP->buffer);
    
    POTLP_ZERO(potQP, pot_qpsolver, 1);
    
    return;
}

extern void potQPDestroy( pot_qpsolver **ppotQP ) {
    
    if ( !ppotQP ) {
        return;
    }
    
    potQPClear(*ppotQP);
    POTLP_FREE(*ppotQP);
    
    return;
}
