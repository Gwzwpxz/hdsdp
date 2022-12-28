/** @file sgm.c
 *  @brief Implement an LP solver based on sharpness and subgradient
 *
 * To allow the routine's use in external code, the SGM solver is written as a driver
 */

#include "sgm.h"
#include "pot_utils.h"
#include "vec_mat.h"

#include <math.h>

static int potSGMIGetUbID( int nCol, double *colUBound, int *ubId ) {
    
    int nUb = 0;
    for ( int i = 0; i < nCol; ++i ) {
        if ( colUBound[i] < SGM_INFINITY ) {
            ubId[nUb] = i;
            nUb += 1;
        }
    }
    
    return nUb;
}

/** @brief Subgradient solver for LP in the potential solver
 *
 */
extern int potSGMSolve( int nCol, int nEqRow, int nIneqRow, int *eqMatBeg, int *eqMatIdx, double *eqMatElem,
                        int *inEqMatBeg, int *inEqMatIdx, double *inEqMatElem, double *eqMatRhs, double *ineqMatRhs,
                        double *colObj, double *colBoundUp, double *colVal, double *rowDual, double *colSlack,
                        double *bdSlack, int maxIter, double relTol, double stepScal ) {
    
    int retcode = SGM_RETCODE_OK;
    
    /* Rename data */
    int meq = nEqRow, mineq = nIneqRow, n = nCol;
    int mall = meq + mineq;
    int *Ap = eqMatBeg, *Ai = eqMatIdx;
    double *Ax = eqMatElem;
    int *Gp = inEqMatBeg, *Gi = inEqMatIdx;
    double *Gx = inEqMatElem;
    double *ub = colBoundUp, *b = eqMatRhs, *g = ineqMatRhs, *c = colObj;
    
    /* Prepare iterations */
    double *x = NULL, *xbest = NULL, *y = NULL, *ybest = NULL, *w = NULL, *wbest = NULL;
    double *dx = NULL, *dy = NULL, *dw = NULL;
    double *Kx = NULL, *sKx = NULL, *KTy = NULL, *sKTy = NULL;
    
    /* Prepare statistics */
    int *ubid = NULL;
    double pres = 0.0, dres = 0.0, cres = 0.0, cpl = 0.0;
    double nrmdx = 0.0, nrmdy = 0.0, nrmdw = 0.0;
    double fval = 0.0, fbest = SGM_INFINITY;
    
    double tol = relTol;
    double aa = stepScal;
    
    /* Allocate working space */
    POTLP_INIT(x, double, n);
    POTLP_INIT(xbest, double, n);
    POTLP_INIT(dx, double, n);
    POTLP_INIT(y, double, mall);
    POTLP_INIT(dy, double, mall);
    POTLP_INIT(ybest, double, mall);
    
    /* K = [A; G] */
    POTLP_INIT(Kx, double, mall);
    POTLP_INIT(sKx, double, mall);
    POTLP_INIT(KTy, double, n);
    POTLP_INIT(sKTy, double, n);
    
    if ( !x || !xbest || !dx || !y || !dy || !ybest || !Kx || !sKx || !KTy || !sKTy ) {
        retcode = SGM_RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    /* Special treatment for w */
    POTLP_INIT(ubid, int, n);
    
    if ( !ubid ) {
        retcode = SGM_RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    int nub = potSGMIGetUbID(n, ub, ubid);
    
    if ( nub > 0 ) {
        POTLP_INIT(w, double, nub);
        POTLP_INIT(wbest, double, nub);
        POTLP_INIT(dw, double, nub);
        if ( !w || !wbest || !dw ) {
            retcode = SGM_RETCODE_FAILED;
            goto exit_cleanup;
        }
    }
    
    /* Initialize */
    POTLP_MEMCPY(x, colVal, double, n);
    POTLP_MEMCPY(xbest, x, double, n);
    POTLP_MEMCPY(y, rowDual, double, meq + mineq);
    POTLP_MEMCPY(ybest, y, double, meq + mineq);
    
    double *yeq = y, *yineq = y + meq;
    double *dyeq = dy, *dyineq = dy + meq;
    double *Kxeq = Kx, *Kxineq = Kx + meq;
    
    /* Start iterating */
    int iter = 0;
    for ( iter = 0; iter < maxIter && fval >= tol; ++iter ) {
        
        /* Kx = K * x - q */
        POTLP_MEMCPY(Kxeq, b, double, meq);
        POTLP_MEMCPY(Kxineq, g, double, mineq);
        spMatAxpy(n, Ap, Ai, Ax, -1.0, x, Kxeq);
        spMatAxpy(n, Gp, Gi, Gx, -1.0, x, Kxineq);
        
        /* Kx(ineqidx) = min(Kx(ineqidx), 0); */
        for ( int i = 0; i < mineq; ++i ) {
            Kxineq[i] = POTLP_MAX(Kxineq[i], 0.0);
        }
        
        for ( int i = 0; i < mall; ++i ) {
            sKx[i] = ( sKx[i] > 0.0 ) ? 1.0 : -1.0;
        }
        
        /* dx = -K' * sign(Kx); */
        POTLP_ZERO(dx, double, n);
        spMatATxpy(n, Ap, Ai, Ax, -1.0, sKx, dx);
        
        /* KTy = c - K' * y; */
        POTLP_MEMCPY(KTy, c, double, n);
        spMatATxpy(n, Ap, Ai, Ax, -1.0, yeq, KTy);
        spMatATxpy(n, Gp, Gi, Gx, -1.0, yineq, KTy);
        
        /* KTy(ubid) = KTy(ubid) + w; */
        for ( int i = 0; i < nub; ++i ) {
            KTy[ubid[i]] += w[i];
        }
        
        /* KTy = min(KTy, 0);
           sKTy = sign(KTy);
         */
        for ( int i = 0; i < n; ++i ) {
            KTy[i] = POTLP_MIN(KTy[i], 0.0);
            sKTy[i] = ( KTy[i] > 0.0 ) ? 1.0 : 0.0;
        }
        
        /* dy = -K * sKTy;
           dw = sKTy(ubid);
         */
        
        POTLP_ZERO(dy, double, mall);
        spMatAxpy(n, Ap, Ai, Ax, -1.0, sKTy, dyeq);
        spMatAxpy(n, Gp, Gi, Gx, -1.0, sKTy, dyineq);
        
        for ( int i = 0; i < nub; ++i ) {
            dw[i] = sKTy[ubid[i]];
        }
    
        /*
         % Compl.
             cpl = c' * x - q' * y + u' * w;
             dx = pr * dx + pdr * c * sign(cpl);
             dy = dr * dy - pdr * q * sign(cpl);
             dw = dr * dw + pdr * u * sign(cpl);
         */
        double cTx = dot(&n, c, &potIntConstantOne, x, &potIntConstantOne);
        double bTyeq = dot(&meq, b, &potIntConstantOne, yeq, &potIntConstantOne);
        double gTyineq = dot(&mineq, g, &potIntConstantOne, yineq, &potIntConstantOne);
        double uTw = 0.0;
        
        for ( int i = 0; i < nub; ++i ) {
            uTw += w[i] * ub[ubid[i]];
        }
        
        cpl = cTx - bTyeq - gTyineq + uTw;
        
        if ( cpl > 0 ) {
            axpy(&n, &potDblConstantOne, c, &potIntConstantOne, dx, &potIntConstantOne);
            axpy(&meq, &potDblConstantMinusOne, b, &potIntConstantOne, dyeq, &potIntConstantOne);
            axpy(&mineq, &potDblConstantMinusOne, g, &potIntConstantOne, dyineq, &potIntConstantOne);
            for ( int i = 0; i < nub; ++i ) {
                dw[i] += ub[ubid[i]];
            }
            
        } else {
            axpy(&n, &potDblConstantMinusOne, c, &potIntConstantOne, dx, &potIntConstantOne);
            axpy(&meq, &potDblConstantOne, b, &potIntConstantOne, dyeq, &potIntConstantOne);
            axpy(&mineq, &potDblConstantOne, g, &potIntConstantOne, dyineq, &potIntConstantOne);
            for ( int i = 0; i < nub; ++i ) {
                dw[i] -= ub[ubid[i]];
            }
        }
        
        /*
         pres = sum(abs(Kx));
         dres = sum(abs(KTy));
         cres = abs(cpl);
         fval = pr * pres + dr * dres + pdr * cres;
         */
        
        pres = dres = cpl = 0.0;
        for ( int i = 0; i < mall; ++i ) {
            pres += fabs(Kx[i]);
        }
        
        for ( int i = 0; i < n; ++i ) {
            dres += fabs(KTy[i]);
        }
        
        cres = fabs(cpl);
        fval = pres + dres + cres;
        
        if ( fval < 1e-08 ) {
            break;
        }
        
        if ( fval < fbest ) {
            fbest = fval;
            POTLP_MEMCPY(xbest, x, double, n);
            POTLP_MEMCPY(ybest, y, double, mall);
            if ( nub > 0 ) {
                POTLP_MEMCPY(wbest, w, double, nub);
            }
        }
        
        /* Polyak */
        nrmdx = nrm2(&n, dx, &potIntConstantOne);
        nrmdy = nrm2(&mall, dy, &potIntConstantOne);
        nrmdw = nrm2(&nub, dw, &potIntConstantOne);
        
        double alpha = aa * fval;
        alpha = POTLP_MAX(alpha, 1.5 * fbest);
        alpha /= - ( nrmdx * nrmdx + nrmdy * nrmdy + nrmdw * nrmdw );
        
        /* Update
             x = x - alpha * dx;
             y = y - alpha * dy;
             w = w - alpha * dw;
        */
        axpy(&n, &alpha, dx, &potIntConstantOne, x, &potIntConstantOne);
        axpy(&mall, &alpha, dy, &potIntConstantOne, y, &potIntConstantOne);
        axpy(&nub, &alpha, dw, &potIntConstantOne, w, &potIntConstantOne);
        
        /* Projection */
        for ( int i = 0; i < n; ++i ) {
            x[i] = POTLP_MAX(x[i], 0.0);
        }
        
        for ( int i = 0; i < nub; ++i ) {
            x[ubid[i]] = POTLP_MIN(x[ubid[i]], ub[ubid[i]]);
        }
        
        for ( int i = 0; i < mineq; ++i ) {
            yineq[i] = POTLP_MAX(yineq[i], 0.0);
        }
        
        for ( int i = 0; i < nub; ++i ) {
            w[i] = POTLP_MAX(w[i], 0.0);
        }
        
        printf("%4d %5.2e %5.2e %5.2e | %5.2e \n", iter, pres, dres, cres, fval);;
        
        if ( iter % 1000 == 1 ) {
            aa *= 0.95;
        }
    }
    
    /* Successful exit */
    POTLP_MEMCPY(colVal, xbest, double, n);
    POTLP_MEMCPY(rowDual, ybest, double, mall);
    POTLP_ZERO(bdSlack, double, n);
    for ( int i = 0; i < nub; ++i ) {
        bdSlack[ubid[i]] = wbest[i];
    }
    
    POTLP_MEMCPY(colSlack, c, double, n);
    axpy(&n, &potDblConstantOne, bdSlack, &potIntConstantOne, colSlack, &potIntConstantOne);
    spMatATxpy(n, Ap, Ai, Ax, -1.0, yeq, colSlack);
    spMatATxpy(n, Gp, Gi, Gx, -1.0, yineq, colSlack);
    
    /* Done */
    
exit_cleanup:
    
    POTLP_FREE(x);
    POTLP_FREE(xbest);
    POTLP_FREE(dx);
    POTLP_FREE(y);
    POTLP_FREE(dy);
    POTLP_FREE(ybest);
    
    POTLP_FREE(Kx);
    POTLP_FREE(sKx);
    POTLP_FREE(KTy);
    POTLP_FREE(sKTy);
    
    POTLP_FREE(ubid);
    POTLP_FREE(w);
    POTLP_FREE(wbest);
    POTLP_FREE(dw);
    
    return retcode;
}
