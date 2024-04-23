/** @file hdsdp\_lpsolve.c
 *  @brief Specialized LP solver
 *
 */

#ifdef HEADERPATH
#include "interface/hdsdp_utils.h"
#include "interface/hdsdp_linsolver.h"
#include "interface/hdsdp_lpsolve.h"
#include "linalg/vec_opts.h"
#else
#include "hdsdp_utils.h"
#include "hdsdp_linsolver.h"
#include "hdsdp_lpsolve.h"
#include "vec_opts.h"
#endif

#include <math.h>

#define PREGULARIZER (1.8190e-12)
#define DREGULARIZER (1.4901e-08)

/* Corrector types */
#define NOCORR_HSD        (-1)
#define PDIPM_HSD         ( 0)
#define MEHROTRA_HSD      ( 1)
#define MEHROTRA          ( 2)
#define PRIMAL            ( 3)

#define HDSDP_DIV(x, y) ((x) / (y))

static double LpNewtonIBarrier( int nCol, double *x, double *s, double kappa, double tau ) {
    
    double mu = kappa * tau;
    for ( int i = 0; i < nCol; ++i ) {
        mu += x[i] * s[i];
    }
    
    return mu / (nCol + 1);
}

static double LpNewtonIRatioTest( int nCol, double *x, double *dx, double *s, double *ds,
                                  double kappa, double dkappa, double tau, double dtau ) {
    
    /* 1 / abs(min([dx./x; ds./s; dkappa / kappa; dtau./tau])) */
    double alphaTmp = HDSDP_INFINITY;
    double ratio;
    
    for ( int i = 0; i < nCol; ++i ) {
        ratio = dx[i] / x[i];
        alphaTmp = HDSDP_MIN(ratio, alphaTmp);
    }
    
    for ( int i = 0; i < nCol; ++i ) {
        ratio = ds[i] / s[i];
        alphaTmp = HDSDP_MIN(ratio, alphaTmp);
    }
    
    ratio = dkappa / kappa;
    alphaTmp = HDSDP_MIN(ratio, alphaTmp);
    
    ratio = dtau / tau;
    alphaTmp = HDSDP_MIN(ratio, alphaTmp);
    
    return fabs(1.0 / alphaTmp);
}

static void LpNewtonUpdate( hdsdp_lpsolver *newton, double *rowDual, double *colVal, double *colDual,
                            double *kappa, double *tau, double spxSize ) {
    
    for ( int i = 0; i < newton->nCol; ++i ) {
        colVal[i] += newton->alpha * newton->dx[i];
        colDual[i] += newton->alpha * newton->ds[i];
    }
    
    for ( int i = 0; i < newton->nRow; ++i ) {
        rowDual[i] += newton->alpha * newton->dy[i];
    }
    
    
    if ( spxSize > 0.0 ) {
        *kappa += newton->alpha * newton->dkappa;
        *tau += newton->alpha * newton->dtau;
        
        double simp = (*kappa) + (*tau);
        
        for ( int i = 0; i < newton->nCol; ++i ) {
            simp += colVal[i];
            simp += colDual[i];
        }
        
        simp = spxSize / simp;
        
        scal(&newton->nRow, &simp, rowDual, &HIntConstantOne);
        scal(&newton->nCol, &simp, colVal, &HIntConstantOne);
        scal(&newton->nCol, &simp, colDual, &HIntConstantOne);
        
        *kappa *= simp;
        *tau *= simp;
    }

    return;
}

/* KKT solver */
static hdsdp_retcode LpNewtonIKKTInit( hdsdp_lpsolver *newton, int nCol, int nRow, int *colMatBeg,
                                 int *colMatIdx, double *colMatElem ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    int nzA = colMatBeg[nCol];
    int ntCol = nRow + nCol;
    int ntNnz = nzA + nCol + nRow;
    
    /*
       [ D^2  A']  +  [ Rp  0  ]
       [ A    0 ]  +  [ 0   Rd ]
     */
    
    /* Elements of D, A and potential regularizers */
    HDSDP_INIT(newton->AugBeg, int, ntCol + 1);
    HDSDP_INIT(newton->AugIdx, int, ntNnz);
    HDSDP_INIT(newton->AugElem, double, ntNnz);
    
    if ( !newton->AugBeg || !newton->AugIdx || !newton->AugElem ) {
        retcode = HDSDP_RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    int *Ap = newton->AugBeg;
    int *Ai = newton->AugIdx;
    double  *Ax = newton->AugElem;
    
    newton->colMatBeg = Ap;
    newton->colMatIdx = Ai;
    newton->colMatElem = Ax;
    
    /* Build up symbolic matrix */
    for ( int i = 0, j; i < nCol; ++i ) {
        Ap[i] = colMatBeg[i] + i;
        Ai[Ap[i]] = i;
        for ( j = colMatBeg[i]; j < colMatBeg[i + 1]; ++j ) {
            Ai[i + j + 1] = colMatIdx[j] + nCol;
            Ax[i + j + 1] = colMatElem[j];
        }
    }
    
    Ap[nCol] = nzA + nCol;
    for ( int i = 0; i < nRow; ++i ) {
        Ap[i + nCol + 1] = Ap[i + nCol] + 1;
        Ai[Ap[nCol] + i] = nCol + i;
    }
    
    /* Create Pardiso solver and run symbolic factorization */
    HDSDP_CALL(HFpLinsysCreate(&newton->kkt, nCol, HDSDP_LINSYS_DENSE_INDEFINITE));
    HFpLinsysSetParam(newton->kkt, -1, -1, newton->nThreads, -1, -1);
    HDSDP_CALL(HFpLinsysSymbolic(newton->kkt, Ap, Ai));
    
exit_cleanup:
    return retcode;
}

static void LpNewtonIKKTLoad( hdsdp_lpsolver *newton, double *colVal, double *colDual ) {
    
    /* Load Newton system */
    if ( colDual ) {
        /* Primal-dual KKT system */
        for ( int iCol = 0; iCol < newton->nCol; ++iCol ) {
            newton->AugElem[newton->AugBeg[iCol]] = HDSDP_DIV(colDual[iCol], colVal[iCol]);
        }
    } else {
        for ( int iCol = 0; iCol < newton->nCol; ++iCol ) {
            newton->AugElem[newton->AugBeg[iCol]] = 1.0;
            for ( int iElem = newton->AugBeg[iCol] + 1; iElem < newton->AugBeg[iCol + 1]; ++iElem ) {
                newton->AugElem[iElem] = newton->colMatElem[iElem - iCol] * colVal[iCol];
            }
        }
    }
    
    return;
}

static void LpNewtonIKKTRegularize( hdsdp_lpsolver *newton, double pReg, double dReg ) {
    
    /* Regularize the augmented system with
     
            [ Rp   0 ]
            [  0  Rd ]
     */
    
    /* Regularize primal. First entries of initial nCol columns */
    if ( pReg > 0.0 ) {
        for ( int i = 0; i < newton->nCol; ++i ) {
            newton->AugElem[newton->AugBeg[i]] += pReg;
        }
    }
    
    /* Regularize dual. Last n entries of the CSC representation of lower-triangular part */
    if ( dReg > 0.0 ) {
        int nNtCol = newton->nCol + newton->nRow;
        for ( int i = 0; i < newton->nRow; ++i ) {
            newton->AugElem[newton->AugBeg[nNtCol - i - 1]] = dReg;
        }
    }
    
    return;
}

static hdsdp_retcode LpNewtonIKKTFactorize( hdsdp_lpsolver *newton ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    HDSDP_CALL(HFpLinsysNumeric(newton->kkt, newton->AugBeg, newton->AugIdx, newton->AugElem));
    
exit_cleanup:
    return retcode;
}

static hdsdp_retcode LpNewtonIKKTPrimalSolve( hdsdp_lpsolver *newton, double *lpObj, double *lpRHS, double *colVal, double *colDual, 
                                              double *resiPrimal, double *resiDual, double *resiMu1 ) {
    
    /* KKT solver for the primal augmented system
     
      [ I   AX' ] [ dx'] = [ -v .* rd / mu + x .* s / mu ]
      [ AX   0  ] [ dy'] = [               rp            ]
     
     After the system is solved, we recover the primal and dual directions by
     
        dx = - v .* dx
        dy = mu * dy'
        ds = -rd - A' * dy
     
    */
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    
    
exit_cleanup:
    return retcode;
}

static hdsdp_retcode LpNewtonIKKTSolve( hdsdp_lpsolver *newton, double *lpObj, double *lpRHS, double *colVal, double *colDual,
                                  double kappa, double tau, double *resiPrimal, double *resiDual, double resiComp,
                                  double *resiMu1, double resiMu2, int isCorr ) {
    
    /* KKT solver for the augmented system
     
     [     A          -b ] [ dy ] = - [  rp ] -> m
     [ -A'      -I     c ] [ dx ] = - [  rd ] -> n
     [  b'   -c'   -1    ] [ ds ] = - [  rk ] -> 1
     [        S     X    ] [ dk ] = - [ rm1 ] -> n
     [              t  k ] [ dt ] = - [ rm2 ] -> 1
     
     The KKT solver reduces the system into
     
     [ X^{-1} S  A' ]
     [   A       0  ]
     
     and performs backward substitution to recover the directions.
     The directions are stored into the solver [dx, dy, ds, dkappa, dtau]
     
     If isCorr is true, then d1 = M \ [-c, b]' is reused in computing the corrector step
     and the directions are stored in [ dxcorr, dycorr, dscorr, dkappacorr, dtaucorr ]
     
    */
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    double *dx = NULL;
    double *dy = NULL;
    double *ds = NULL;
    double *dkappa = NULL;
    double *dtau = NULL;
    
    if ( isCorr ) {
        dx = newton->dxcorr;
        dy = newton->dycorr;
        ds = newton->dscorr;
        dkappa = &newton->dkappacorr;
        dtau = &newton->dtaucorr;
    } else {
        dx = newton->dx;
        dy = newton->dy;
        ds = newton->ds;
        dkappa = &newton->dkappa;
        dtau = &newton->dtau;
    }
    
    /* Set up d1 */
    if ( !isCorr ) {
        for ( int i = 0; i < newton->nCol; ++i ) {
            newton->d1[i] = - lpObj[i];
        }
        HDSDP_MEMCPY(newton->d1 + newton->nCol, lpRHS, double, newton->nRow);
    }
    
    /* Set up d2 */
    HDSDP_ZERO(newton->d2, double, newton->nCol + newton->nRow);
    
    if ( resiDual ) {
        HDSDP_MEMCPY(newton->d2, resiDual, double, newton->nCol);
    }
    
    if ( resiMu1 ) {
        for ( int i = 0; i < newton->nCol; ++i ) {
            newton->d2[i] += HDSDP_DIV(resiMu1[i], colVal[i]);
        }
    }
    
    if ( resiPrimal ) {
        HDSDP_MEMCPY(newton->d2 + newton->nCol, resiPrimal, double, newton->nRow);
    }
    
    /* Solve for d1 and d2 */
    if ( !isCorr ) {
        HDSDP_CALL(HFpLinsysSolve(newton->kkt, 1, newton->d1, NULL));
    }
    HDSDP_CALL(HFpLinsysSolve(newton->kkt, 1, newton->d2, NULL));
    
    /* Retreive the directions */
    double cbTd1 = 0.0;
    double cbTd2 = 0.0;
    
    for ( int i = 0; i < newton->nCol; ++i ) {
        cbTd1 += lpObj[i] * newton->d1[i];
        cbTd2 += lpObj[i] * newton->d2[i];
    }
    
    for ( int i = 0; i < newton->nRow; ++i ) {
        cbTd1 += lpRHS[i] * newton->d1[newton->nCol + i];
        cbTd2 += lpRHS[i] * newton->d2[newton->nCol + i];
    }
    
    double dTauDenom = kappa - cbTd1 * tau;
    double dTauNumer = - (cbTd2 * tau + resiComp * tau + resiMu2);
    
    if ( fabs(dTauDenom) > 1e-10 ) {
        *dtau = dTauNumer / dTauDenom;
    } else {
        if ( dTauDenom < 0 ) {
            *dtau = -HDSDP_DIV(dTauNumer, -dTauDenom);
        } else {
            *dtau = HDSDP_DIV(dTauNumer, dTauDenom);
        }
    }
    
    /* Recover dx and dy */
    for ( int i = 0; i < newton->nCol; ++i ) {
        dx[i] = newton->d1[i] * (*dtau) - newton->d2[i];
        ds[i] = - dx[i] * HDSDP_DIV(colDual[i], colVal[i]);
    }
    
    if ( resiMu1 ) {
        for ( int i = 0; i < newton->nCol; ++i ) {
            ds[i] -= HDSDP_DIV(resiMu1[i], colVal[i]);
        }
    }
    
    for ( int i = 0; i < newton->nRow; ++i ) {
        dy[i] = newton->d2[i + newton->nCol] - newton->d1[i + newton->nCol] * (*dtau);
    }
    
    *dkappa = - kappa * (*dtau) - resiMu2;
    *dkappa = HDSDP_DIV(*dkappa, tau);
    
    /* Overflow from pardiso */
    if ( *dkappa > HDSDP_INFINITY || *dkappa < -HDSDP_INFINITY || (*dkappa != *dkappa) ) {
        retcode = HDSDP_RETCODE_FAILED;
    }
    
exit_cleanup:
    
    return retcode;
}

static void LpNewtonIRegularize( int nCol, int nRow, int nNtCol, int *augMatBeg,
                                int *augMatIdx, double *augMatElem, double pReg, double dReg ) {
    
    /* Regularize the augmented system with
     
            [ Rp   0 ]
            [  0  Rd ]
     */
    
    /* Regularize primal. First entries of initial nCol columns */
    if ( pReg > 0 ) {
        for ( int i = 0; i < nCol; ++i ) {
            augMatElem[augMatBeg[i]] += pReg;
        }
    }
    
    /* Regularize dual. Last n entries of the CSC representation of lower-triangular part */
    if ( dReg > 0 ) {
        for ( int i = 0; i < nRow; ++i ) {
            augMatElem[augMatBeg[nNtCol] - i - 1] = dReg;
        }
    }
    
    return;
}

extern hdsdp_retcode LpNewtonCreate( hdsdp_lpsolver **pnewton, int nThreads ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    if ( !pnewton ) {
        retcode = HDSDP_RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    hdsdp_lpsolver *nt = NULL;
    HDSDP_INIT(nt, hdsdp_lpsolver, 1);
    
    if ( !nt ) {
        retcode = HDSDP_RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    HDSDP_ZERO(nt, hdsdp_lpsolver, 1);
    
    nt->beta = 0.99;
    nt->gamma = 0.7;
    nt->badNewton = 0;
    nt->pReg = PREGULARIZER;
    nt->dReg = DREGULARIZER;
    nt->nThreads = nThreads;
    
    *pnewton = nt;
    
exit_cleanup:
    return retcode;
}

extern hdsdp_retcode LpNewtonInit( hdsdp_lpsolver *newton, int nCol, int nRow, int *colMatBeg, int *colMatIdx, double *colMatElem ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    newton->nCol = nCol;
    newton->nRow = nRow;
    
    /* Initialize KKT solver. We place KKT solver within the LP Newton structure */
    HDSDP_CALL(LpNewtonIKKTInit(newton, nCol, nRow, colMatBeg, colMatIdx, colMatElem));
    
    /* Initialize vectors */
    HDSDP_INIT(newton->dd, double, nCol);
    HDSDP_INIT(newton->d1, double, nCol + nRow);
    HDSDP_INIT(newton->d2, double, nCol + nRow);
    HDSDP_INIT(newton->daux, double, nCol + nRow);
    
    HDSDP_INIT(newton->dx, double, 2 * nCol + nRow);
    HDSDP_INIT(newton->dxcorr, double, 2 * nCol + nRow);
    
    newton->dy = newton->dx + nCol;
    newton->ds = newton->dy + nRow;
    
    newton->dycorr = newton->dxcorr + nCol;
    newton->dscorr = newton->dycorr + nRow;
    
    newton->iterType = MEHROTRA_HSD;
    
    if ( !newton->dd || !newton->d1 || !newton->d2 || !newton->daux || !newton->dx ||
         !newton->dycorr || !newton->dscorr ) {
        retcode = HDSDP_RETCODE_FAILED;
        goto exit_cleanup;
    }
    
exit_cleanup:
    return retcode;
}

/** @brief Implement one Newton's step
 * Currently no centering step is taken. On exit, colVal, rowDual, colDual, kappa, tau are modified
 */
extern hdsdp_retcode LpNewtonOneStep( hdsdp_lpsolver *newton, double *lpObj, double *lpRHS, int *colMatBeg, int *colMatIdx, double *colMatElem,
                                double *colVal, double *rowDual, double *colDual, double *kappa, double *tau,
                                double *pRes, double *dRes, double pObjVal, double dObjVal, double spxSize ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    if ( !newton ) {
        goto exit_cleanup;
    }
    
    int nCol = newton->nCol;
    int nRow = newton->nRow;
    
    /* Prepare iteration */
    double kval = *kappa, tval = *tau;
    double *x = colVal, *s = colDual;
    double alpha = 0.0, beta = newton->beta, gamma = newton->gamma;
    double *daux = newton->daux;
    double *dx = newton->dx, *dy = newton->dy, *ds = newton->ds;
    double dkappa = 0.0, dtau = 0.0;
    double mu = LpNewtonIBarrier(nCol, x, s, kval, tval);
    newton->mu = mu;
    
    /* Add regularization */
    double pReg = HDSDP_MIN(mu * 1e-04, PREGULARIZER);
    
    if ( pReg < 1e-15 ) {
        pReg = 0.0;
    }
    
    double dReg = pReg;
    newton->pReg = pReg;
    newton->dReg = dReg;
    dReg = HDSDP_MIN(dReg, DREGULARIZER);
    
    if ( newton->iterType == PRIMAL ) {
        LpNewtonIKKTLoad(newton, colVal, NULL);
        pReg = 0.0;
        dReg = 0.0;
    } else {
        LpNewtonIKKTLoad(newton, colVal, colDual);
    }
    
    LpNewtonIKKTRegularize(newton, pReg, 0.0);
    
    /* Factorize the augmented system */
    HDSDP_CALL(LpNewtonIKKTFactorize(newton));
    
    if ( newton->iterType == PDIPM_HSD ) {
        
        double xsi = 0.0;
        double xinvsi = 0.0;
        double minXSi = kval * tval;
        double minXinvSi = HDSDP_INFINITY;
        
        for ( int i = 0; i < nCol; ++i ) {
            xsi = x[i] * s[i];
            xinvsi = x[i] / s[i];
            minXinvSi = HDSDP_MIN(xinvsi, minXinvSi);
            minXSi = HDSDP_MIN(minXSi, xsi);
        }
        
        double ksi = minXSi / mu;
        newton->gamma = (1 - ksi) / ksi;
        newton->gamma = 0.05 * newton->gamma;
        newton->gamma = 0.1 * newton->gamma * newton->gamma * newton->gamma;
        newton->gamma = HDSDP_MIN(0.8, newton->gamma);
        newton->gamma = HDSDP_MAX(0.1, newton->gamma);
        
        double mugamma = mu * gamma;
        double resiMu2 = kval * tval - mugamma;
        
        for ( int i = 0; i < nCol; ++i ) {
            daux[i] = x[i] * s[i] - mugamma;
        }
        
        HDSDP_CALL(LpNewtonIKKTSolve(newton, lpObj, lpRHS, colVal, colDual, kval, tval, pRes,
                                   dRes, dObjVal - pObjVal - kval, daux, resiMu2, 0));
        
        
    } else if ( newton->iterType == MEHROTRA_HSD ) {
        
        /* Predictor step */
        double resiMu2 = kval * tval;
        for ( int i = 0; i < nCol; ++i ) {
            daux[i] = x[i] * s[i];
        }
        
        HDSDP_CALL(LpNewtonIKKTSolve(newton, lpObj, lpRHS, colVal, colDual, kval, tval, pRes,
                                   dRes, dObjVal - pObjVal - kval, daux, resiMu2, 0));
        alpha = LpNewtonIRatioTest(nCol, x, newton->dx, s, newton->ds, kval,
                                   newton->dkappa, tval, newton->dtau);
        
        double muAff = ( kval + alpha * newton->dkappa ) * ( tval + alpha * newton->dtau );
        
        /* Compute sigma */
        for ( int i = 0; i < nCol; ++i ) {
            muAff += ( x[i] + alpha * newton->dx[i] ) * ( s[i] + alpha * newton->ds[i] );
        }
        
        muAff = muAff / ( nCol + 1 );
        gamma = HDSDP_DIV(muAff, mu);
        gamma = gamma * gamma * gamma;
        
        if ( mu < 1e-10 ) {
            gamma = HDSDP_MAX(gamma, 0.3);
        }
        
#ifdef IPM_DEBUG
        printf("Mehrotra gamma: %f ", gamma);
#endif
        
        double mugamma = mu * gamma;
        
        /* Corrector step */
        resiMu2 = newton->dkappa * newton->dtau - mugamma;
        for ( int i = 0; i < nCol; ++i ) {
            daux[i] = newton->dx[i] * newton->ds[i] - mugamma;
        }
        
        HDSDP_CALL(LpNewtonIKKTSolve(newton, lpObj, lpRHS, colVal, colDual, kval, tval,
                                   NULL, NULL, 0.0, daux, resiMu2, 1));
        
        /* Combine predictor and corrector */
        newton->dkappa += newton->dkappacorr;
        newton->dtau += newton->dtaucorr;
        for ( int i = 0; i < nCol; ++i ) {
            dx[i] += newton->dxcorr[i];
            ds[i] += newton->dscorr[i];
        }
        
        for ( int i = 0; i < nRow; ++i ) {
            dy[i] += newton->dycorr[i];
        }
        
    } else if ( newton->iterType == NOCORR_HSD ) {
                
        newton->gamma = 0.8;
        newton->beta = 0.1;
        
        double mugamma = mu * gamma;
        double resiMu2 = kval * tval - mugamma;
        
        for ( int i = 0; i < nCol; ++i ) {
            daux[i] = x[i] * s[i] - mugamma;
        }
        
        HDSDP_CALL(LpNewtonIKKTSolve(newton, lpObj, lpRHS, colVal, colDual, kval, tval, pRes,
                                   dRes, dObjVal - pObjVal - kval, daux, resiMu2, 0));
    }
    
    dkappa = newton->dkappa;
    dtau = newton->dtau;
    alpha = LpNewtonIRatioTest(nCol, x, dx, s, ds, kval, dkappa, tval, dtau);
    
    if ( alpha < 1e-04 ) {
        if ( newton->badNewton >= 2 ) {
            retcode = HDSDP_RETCODE_FAILED;
            goto exit_cleanup;
        } else {
            /* Try to rescue the system using scaling and matching */
            HFpLinsysSwitchToBackUp(newton->kkt);
            newton->badNewton += 1;
        }
    }
    
    alpha = alpha * beta;
    alpha = HDSDP_MIN(alpha, 1.0);
    newton->alpha = alpha;
    LpNewtonUpdate(newton, rowDual, colVal, colDual, kappa, tau, spxSize);
    
#ifdef IPM_DEBUG
    printf("Mu: %10.3e Barrier stepsize: %f \n", mu, alpha);
#endif
    
exit_cleanup:
    return retcode;
}

extern void LpNewtonClear( hdsdp_lpsolver *newton ) {
    
    if ( !newton ) {
        return;
    }
    
    HFpLinsysDestroy(&newton->kkt);
    
    HDSDP_FREE(newton->AugBeg);
    HDSDP_FREE(newton->AugIdx);
    HDSDP_FREE(newton->AugElem);
    HDSDP_FREE(newton->colBackup);
    
    HDSDP_FREE(newton->dd);
    HDSDP_FREE(newton->xse);
    HDSDP_FREE(newton->d1);
    HDSDP_FREE(newton->d2);
    HDSDP_FREE(newton->daux);
    
    /* dy and ds are following dx and we do need to free them */
    HDSDP_FREE(newton->dx);
    HDSDP_FREE(newton->dxcorr);
    
    HDSDP_ZERO(newton, hdsdp_lpsolver, 1);
    return;
}

extern void LpNewtonDestroy( hdsdp_lpsolver **pnewton ) {
    
    if ( !pnewton ) {
        return;
    }
    
    LpNewtonClear(*pnewton);
    HDSDP_FREE(*pnewton);
    
    return;
}
