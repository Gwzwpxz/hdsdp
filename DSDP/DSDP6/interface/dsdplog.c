#include <stdio.h>
#include "dsdplog.h"
#include "hsd.h"
#define dsdplog printf

#ifndef isnan
#define isnan(x) ((x) != (x))
#endif

#ifndef EVENT
#define E(x) dsdpSolver->eventMonitor(x);
#endif

static char etype[] = "DSDP logging";

extern inline double my_clock(void) {
    struct timeval t;
    gettimeofday(&t, NULL);
    return (1.0e-6 * t.tv_usec + t.tv_sec);
}

static void dsdpstatus( DSDP_INT code, char *word ) {
    
    // Convert status code
    switch (code) {
        case DSDP_UNKNOWN:
            strcpy(word, "DSDP_UKNOWN");
            break;
        case DSDP_OPTIMAL:
            strcpy(word, "DSDP_OPTIMAL");
            break;
        case DSDP_MAXITER:
            strcpy(word, "DSDP_MAXITER");
            break;
        case DSDP_INTERNAL_ERROR:
            strcpy(word, "DSDP_INTERNAL_ERROR");
            break;
        case DSDP_PD_FEASIBLE:
            strcpy(word, "DSDP_PRIMAL_DUAL_FEASIBLE");
            break;
        case DSDP_PFEAS_DINFEAS:
            strcpy(word, "DSDP_PFEASIBLE_DINFEASIBLE");
            break;
        case DSDP_PUNKNOWN_DFEAS:
            strcpy(word, "DSDP_PUNKNOWN_DFEASIBLE");
            break;
        case DSDP_PUNKNOWN_DINFEAS:
            strcpy(word, "DSDP_PUNKNOWN_DINFEASIBLE");
            break;
        case DSDP_PINFEAS_DFEAS:
            strcpy(word, "DSDP_PINFEASIBLE_DFEASIBLE");
            break;
        default:
            break;
    }
}

extern void dsdpshowdash( void ) {
    printf("---------------------------------------"
           "---------------------------------------"
           "----------------------\n");
}

extern void dsdpCheckNan( HSDSolver *dsdpSolver ) {
    
    // Check nan in some iterations
    if (isnan(dsdpSolver->y->x[0])) {
        dsdpSolver->eventMonitor[EVENT_NAN_IN_ITER] = TRUE;
        return;
    }
    
    if (isnan(dsdpSolver->tau)   || isnan(dsdpSolver->kappa) ||
        isnan(dsdpSolver->alpha) || isnan(dsdpSolver->Ry)    ||
        isnan(dsdpSolver->dObjVal|| isnan(dsdpSolver->pObjVal))) {
        dsdpSolver->eventMonitor[EVENT_NAN_IN_ITER] = TRUE;
    }
}

extern void dsdpprintPhaseAheader( void ) {
    
    printf("| DSDP Phase A starts. Eliminating dual infeasibility %45s|\n", "");
    dsdpshowdash();
    /*      --------------------------------------------------------------------------------*/
    /*      | niter |   pObj   |   dObj   |  dInf  |  k/t  |  mu  |  alpha  |  pNrm  |  E  |*/
    /*      --------------------------------------------------------------------------------*/
    printf("| %4s | %12s | %12s | %9s | %8s | %8s | %8s | %8s | %3s |\n",
            "Iter", "pObj", "dObj", "dInf", "k/t", "mu", "alpha", "pNrm", "E");
    dsdpshowdash();
}

extern void dsdpprintPhaseBheader( void ) {
    
    dsdpshowdash();
    /*      --------------------------------------------------------------------------------*/
    /*      | niter |     pObj     |     dObj     |  mu  |  alpha  |  pNrm  |  dPot  |  E  |*/
    /*      --------------------------------------------------------------------------------*/
    printf("| %4s | %17s | %17s | %10s | %8s | %8s | %8s | %3s |\n",
            "Iter", "pObj", "dObj", "mu", "alpha", "pNrm", "dPot", "E");
    dsdpshowdash();
}

extern void DSDPResetPhaseAMonitor( HSDSolver *dsdpSolver ) {
    
    // Reset event and iteration monitor in Phase A
    DSDP_INT islogged = dsdpSolver->iterProgress[ITER_LOGGING];
    DSDP_INT isdobj = dsdpSolver->iterProgress[ITER_DUAL_OBJ];
    
    memset(dsdpSolver->iterProgress, 0, sizeof(DSDP_INT) * IterStep);
    dsdpSolver->iterProgress[ITER_INITIALIZE] = TRUE;
    dsdpSolver->iterProgress[ITER_LOGGING] = islogged;
    dsdpSolver->iterProgress[ITER_DUAL_OBJ] = isdobj;
    
    DSDP_INT *monitor = dsdpSolver->eventMonitor;
    assert( monitor[EVENT_IN_PHASE_A] );
    monitor[EVENT_PFEAS_FOUND] = FALSE;
    monitor[EVENT_LARGE_NORM] = FALSE;
    monitor[EVENT_SMALL_STEP] = FALSE;
}

extern void DSDPResetPhaseBMonitor( HSDSolver *dsdpSolver ) {
    
    // Reset event and iteration monitor in Phase B
    DSDP_INT islogged = dsdpSolver->iterProgress[ITER_LOGGING];
    DSDP_INT isdobj = dsdpSolver->iterProgress[ITER_DUAL_OBJ];
    
    memset(dsdpSolver->iterProgress, 0, sizeof(DSDP_INT) * IterStep);
    dsdpSolver->iterProgress[ITER_INITIALIZE] = TRUE;
    dsdpSolver->iterProgress[ITER_RESIDUAL] = TRUE;
    dsdpSolver->iterProgress[ITER_LOGGING] = islogged;
    dsdpSolver->iterProgress[ITER_DUAL_OBJ] = isdobj;
    
    DSDP_INT *monitor = dsdpSolver->eventMonitor;
    assert( monitor[EVENT_IN_PHASE_B] );
    
    monitor[EVENT_PFEAS_FOUND] = FALSE;
    monitor[EVENT_LARGE_NORM]  = FALSE;
    monitor[EVENT_SMALL_STEP]  = FALSE;
}

extern DSDP_INT DSDPPhaseALogging( HSDSolver *dsdpSolver ) {
    
    /*  --------------------------------------------------------------------------------*/
    /*  | niter |   pObj   |   dObj   |  dInf  |  k/t  |  mu  |  alpha  |  pNrm  |  E  |*/
    /*  --------------------------------------------------------------------------------*/
    
    // Implement the logging procedure of DSDP phase A
    DSDP_INT retcode = DSDP_RETCODE_OK;
    retcode = checkIterProgress(dsdpSolver, ITER_LOGGING);
    
    char event[3] = " ";
    double nRy = fabs(dsdpSolver->Ry) * sqrt(dsdpSolver->n);
    double tau = dsdpSolver->tau;
    double kovert = dsdpSolver->kappa / tau;
    
    DSDP_INT *moniter = dsdpSolver->eventMonitor;
    
    assert( moniter[EVENT_IN_PHASE_A] );
    
    if (moniter[EVENT_PFEAS_FOUND]) {
        strcpy(&event[0], "P");
    }
    
    if (moniter[EVENT_NO_RY]) {
        strcpy(&event[0], "*");
        goto print_log;
    }
    
    if (moniter[EVENT_MU_QUALIFIES] &&
        moniter[EVENT_KT_QUALIFIES]) {
        strcpy(&event[0], "I");
        goto print_log;
    }
    
    if (dsdpSolver->smallIter == 3 ||
        moniter[EVENT_NAN_IN_ITER] ||
        moniter[EVENT_BAD_SCHUR]) {
        strcpy(&event[0], "F");
        goto print_log;
    }
    
    if (moniter[EVENT_LARGE_NORM]) {
        strcpy(&event[1], "N");
    }
    
    if (moniter[EVENT_MAX_ITERATION]) {
        strcpy(&event[1], "M");
    }
    
print_log:
    printf("| %4d | %12.3e | %12.3e | %9.2e | %8.2e | %8.2e | %8.2e | %8.2e | %3s |\n",
            dsdpSolver->iterA + 1, dsdpSolver->pObjVal, dsdpSolver->dObjVal / tau,
            nRy / tau, kovert, dsdpSolver->mu, dsdpSolver->alpha, dsdpSolver->Pnrm,
            event);
    
    dsdpSolver->iterProgress[ITER_LOGGING] = TRUE;
    return retcode;
}

extern DSDP_INT DSDPPhaseBLogging( HSDSolver *dsdpSolver ) {
    
    /*
     ----------------------------------------------------------------------------------------
    | Iter |         pObj |         dObj |        mu |    alpha |     pNrm |     dPot |   E |
    -----------------------------------------------------------------------------------------
     */
    
    // Implement the logging procedure of DSDP phase A
    DSDP_INT retcode = DSDP_RETCODE_OK;
    retcode = checkIterProgress(dsdpSolver, ITER_LOGGING);
    
    char event[3] = " ";
    DSDP_INT *moniter = dsdpSolver->eventMonitor;
    
    assert( moniter[EVENT_IN_PHASE_B] );
    
    if (moniter[EVENT_PFEAS_FOUND]) {
        strcpy(&event[0], "P");
    }
    
    if (dsdpSolver->solStatus == DSDP_OPTIMAL) {
        strcpy(&event[0], "*");
        goto print_log;
    }
    
    if (dsdpSolver->solStatus == DSDP_PINFEAS_DFEAS) {
        strcpy(&event[0], "I");
        goto print_log;
    }
    
    if (dsdpSolver->smallIter == 3 ||
        moniter[EVENT_NAN_IN_ITER] ||
        moniter[EVENT_BAD_SCHUR]   ||
        dsdpSolver->solStatus == DSDP_INTERNAL_ERROR) {
        strcpy(&event[0], "F");
        goto print_log;
    }
    
    if (moniter[EVENT_LARGE_NORM]) {
        strcpy(&event[1], "N");
    }
    
    if (moniter[EVENT_MAX_ITERATION]) {
        strcpy(&event[1], "M");
    }
    
    if (moniter[EVENT_INVALID_GAP]) {
        strcpy(&event[1], "!");
    }
    
print_log:
    printf("| %4d | %17.10e | %17.10e | %10.3e | %8.2e | %8.2e | %8.1e | %3s |\n",
            dsdpSolver->iterB + 1, dsdpSolver->pObjVal, dsdpSolver->dObjVal,
            dsdpSolver->mu, dsdpSolver->alpha, dsdpSolver->Pnrm, dsdpSolver->dPotential,
            event);
    
    dsdpSolver->iterProgress[ITER_LOGGING] = TRUE;
    return retcode;
}

extern DSDP_INT DSDPCheckPhaseAConvergence( HSDSolver *dsdpSolver, DSDP_INT *isOK ) {
    
    /* Check convergence of DSDP */
    DSDP_INT retcode = DSDP_RETCODE_OK;
    DSDP_INT *monitor = dsdpSolver->eventMonitor;
    *isOK = FALSE;
    
    assert( dsdpSolver->eventMonitor[EVENT_IN_PHASE_A] );
    
    // Dual infeasibility eliminated
    if (fabs(dsdpSolver->Ry) * sqrt(dsdpSolver->n) < \
        dsdpSolver->param->absOptTol * 0.1 * dsdpSolver->tau) {
        monitor[EVENT_NO_RY] = TRUE;
        
        if (dsdpSolver->pObjVal != DSDP_INFINITY) {
            dsdpSolver->solStatus = DSDP_PD_FEASIBLE;
        } else {
            dsdpSolver->solStatus = DSDP_PUNKNOWN_DFEAS;
        }
        
        *isOK = TRUE;
        return retcode;
    }
    
    // Dual infeasibility certificate found
    if (dsdpSolver->tau < 0.001 * dsdpSolver->kappa) {
        monitor[EVENT_KT_QUALIFIES] = TRUE;
    }
    
    if (dsdpSolver->mu < dsdpSolver->param->absOptTol) {
        monitor[EVENT_MU_QUALIFIES] = TRUE;
    }
    
    if (monitor[EVENT_KT_QUALIFIES] &&
        monitor[EVENT_MU_QUALIFIES]) {
        if (dsdpSolver->pObjVal != DSDP_INFINITY) {
            dsdpSolver->solStatus = DSDP_PFEAS_DINFEAS;
        } else {
            dsdpSolver->solStatus = DSDP_PUNKNOWN_DINFEAS;
        }
        *isOK = TRUE;
        return retcode;
    }

    // Maximum iteration
    if (dsdpSolver->iterA >= dsdpSolver->param->AmaxIter) {
        monitor[EVENT_MAX_ITERATION] = TRUE;
        dsdpSolver->solStatus = DSDP_MAXITER;
        *isOK = TRUE;
        return retcode;
    }
    
    // Small step
    if (dsdpSolver->iterA > 0 && dsdpSolver->alpha < 1e-04 &&
        !monitor[EVENT_PFEAS_FOUND]) {
        monitor[EVENT_SMALL_STEP] = TRUE;
        dsdpSolver->smallIter += 1;
    }
    
    // NAN in iteration
    if (monitor[EVENT_NAN_IN_ITER]) {
        dsdpSolver->solStatus = DSDP_INTERNAL_ERROR;
        *isOK = TRUE;
        return retcode;
    }
    
    if (monitor[EVENT_BAD_SCHUR]) {
        dsdpSolver->solStatus = DSDP_INTERNAL_ERROR;
        *isOK = TRUE;
        return retcode;
    }

    return retcode;
}

extern DSDP_INT printPhaseASummary( HSDSolver *dsdpSolver, double time ) {
    
    // Summarize Phase A iteration
    DSDP_INT retcode = DSDP_RETCODE_OK;
    char sAlog[80] = "| DSDP Phase A ends with status: ";
    dsdpshowdash();
    dsdpstatus(dsdpSolver->solStatus, &sAlog[33]);
    printf("%-80s %18s|\n", sAlog, "");
    printf("| Elapsed Time: %6.5f seconds %65s |\n", time, "");
    dsdpshowdash();
    return retcode;
}

extern DSDP_INT printPhaseBSummary( HSDSolver *dsdpSolver, double time ) {
    
    // Summarize Phase B iteration
    DSDP_INT retcode = DSDP_RETCODE_OK;
    char sBlog[80] = "| DSDP Phase B ends with status: ";
    dsdpshowdash();
    dsdpstatus(dsdpSolver->solStatus, &sBlog[33]);
    printf("%-80s %18s|\n", sBlog, "");
    printf("| Elapsed Time: %6.5f seconds %65s |\n", time, "");
    dsdpshowdash();
    return retcode;
}

extern DSDP_INT DSDPCheckPhaseBConvergence( HSDSolver *dsdpSolver, DSDP_INT *isOK ) {
    
    /* Check convergence of DSDP */
    DSDP_INT retcode = DSDP_RETCODE_OK;
    DSDP_INT *monitor = dsdpSolver->eventMonitor;
    
    assert( dsdpSolver->eventMonitor[EVENT_IN_PHASE_B] );
    *isOK = FALSE;
    
    // Dual objective exceeds primal
    /*
       This may happen due to the perturbation on Schur matrix but once instability
       is detected the perturbation is swiched off.
     */
    
    // Primal infeasibility detected
    if (monitor[EVENT_PINFEAS_DETECTED]) {
        dsdpSolver->solStatus = DSDP_PINFEAS_DFEAS;
        *isOK = TRUE;
        return retcode;
    }
    
    if (dsdpSolver->dObjVal > dsdpSolver->pObjVal) {
        if (!dsdpSolver->gapBroken) {
            dsdpSolver->gapBroken = TRUE;
            dsdpSolver->pObjVal = dsdpSolver->dObjVal + 1e+08;
            monitor[EVENT_INVALID_GAP] = TRUE;
        } else {
            dsdpSolver->solStatus = DSDP_INTERNAL_ERROR;
            *isOK = TRUE;
            return retcode;
        }
    }
    
    // Gap is sufficiently small
    double gap = fabs(dsdpSolver->pObjVal - dsdpSolver->dObjVal);
    gap = gap / (fabs(dsdpSolver->pObjVal) + fabs(dsdpSolver->dObjVal) + 1);
    
    if (dsdpSolver->dObjVal > 1e+08) {
        dsdpSolver->eventMonitor[EVENT_LARGE_DOBJ] = TRUE;
    }
    
    if (gap < 1e-06) {
        monitor[EVENT_MU_QUALIFIES] = TRUE;
        dsdpSolver->solStatus = DSDP_OPTIMAL;
        *isOK = TRUE;
        return retcode;
    }
    
    // Maximum iteration
    if (dsdpSolver->iterB >= dsdpSolver->param->BmaxIter) {
        monitor[EVENT_MAX_ITERATION] = TRUE;
        dsdpSolver->solStatus = DSDP_MAXITER;
        *isOK = TRUE;
    }
    
    // Small step
    if (dsdpSolver->iterB > 0 && dsdpSolver->alpha < 1e-04 &&
        !monitor[EVENT_PFEAS_FOUND]) {
        monitor[EVENT_SMALL_STEP] = TRUE;
        dsdpSolver->smallIter += 1;
    }
    
    if (dsdpSolver->smallIter == 3) {
        dsdpSolver->solStatus = DSDP_INTERNAL_ERROR;
        return retcode;
    }
    
    // NAN in iteration
    if (monitor[EVENT_NAN_IN_ITER]) {
        dsdpSolver->solStatus = DSDP_INTERNAL_ERROR;
        *isOK = TRUE;
        return retcode;
    }
    
    if (monitor[EVENT_BAD_SCHUR]) {
        dsdpSolver->solStatus = DSDP_INTERNAL_ERROR;
        *isOK = TRUE;
        return retcode;
    }

    return retcode;
}

extern DSDP_INT printPhaseABConvert( HSDSolver *dsdpSolver, DSDP_INT *goPb ) {
    
    // Print the conversion logging between phase A and B
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    assert( dsdpSolver->eventMonitor[EVENT_IN_PHASE_A] );
    DSDP_INT status = dsdpSolver->solStatus;
    
    *goPb = TRUE;
    
    switch (status) {
        case DSDP_PD_FEASIBLE:
            printf("| DSDP Phase A certificates primal-dual feasibility %47s|\n", "");
            break;
        case DSDP_PUNKNOWN_DFEAS:
            printf("| DSDP Phase A certificates dual feasibility %54s|\n", "");
            break;
        case DSDP_PFEAS_DINFEAS:
            printf("| DSDP Phase A certificates dual infeasibility "
                   "and primal feasibility \n");
            *goPb = FALSE;
            break;
        case DSDP_PUNKNOWN_DINFEAS:
            printf("| DSDP Phase A certificates dual infeasibility \n");
            *goPb = FALSE;
            break;
        case DSDP_INTERNAL_ERROR:
            printf("| DSDP Phase A encounters internal error \n");
            *goPb = FALSE;
            break;
        case DSDP_MAXITER:
            printf("| DSDP Phase A reaches maximum iteration \n");
            *goPb = FALSE;
            break;
        default:
            error(etype, "Invalid status code at the end of Phase A. \n");
            break;
    }
    
    return retcode;
}
