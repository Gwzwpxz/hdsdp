#include "pot_utils.h"
#include "lp_filter.h"

extern pot_int ConeFilterCreate( col_filter **pcolFilter ) {
    
    pot_int retcode = RETCODE_OK;
    
    if ( !pcolFilter ) {
        retcode = RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    col_filter *colFilter = NULL;
    POTLP_INIT(colFilter, col_filter, 1);
    
    if ( !colFilter ) {
        retcode = RETCODE_OK;
        goto exit_cleanup;
    }
    
    POTLP_ZERO(colFilter, col_filter, 1);
    *pcolFilter = colFilter;
    
exit_cleanup:
    
    return retcode;
}

extern pot_int ConeFilterInit( col_filter *colFilter, pot_int nCone, double *coneData ) {
    
    pot_int retcode = RETCODE_OK;
    
    if ( !coneData ) {
        retcode = RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    colFilter->nCone = nCone;
    
    colFilter->coneVals = coneData;
    POTLP_INIT(colFilter->coneFilter, double, nCone);
    POTLP_INIT(colFilter->coneOrder, int, nCone);
    
exit_cleanup:
    
    return retcode;
}

extern void ConeFilterClear( col_filter *colFilter ) {
    
    if ( !colFilter ) {
        return;
    }
    
    POTLP_FREE(colFilter->coneFilter);
    POTLP_FREE(colFilter->coneOrder);
    POTLP_ZERO(colFilter, col_filter, 1);
    
    return;
}

#define FILTER_RATIO  (0.300)
#define FILTER_THRESH (1e-03)
extern void ConeFilterBuildUp( col_filter *colFilter ) {
    /* Build up filter for the current iterate */
    double *xVal = colFilter->coneVals;
    
    int nCone = colFilter->nCone;
    int nFilters = nCone * FILTER_RATIO;
    int nActive = 0;
    int *coneOrder = colFilter->coneOrder;
    
    for ( int i = 0; i < nCone; ++i ) {
        coneOrder[i] = i;
    }
    
    potUtilSortbyDbl(coneOrder, xVal, 0, nCone);
    for ( nActive = 0; ( xVal[nActive] < FILTER_THRESH ) &&
                       ( nActive <= nFilters ); ++nActive );
    
    return;
}

extern void ConeFilterZeroOut( col_filter *colFilter, double *vVal ) {
    
    for ( int i = 0; i < colFilter->nActCol; ++i ) {
        vVal[colFilter->coneOrder[i]] = 0.0;
    }
    
    return;
}

extern void ConeFilterDestroy( col_filter **pcolFilter ) {
    
    if ( !pcolFilter ) {
        return;
    }
    
    ConeFilterClear(*pcolFilter);
    POTLP_FREE(*pcolFilter);
    
    return;
}
