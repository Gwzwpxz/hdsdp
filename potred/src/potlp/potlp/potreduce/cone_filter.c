#include "pot_utils.h"
#include "cone_filter.h"

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

#ifdef FILTER_DEBUG
#undef FILTER_DEBUG
#define FILTER_DEBUG printf
#else
#define FILTER_DEBUG(...)
#endif
#define FILTER_RATIO  (0.1)
extern void ConeFilterBuildUp( col_filter *colFilter, double threshVal ) {
    /* Build up filter for the current iterate */
    
    if ( !colFilter ) {
        return;
    }
    
    double *xVal = colFilter->coneVals;
    double *coneFilter = colFilter->coneFilter;
    
    int nCone = colFilter->nCone;
    int nFilters = nCone * FILTER_RATIO;
    int nActive = 0;
    int *coneOrder = colFilter->coneOrder;
    
    for ( int i = 0; i < nCone; ++i ) {
        coneOrder[i] = i;
    }
    
    POTLP_MEMCPY(coneFilter, xVal, double, nCone);
    potUtilSortbyDbl(coneOrder, coneFilter, 0, nCone - 1);
    for ( nActive = 0; nActive <= nFilters; ++nActive ) {
        if ( coneFilter[nActive] > threshVal ) {
            break;
        }
    }
    
    colFilter->nActCol = nActive;
    FILTER_DEBUG("%d variables not in support. \n", nActive);
    
    return;
}

extern void ConeFilterZeroOut( col_filter *colFilter, double *vVal ) {
    
    if ( !colFilter ) {
        return;
    }
    
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
