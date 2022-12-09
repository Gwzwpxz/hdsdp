#ifndef lp_filter_h
#define lp_filter_h

#include "pot_def.h"

typedef struct {
    
    int nCone;
    int nActCol;
    
    int *coneOrder;
    double *coneFilter;
    
    /* Data */
    double *coneVals;
    
} col_filter;

#ifdef __cplusplus
extern "C" {
#endif

extern pot_int ConeFilterCreate( col_filter **pcolFilter );
extern pot_int ConeFilterInit( col_filter *colFilter, pot_int nCone, double *coneData );
extern void ConeFilterClear( col_filter *colFilter );
extern void ConeFilterBuildUp( col_filter *colFilter );
extern void ConeFilterZeroOut( col_filter *colFilter, double *xVal );
extern void ConeFilterDestroy( col_filter **pcolFilter );

#ifdef __cplusplus
}
#endif

#endif /* lp_filter_h */
