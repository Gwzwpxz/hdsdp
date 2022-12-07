/** @file pot\_utils.c
 *  @brief Implement the potential reduction utilities
 *
 * @TODO: Add more detailed comments
 */

#ifndef pot_utils_h
#define pot_utils_h

#include "pot_def.h"
#include "pot_structs.h"

#define error_traceback(info)                                               \
if ( retcode != RETCODE_OK ) {                                              \
    printf("Error:File:%s -> Line:%d -> %s \n", __FILE__, __LINE__, info); \
}

extern double potUtilGetTimeStamp( void );
extern void potUtilGetDefaultParams( double dblParams[NUM_DBL_PARAM], int intParams[NUM_INT_PARAM] );
extern void potUtilPrintParams( double dblParams[NUM_DBL_PARAM], int intParams[NUM_INT_PARAM] );
extern void potUtilPrintIParams( double dblParams[NUM_DBL_PARAM], int intParams[NUM_INT_PARAM] );

/* Debugging */
extern void potUtilPrintDblContent( int n, double *d );
extern void potUtilDumpDblMatrix( int m, int n, double *d );
extern void potUtilDumpSymMatrix( int m, int n, double *d );
extern void potUtilPrintIntContent( int n, int *d );
extern void potUtilPrintDblMin( int n, double *d );
extern void potUtilPrintDblSum( int n, double *d );
extern int potUtilVerifyNeighbour( int n, double *d, double r );

#endif /* pot_utils_h */
