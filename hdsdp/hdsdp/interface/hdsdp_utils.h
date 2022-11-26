/** @file hdsdp\_utils.h
 *  @brief Define utilities for HDSDP
 *  @date 11/24/2022
 */

#ifndef hdsdp_utils_h
#define hdsdp_utils_h

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "hdsdp.h"

/* Define macros */
#define HDSDP_FREE(var) do {free((var)); (var) = NULL;} while (0)
#define HDSDP_INIT(var, type, size) (var) = (type *) calloc(size, sizeof(type))
#define HDSDP_MEMCPY(dst, src, type, size) memcpy(dst, src, sizeof(type) * (size))
#define HDSDP_ZERO(var, type, size) memset(var, 0, sizeof(type) * (size))

#define HDSDP_CALL(func) retcode = (func);                      \
                         if (retcode != HDSDP_RETCODE_OK) {     \
                             goto exit_cleanup;                 \
                         }

#define HDSDP_MAX(x, y) (x) >= (y) ? (x) : (y);
#define HDSDP_MIN(x, y) (x) <= (y) ? (x) : (y);

#ifdef __cplusplus
extern "C" {
#endif

extern double HUtilGetTimeStamp( void );
extern void HUtilPrintDblContent( int n, double *d );
extern void HUtilPrintIntContent( int n, int *d );
extern void HUtilPrintDblSum( int n, double *d );

#ifdef __cplusplus
}
#endif

#endif /* hdsdp_utils_h */
