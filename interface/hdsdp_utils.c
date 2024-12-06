/** @file hdsdp\_utils.c
 *  @brief Implement the utilities for HDSDP
 *  @date 11/24/2022
 */
#ifdef HEADERPATH
#include "interface/hdsdp_utils.h"
#include "interface/hdsdp_schur.h"
#else
#include "hdsdp_utils.h"
#include "hdsdp_schur.h"
#endif

#include <math.h>
#include <sys/time.h>
#include <signal.h>

static int isCtrlC = 0;
static void monitorCtrlC( int sigNum ) {
    isCtrlC = 1;
    return;
}

static struct sigaction act;

/* TODO: Add compatibility for Windows platform */
static double my_clock( void ) {
    
    struct timeval t;
    gettimeofday(&t, NULL);
    return (1e-06 * t.tv_usec + t.tv_sec);
}

static int dpartitioni( int *ind, double *val, int l, int h ) {
    
    double tmp2 = 0.0, p = val[l];
    int tmp = l, tmp3 = 0;
    
    while ( l < h ) {
        
        while ( l < h && val[h] >= p ) { --h; }
        while ( l < h && val[l] <= p ) { ++l; }
        
        if ( l < h ) {
            tmp2 = val[l]; val[l] = val[h]; val[h] = tmp2;
            tmp3 = ind[l]; ind[l] = ind[h]; ind[h] = tmp3;
        }
    }
    
    tmp2 = val[l]; val[l] = val[tmp]; val[tmp] = tmp2;
    tmp3 = ind[l]; ind[l] = ind[tmp]; ind[tmp] = tmp3;
    
    return l;
}

static int ipartitiond( double *ind, int *val, int l, int h ) {
    
    int tmp2 = 0, p = val[l], tmp = l;
    double tmp3;
    
    while ( l < h ) {
        
        while ( l < h && val[h] >= p ) { --h; }
        while ( l < h && val[l] <= p ) { ++l; }
        
        if ( l < h ) {
            tmp2 = val[l]; val[l] = val[h]; val[h] = tmp2;
            tmp3 = ind[l]; ind[l] = ind[h]; ind[h] = tmp3;
        }
    }
    
    tmp2 = val[l]; val[l] = val[tmp]; val[tmp] = tmp2;
    tmp3 = ind[l]; ind[l] = ind[tmp]; ind[tmp] = tmp3;
    
    return l;
}

static int ipartitioni( int *ind, int *val, int l, int h ) {
    
    int tmp = l, tmp2 = 0, tmp3, p = val[l];
    
    while ( l < h ) {
        
        while ( l < h && val[h] <= p ) { --h; }
        while ( l < h && val[l] >= p ) { ++l; }
        
        if ( l < h ) {
            tmp2 = val[l]; val[l] = val[h]; val[h] = tmp2;
            tmp3 = ind[l]; ind[l] = ind[h]; ind[h] = tmp3;
        }
    }
    
    tmp2 = val[l]; val[l] = val[tmp]; val[tmp] = tmp2;
    tmp3 = ind[l]; ind[l] = ind[tmp]; ind[tmp] = tmp3;
    
    return l;
}

static int dpartitiond( double *ind, double *val, int l, int h ) {
    
    int tmp = l;
    double tmp2 = 0.0, tmp3, p = val[l];
    
    while ( l < h ) {
        
        while ( l < h && val[h] >= p ) { --h; }
        while ( l < h && val[l] <= p ) { ++l; }
        
        if ( l < h ) {
            tmp2 = val[l]; val[l] = val[h]; val[h] = tmp2;
            tmp3 = ind[l]; ind[l] = ind[h]; ind[h] = tmp3;
        }
    }
    
    tmp2 = val[l]; val[l] = val[tmp]; val[tmp] = tmp2;
    tmp3 = ind[l]; ind[l] = ind[tmp]; ind[tmp] = tmp3;
    
    return l;
}

#define LINEBUFFER  (1024)
static hdsdp_retcode HUtilIGetIntVecData( char *fname, int nElem, int *iData ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    char thisLine[128] = "";
    FILE *file = NULL;
    
    file = fopen(fname, "r");
    
    if (!file) {
        printf("Failed to read %s \n", fname);
        retcode = HDSDP_RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    int iLine = 0;
    int nGet = 0;
    
    for ( iLine = 0; !feof(file); ) {
        
        fgets(thisLine, LINEBUFFER, file);
        
        if ( iLine >= nElem ) {
            break;
        }
        
        nGet = sscanf(thisLine, "%d", &iData[iLine]);
        
        if ( nGet != 1 ) {
            printf("Error at line %d of %s \n", iLine + 1, fname);
            retcode = HDSDP_RETCODE_FAILED;
            goto exit_cleanup;
        }
    
        iLine += 1;
    }
    
exit_cleanup:
    return retcode;
}

static hdsdp_retcode HUtilIGetDblVecData( char *fname, int nElem, double *dData ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    char thisLine[128] = "";
    FILE *file = NULL;
    
    file = fopen(fname, "r");
    
    if (!file) {
        printf("Failed to read %s \n", fname);
        retcode = HDSDP_RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    int iLine = 0;
    int nGet = 0;
    
    for ( iLine = 0; !feof(file); ) {
        
        fgets(thisLine, LINEBUFFER, file);
        
        if ( iLine >= nElem ) {
            break;
        }
        
        nGet = sscanf(thisLine, "%lg", &dData[iLine]);
        
        if ( nGet != 1 ) {
            printf("Error at line %d of %s \n", iLine + 1, fname);
            retcode = HDSDP_RETCODE_FAILED;
            goto exit_cleanup;
        }
    
        iLine += 1;
    }
    
exit_cleanup:
    return retcode;
}

static hdsdp_retcode HUtilIGetSparseMatData( int nRow, int nCol, char *path, char *pFile, char *iFile, char *xFile,
                             int **pAp, int **pAi, double **pAx ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    char filename[100] = "";
    
    int *Ap = NULL;
    int *Ai = NULL;
    double *Ax = NULL;
    
    HDSDP_INIT(Ap, int, nCol + 1);
    
    if ( !Ap ) {
        retcode = HDSDP_RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    strcpy(filename, path);
    strcat(filename, pFile);
    HDSDP_CALL(HUtilIGetIntVecData(filename, nCol + 1, Ap));
    
    int nNz = Ap[nCol];
    HDSDP_INIT(Ai, int, nNz);
    HDSDP_INIT(Ax, double, nNz);
    
    if ( !Ai || !Ax ) {
        retcode = HDSDP_RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    strcpy(filename, path);
    strcat(filename, iFile);
    HDSDP_CALL(HUtilIGetIntVecData(filename, nNz, Ai));
    
    strcpy(filename, path);
    strcat(filename, xFile);
    HDSDP_CALL(HUtilIGetDblVecData(filename, nNz, Ax));
    
    *pAp = Ap;
    *pAi = Ai;
    *pAx = Ax;
    
exit_cleanup:
    
    if ( retcode != HDSDP_RETCODE_OK ) {
        HDSDP_FREE(Ap);
        HDSDP_FREE(Ai);
        HDSDP_FREE(Ax);
    }
    
    return retcode;
}

#define DIMENSION   ("dim.txt")
#define ABEG        ("Ap.txt")
#define AIDX        ("Ai.txt")
#define AELEM       ("Ax.txt")
#define ARHS        ("b.txt")
extern hdsdp_retcode HUtilGetSparseMatrix( char *path, int *pnRow, int *pnCol, int **pcolMatBeg,
                                           int **pColMatIdx, double **pColMatElem, double **pRhs ) {
    /* Read sparse matrix data from csv files*/
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    FILE *file = NULL;
    char filename[100] = "";
    char line[128] = "";
    
    int nRow = 0;
    int nCol = 0;
    
    int *Ap = NULL;
    int *Ai = NULL;
    double *Ax = NULL;
    double *dRhs = NULL;
    
    strcpy(filename, path);
    strcat(filename, DIMENSION);
    
    file = fopen(filename, "r");
    
    if ( !file ) {
        printf("Failed to read %s \n", filename);
        retcode = HDSDP_RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    fgets(line, LINEBUFFER, file);
    sscanf(line, "%d, %d", &nRow, &nCol);
    
    HDSDP_INIT(Ap, int, nCol + 1);
    HDSDP_INIT(dRhs, double, nCol);
    
    HDSDP_MEMCHECK(Ap);
    HDSDP_MEMCHECK(dRhs);
    
    strcpy(filename, path);
    strcat(filename, ABEG);
    HDSDP_CALL(HUtilIGetIntVecData(filename, nCol + 1, Ap));
    
    int nNz = Ap[nCol];
    HDSDP_INIT(Ai, int, nNz);
    HDSDP_INIT(Ax, double, nNz);
    HDSDP_MEMCHECK(Ai);
    HDSDP_MEMCHECK(Ax);
    
    strcpy(filename, path);
    strcat(filename, AIDX);
    HDSDP_CALL(HUtilIGetIntVecData(filename, nNz, Ai));
    
    strcpy(filename, path);
    strcat(filename, AELEM);
    HDSDP_CALL(HUtilIGetDblVecData(filename, nNz, Ax));
    
    strcpy(filename, path);
    strcat(filename, ARHS);
    HDSDP_CALL(HUtilIGetDblVecData(filename, nCol, dRhs));
    
    *pnRow = nRow;
    *pnCol = nCol;
    
    *pcolMatBeg = Ap;
    *pColMatIdx = Ai;
    *pColMatElem = Ax;
    *pRhs = dRhs;
    
exit_cleanup:
    
    if ( retcode != HDSDP_RETCODE_OK ) {
        HDSDP_FREE(Ap);
        HDSDP_FREE(Ai);
        HDSDP_FREE(Ax);
        HDSDP_FREE(dRhs);
    }
    
    return retcode;
}

extern double HUtilGetTimeStamp( void ) {
    
    return my_clock();
}

/** @brief Symmetrize an n by n matrix whose lower triangular is filled
 *
 */
extern void HUtilMatSymmetrize( int n, double *v ) {
    
    for ( int i = 0, j; i < n; ++i ) {
        for ( j = i + 1; j < n; ++j ) {
            FULL_ENTRY(v, n, i, j) = FULL_ENTRY(v, n, j, i);
        }
    }
    
    return;
}

extern void HUtilMatTranspose( int n, double *A ) {
    
    double tmp = 0.0;
    for ( int i = 0, j; i < n; ++i ) {
        for ( j = i + 1; j < n; ++j ) {
            tmp = A[i * n + j];
            A[i * n + j] = A[j * n + i];
            A[j * n + i] = tmp;
        }
    }
}

/* Debugging */
extern void HUtilPrintDblContent( int n, double *d ) {
    
    for ( int i = 0; i < n; ++i ) {
        printf("%5.3e, ", d[i]);
    }
    printf("\n");
    return;
}

extern void HUtilPrintIntContent( int n, int *d ) {
    
    for ( int i = 0; i < n; ++i ) {
        printf("%5d, ", d[i]);
    }
    printf("\n");
    return;
}

extern double HUtilGetDblMinimum( int n, double *d ) {
    
    double dMin = HDSDP_INFINITY;
    
    for ( int i = 0; i < n; ++i ) {
        dMin = HDSDP_MIN(dMin, d[i]);
    }
    
    return dMin;
}

extern double HUtilPrintDblSum( int n, double *d ) {
    
    double ds = 0.0;
    
    for ( int i = 0; i < n; ++i ) {
        ds += d[i];
    }
    
    return ds;
}

extern double HUtilPrintDblAbsSum( int n, double *d ) {
    
    double ds = 0.0;
    
    for ( int i = 0; i < n; ++i ) {
        ds += fabs(d[i]);
    }
    
    return ds;
}

extern void HUtilTraceBack( void ) {
    
    hdsdp_printf("Error Traceback \n");
    
    return;
}

/* Sorting */
extern int HUtilCheckIfAscending( int n, int *idx ) {
    /* Check is an integer array is ascending. */
    
    for ( int i = 0; i < n - 1; ++i ) {
        if ( idx[i] > idx[i + 1] ) {
            return 0;
        }
    }

    return 1;
}

extern void HUtilSortIntbyDbl( int *data, double *ref, int low, int up ) {
    
    if ( low < up ) {
        int p = dpartitioni(data, ref, low, up);
        HUtilSortIntbyDbl(data, ref, low, p - 1);
        HUtilSortIntbyDbl(data, ref, p + 1, up);
    }
    
    return;
}

extern void HUtilDescendSortIntByInt( int *data, int *ref, int low, int up ) {
    
    if ( low < up ) {
        int p = ipartitioni(data, ref, low, up);
        HUtilDescendSortIntByInt(data, ref, low, p - 1);
        HUtilDescendSortIntByInt(data, ref, p + 1, up);
    }
    
    return;
}

extern void HUtilAscendSortDblByInt( double *data, int *ref, int low, int up ) {
    
    if ( low < up ) {
        int p = ipartitiond(data, ref, low, up);
        HUtilAscendSortDblByInt(data, ref, low, p - 1);
        HUtilAscendSortDblByInt(data, ref, p + 1, up);
    }
    
    return;
}

extern void HUtilSortDblByDbl( double *data, double *ref, int low, int up ) {
    
    if ( low < up ) {
        int p = dpartitiond(data, ref, low, up);
        HUtilSortDblByDbl(data, ref, low, p - 1);
        HUtilSortDblByDbl(data, ref, p + 1, up);
    }
    
    return;
}

extern void HUtilStartCtrlCCheck( void ) {
    
    act.sa_handler = monitorCtrlC;
    sigaction(SIGINT, &act, NULL);
    
    return;
}

extern int HUtilCheckCtrlC( void ) {
    
    return isCtrlC;
}

extern void HUtilResetCtrl( void ) {
    
    isCtrlC = 0;
}

extern void HUtilWriteDblArray( char *outFileName, int nLen, double *dblContent ) {
    
    FILE *outFile = fopen(outFileName, "w");
    
    if ( !outFile ) {
        return;
    }
    
    for ( int i = 0; i < nLen; ++i ) {
        fprintf(outFile, "%20.20e,", dblContent[i]);
    }
    
    fclose(outFile);
}

hdsdp_retcode HUtilKKTCheck( void *Hkkt ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;

    hdsdp_kkt *kkt = (hdsdp_kkt *) Hkkt;
    int isValid = 1;
    
    double *kktBuffer2345 = NULL;
    double *kktBuffer3 = NULL;
    double *kktBuffer4 = NULL;
    double *kktASinvVecBuffer2345 = NULL;
    double *kktASinvVecBuffer3 = NULL;
    double *kktASinvVecBuffer4 = NULL;
    
    double *kktASinvRdSinvVecBuffer2345 = NULL;
    double *kktASinvRdSinvVecBuffer3 = NULL;
    double *kktASinvRdSinvVecBuffer4 = NULL;
    
    int nKKTNzs = 0;
    
    if ( kkt->isKKTSparse ) {
        nKKTNzs = kkt->kktMatBeg[kkt->nRow];
    } else {
        nKKTNzs = kkt->nRow * kkt->nRow;
    }
    
    HDSDP_INIT(kktBuffer2345, double, nKKTNzs);
    HDSDP_MEMCHECK(kktBuffer2345);
    
    HDSDP_INIT(kktBuffer3, double, nKKTNzs);
    HDSDP_MEMCHECK(kktBuffer3);
    
    HDSDP_INIT(kktBuffer4, double, nKKTNzs);
    HDSDP_MEMCHECK(kktBuffer4);
    
    HDSDP_INIT(kktASinvVecBuffer3, double, kkt->nRow);
    HDSDP_MEMCHECK(kktASinvVecBuffer3);
    
    HDSDP_INIT(kktASinvVecBuffer4, double, kkt->nRow);
    HDSDP_MEMCHECK(kktASinvVecBuffer4);
    
    HDSDP_INIT(kktASinvVecBuffer2345, double, kkt->nRow);
    HDSDP_MEMCHECK(kktASinvVecBuffer2345);
    
    HDSDP_INIT(kktASinvRdSinvVecBuffer3, double, kkt->nRow);
    HDSDP_MEMCHECK(kktASinvRdSinvVecBuffer3);
    
    HDSDP_INIT(kktASinvRdSinvVecBuffer4, double, kkt->nRow);
    HDSDP_MEMCHECK(kktASinvRdSinvVecBuffer4);
    
    HDSDP_INIT(kktASinvRdSinvVecBuffer2345, double, kkt->nRow);
    HDSDP_MEMCHECK(kktASinvRdSinvVecBuffer2345);
    
    /* Use hybrid strategy as a benchmark */
    printf("KKT check starts \n");
    HDSDP_CALL(HKKTBuildUpFixed(kkt, KKT_TYPE_INFEASIBLE, KKT_M3));
    HKKTExport(kkt, kktASinvVecBuffer3, kktASinvRdSinvVecBuffer3,
               NULL, NULL, NULL, NULL, NULL);
    HDSDP_MEMCPY(kktBuffer3, kkt->kktMatElem, double, nKKTNzs);
    printf("KKT check 1/3 \n");
    
    HDSDP_CALL(HKKTBuildUpFixed(kkt, KKT_TYPE_INFEASIBLE, KKT_M4));
    HKKTExport(kkt, kktASinvVecBuffer4, kktASinvRdSinvVecBuffer4,
               NULL, NULL, NULL, NULL, NULL);
    HDSDP_MEMCPY(kktBuffer4, kkt->kktMatElem, double, nKKTNzs);
    printf("KKT check 2/3 \n");
    
    HDSDP_CALL(HKKTBuildUp(kkt, KKT_TYPE_INFEASIBLE));
    HDSDP_MEMCPY(kktBuffer2345, kkt->kktMatElem, double, nKKTNzs);
    HKKTExport(kkt, kktASinvVecBuffer2345, kktASinvRdSinvVecBuffer2345,
               NULL, NULL, NULL, NULL, NULL);
    printf("KKT check 3/3 \n");
    
    double err3_4 = 0.0;
    double err3_2345 = 0.0;
    double err4_2345 = 0.0;
    
    double aerr3_4 = 0.0;
    double aerr3_2345 = 0.0;
    double aerr4_2345 = 0.0;
    
    double e3_4 = 0.0;
    double e3_2345 = 0.0;
    double e4_2345 = 0.0;
    
    for ( int iElem = 0; iElem < nKKTNzs; ++iElem ) {
        
        e3_4 = fabs(kktBuffer3[iElem] - kktBuffer4[iElem]) / (fabs(kktBuffer3[iElem]) + 1e-04);
        e3_2345 = fabs(kktBuffer3[iElem] - kktBuffer2345[iElem]) / (fabs(kktBuffer3[iElem]) + 1e-04);
        e4_2345 = fabs(kktBuffer4[iElem] - kktBuffer2345[iElem]) / (fabs(kktBuffer3[iElem]) + 1e-04);
        
        err3_4 += e3_4;
        err3_2345 += e3_2345;
        err4_2345 += e4_2345;
        
        aerr3_4 = HDSDP_MAX(aerr3_4, e3_4);
        aerr3_2345 = HDSDP_MAX(aerr3_2345, e3_2345);
        aerr4_2345 = HDSDP_MAX(aerr4_2345, e4_2345);
        
        if ( aerr3_4 >= 1e-08 || aerr3_2345 >= 1e-08 || aerr4_2345 >= 1e-08 ) {
            if ( isValid ) {
                printf("Warning. KKT consistency check failed at element %d \n", iElem);
                isValid = 0;
                retcode = HDSDP_RETCODE_FAILED;
            }
        }
    }
    
    for ( int iRow = 0; iRow < kkt->nRow; ++iRow ) {
        
        e3_4 = fabs(kktASinvVecBuffer3[iRow] - kktASinvVecBuffer4[iRow]) / (fabs(kktASinvVecBuffer3[iRow]) + 1e-04);
        e3_2345 = fabs(kktASinvVecBuffer3[iRow] - kktASinvVecBuffer2345[iRow]) / (fabs(kktASinvVecBuffer3[iRow]) + 1e-04);
        e4_2345 = fabs(kktASinvVecBuffer4[iRow] - kktASinvVecBuffer2345[iRow]) / (fabs(kktASinvVecBuffer3[iRow]) + 1e-04);
        
        err3_4 += e3_4;
        err3_2345 += e3_2345;
        err4_2345 += e4_2345;
        
        aerr3_4 = HDSDP_MAX(aerr3_4, e3_4);
        aerr3_2345 = HDSDP_MAX(aerr3_2345, e3_2345);
        aerr4_2345 = HDSDP_MAX(aerr4_2345, e4_2345);
        
        if ( aerr3_4 >= 1e-08 || aerr3_2345 >= 1e-08 || aerr4_2345 >= 1e-08 ) {
            if ( isValid ) {
                printf("Warning. KKT consistency check failed at ASinv buffer, row %d \n", iRow);
                isValid = 0;
                retcode = HDSDP_RETCODE_FAILED;
            }
        }
    }
    
    for ( int iRow = 0; iRow < kkt->nRow; ++iRow ) {
        
        e3_4 = fabs(kktASinvRdSinvVecBuffer3[iRow] - kktASinvRdSinvVecBuffer4[iRow]) / (fabs(kktASinvRdSinvVecBuffer3[iRow]) + 1e-04);
        e3_2345 = fabs(kktASinvRdSinvVecBuffer3[iRow] - kktASinvRdSinvVecBuffer2345[iRow]) / (fabs(kktASinvRdSinvVecBuffer3[iRow]) + 1e-04);
        e4_2345 = fabs(kktASinvRdSinvVecBuffer4[iRow] - kktASinvRdSinvVecBuffer2345[iRow]) / (fabs(kktASinvRdSinvVecBuffer3[iRow]) + 1e-04);
        
        err3_4 += e3_4;
        err3_2345 += e3_2345;
        err4_2345 += e4_2345;
        
        aerr3_4 = HDSDP_MAX(aerr3_4, e3_4);
        aerr3_2345 = HDSDP_MAX(aerr3_2345, e3_2345);
        aerr4_2345 = HDSDP_MAX(aerr4_2345, e4_2345);
        
        if ( aerr3_4 >= 1e-08 || aerr3_2345 >= 1e-08 || aerr4_2345 >= 1e-08 ) {
            if ( isValid ) {
                printf("Warning. KKT consistency check failed at ASinvRdSinv buffer, row %d \n", iRow);
                isValid = 0;
                retcode = HDSDP_RETCODE_FAILED;
            }
        }
    }
    
    
    printf("KKT consistency check: | 3-4: %6.3e | 3-2345: %6.3e | 4-2345: %6.3e \n",
           aerr3_4, aerr3_2345, aerr4_2345);
    
exit_cleanup:
    
    HDSDP_FREE(kktBuffer2345);
    HDSDP_FREE(kktBuffer3);
    HDSDP_FREE(kktBuffer4);
    HDSDP_FREE(kktASinvVecBuffer2345);
    HDSDP_FREE(kktASinvVecBuffer3);
    HDSDP_FREE(kktASinvVecBuffer4);
    HDSDP_FREE(kktASinvRdSinvVecBuffer2345);
    HDSDP_FREE(kktASinvRdSinvVecBuffer3);
    HDSDP_FREE(kktASinvRdSinvVecBuffer4);
    
    return retcode;
}

#ifdef LINSYS_PARDISO
extern int MKL_Get_Max_Threads( void );
extern int MKL_Set_Num_Threads( int nth );
#else
extern int MKL_Get_Max_Threads( void ) {
    return 1;
}
extern int MKL_Set_Num_Threads( int nth ) {
    return 0;
}
#endif

extern int HUtilGetGlobalMKLThreads( void ) {
    
    return MKL_Get_Max_Threads();
}

extern void HUtilSetGlobalMKLThreads( int nTargetThreads ) {
    
    MKL_Set_Num_Threads(nTargetThreads);
    
    return;
}
