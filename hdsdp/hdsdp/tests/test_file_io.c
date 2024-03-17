#include <stdio.h>

#ifdef HEADERPATH
#include "interface/hdsdp.h"
#include "interface/hdsdp_utils.h"
#include "interface/hdsdp_user_data.h"
#include "interface/hdsdp_file_io.h"
#include "interface/hdsdp_conic.h"
#include "interface/hdsdp_schur.h"
#else
#include "hdsdp.h"
#include "hdsdp_utils.h"
#include "hdsdp_user_data.h"
#include "hdsdp_file_io.h"
#include "hdsdp_conic.h"
#include "hdsdp_schur.h"
#endif

#include <math.h>

#define LINEBUFFER  (1024)
static int getIntVecData( char *fname, int nElem, int *iData ) {
    
    int retcode = HDSDP_RETCODE_OK;
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

static int getDblVecData( char *fname, int nElem, double *dData ) {
    
    int retcode = HDSDP_RETCODE_OK;
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

static int getSparseMatData( int nRow, int nCol, char *path, char *pFile, char *iFile, char *xFile,
                             int **pAp, int **pAi, double **pAx ) {
    
    int retcode = HDSDP_RETCODE_OK;
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
    HDSDP_CALL(getIntVecData(filename, nCol + 1, Ap));
    
    int nNz = Ap[nCol];
    HDSDP_INIT(Ai, int, nNz);
    HDSDP_INIT(Ax, double, nNz);
    
    if ( !Ai || !Ax ) {
        retcode = HDSDP_RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    strcpy(filename, path);
    strcat(filename, iFile);
    HDSDP_CALL(getIntVecData(filename, nNz, Ai));
    
    strcpy(filename, path);
    strcat(filename, xFile);
    HDSDP_CALL(getDblVecData(filename, nNz, Ax));
    
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

int test_file_io( char *fname ) {
    
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    int nConstrs = 0;
    int nBlks = 0;
    int *BlkDims = NULL;
    double *rowRHS = NULL;
    int **coneMatBeg = NULL;
    int **coneMatIdx = NULL;
    double **coneMatElem = NULL;
    int nLpCols = 0;
    int *LpMatBeg = NULL;
    int *LpMatIdx = NULL;
    double *LpMatElem = NULL;
    int nCols = 0;
    int nElem = 0;
    double *rowDual = NULL;
    double *rowDualStep = NULL;
    double logdet = 0.0;
    
    double *kktLhsBuffer = NULL;
    
    hdsdp_cone **SDPCones = NULL;
    user_data *SDPData = NULL;
    hdsdp_cone *SDPCone = NULL;
    hdsdp_kkt *kkt = NULL;
    
    double timeStart = HUtilGetTimeStamp();
    printf("Filename: %s\n", fname);
    
    HDSDP_CALL(HReadSDPA(fname, &nConstrs, &nBlks, &BlkDims, &rowRHS, &coneMatBeg,
                         &coneMatIdx, &coneMatElem, &nCols, &nLpCols, &LpMatBeg,
                         &LpMatIdx, &LpMatElem, &nElem));
    
    printf("Reading SDPA file in %f seconds \n", HUtilGetTimeStamp() - timeStart);
    
    HDSDP_CALL(HUserDataCreate(&SDPData));
    
    HDSDP_INIT(rowDual, double, nConstrs);
    HDSDP_MEMCHECK(rowDual);
    
    HDSDP_INIT(rowDualStep, double, nConstrs);
    HDSDP_MEMCHECK(rowDualStep);
    
    HDSDP_INIT(kktLhsBuffer, double, nConstrs);
    HDSDP_MEMCHECK(kktLhsBuffer);
    
    HDSDP_INIT(SDPCones, hdsdp_cone *, nBlks);
    HDSDP_MEMCHECK(SDPCones);
    
    for ( int iBlk = 0; iBlk < nBlks; ++iBlk ) {
        HUserDataSetConeData(SDPData, HDSDP_CONETYPE_DENSE_SDP, nConstrs, BlkDims[iBlk],
                             coneMatBeg[iBlk], coneMatIdx[iBlk], coneMatElem[iBlk]);
        cone_type cone = HUserDataChooseCone(SDPData);
        HDSDP_CALL(HConeCreate(&SDPCone));
        
        SDPCones[iBlk] = SDPCone;
        
        HDSDP_CALL(HConeSetData(SDPCone, SDPData));
        HDSDP_CALL(HConeProcData(SDPCone));
//        HConeView(SDPCone);
        HDSDP_CALL(HConePresolveData(SDPCone));
        HConeView(SDPCone);
        
        for ( int i = 0; i < nConstrs; ++i ) {
            rowDual[i] = 0.0 * (double) (i + 1) / nConstrs;
            rowDualStep[i] = (double) (i + 1) / nConstrs;
        }
        
        HConeSetStart(SDPCone, -1e+03);
        HConeUpdate(SDPCone, 1.0, rowDual);
//        HConeView(SDPCone);
        
        HDSDP_CALL(HConeGetLogBarrier(SDPCone, 1.5, rowDual, BUFFER_DUALVAR, &logdet));
//        printf("- Conic log det (S) = %e. \n", logdet);
        
        double ratio = 0.0;
        HDSDP_CALL(HConeRatioTest(SDPCone, 1.0, rowDualStep, 0.9, BUFFER_DUALVAR, &ratio));
        printf("- Conic ratio test: %e. \n", ratio);
        HDSDP_CALL(HConeRatioTest(SDPCone, 1.0, rowDualStep, 0.9, BUFFER_DUALVAR, &ratio));
        printf("- Conic ratio test again: %e. \n", ratio);
        
        HUserDataClear(SDPData);
    }
    
    /* KKT setup */
    HDSDP_CALL(HKKTCreate(&kkt));
    HDSDP_CALL(HKKTInit(kkt, nConstrs, nBlks, SDPCones));
    
    HDSDP_CALL(HKKTBuildUp(kkt, KKT_TYPE_INFEASIBLE));
//    HDSDP_PROFILER(HKKTBuildUp(kkt, KKT_TYPE_INFEASIBLE), 100);
    
    /* KKT solve */
    double dCSinv = 0.0;
    double dCSinvRdCSinv = 0.0;
    double dCSinvCSinv = 0.0;
    
    HKKTExport(kkt, kktLhsBuffer, NULL, NULL, &dCSinvCSinv, &dCSinv, &dCSinvRdCSinv, NULL);
    HDSDP_CALL(HKKTFactorize(kkt));
    HDSDP_CALL(HKKTSolve(kkt, kktLhsBuffer, NULL));
    
    HKKTExport(kkt, NULL, kktLhsBuffer, NULL, &dCSinvCSinv, &dCSinv, &dCSinvRdCSinv, NULL);
    HDSDP_CALL(HKKTFactorize(kkt));
    HDSDP_CALL(HKKTSolve(kkt, kktLhsBuffer, NULL));
    
    HKKTExport(kkt, NULL, NULL, kktLhsBuffer, &dCSinvCSinv, &dCSinv, &dCSinvRdCSinv, NULL);
    HDSDP_CALL(HKKTFactorize(kkt));
    HDSDP_CALL(HKKTSolve(kkt, kktLhsBuffer, NULL));
    
    /* KKT consistency */
    HDSDP_CALL(HUtilKKTCheck(kkt));
    
exit_cleanup:
    
    HDSDP_FREE(kktLhsBuffer);
    
    HUserDataDestroy(&SDPData);
    HKKTDestroy(&kkt);
    
    for ( int iBlk = 0; iBlk < nBlks; ++iBlk ) {
        SDPCone = SDPCones[iBlk];
        HConeDestroy(&SDPCone);
    }
    
    HDSDP_FREE(BlkDims);
    HDSDP_FREE(rowRHS);
    HDSDP_FREE(rowDual);
    HDSDP_FREE(rowDualStep);
    HDSDP_FREE(SDPCones);
    
    for ( int iBlk = 0; iBlk < nBlks; ++iBlk ) {
        HDSDP_FREE(coneMatBeg[iBlk]);
        HDSDP_FREE(coneMatIdx[iBlk]);
        HDSDP_FREE(coneMatElem[iBlk]);
    }
    
    HDSDP_FREE(coneMatBeg);
    HDSDP_FREE(coneMatIdx);
    HDSDP_FREE(coneMatElem);
    
    if ( nLpCols > 0 ) {
        HDSDP_FREE(LpMatBeg);
        HDSDP_FREE(LpMatIdx);
        HDSDP_FREE(LpMatElem);
    }
    
    return retcode;
}

int test_solver( char *fname ) {
    
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    int nConstrs = 0;
    int nCones = 0;
    int *BlkDims = NULL;
    double *rowRHS = NULL;
    int **coneMatBeg = NULL;
    int **coneMatIdx = NULL;
    double **coneMatElem = NULL;
    int nLpCols = 0;
    int *LpMatBeg = NULL;
    int *LpMatIdx = NULL;
    double *LpMatElem = NULL;
    int nCols = 0;
    int nElem = 0;
    
    hdsdp *hsolve = NULL;
    user_data **SDPDatas = NULL;
    user_data *SDPData = NULL;
    user_data *LPData = NULL;
    
    double timeStart = HUtilGetTimeStamp();
    printf("Filename: %s\n", fname);
    
    HDSDP_CALL(HReadSDPA(fname, &nConstrs, &nCones, &BlkDims, &rowRHS, &coneMatBeg,
                         &coneMatIdx, &coneMatElem, &nCols, &nLpCols, &LpMatBeg,
                         &LpMatIdx, &LpMatElem, &nElem));
    
    printf("Reading SDPA file in %f seconds \n", HUtilGetTimeStamp() - timeStart);
    
    HDSDP_INIT(SDPDatas, user_data *, nCones);
    HDSDP_MEMCHECK(SDPDatas);
    
    HDSDP_CALL(HDSDPCreate(&hsolve));
    
    if ( nLpCols > 0 ) {
        HDSDP_CALL(HDSDPInit(hsolve, nConstrs, nCones + 1));
    } else {
        HDSDP_CALL(HDSDPInit(hsolve, nConstrs, nCones));
    }
    
    for ( int iCone = 0; iCone < nCones; ++iCone ) {
        HDSDP_CALL(HUserDataCreate(&SDPData));
        HUserDataSetConeData(SDPData, HDSDP_CONETYPE_DENSE_SDP, nConstrs, BlkDims[iCone],
                             coneMatBeg[iCone], coneMatIdx[iCone], coneMatElem[iCone]);
        HDSDP_CALL(HDSDPSetCone(hsolve, iCone, SDPData));
        SDPDatas[iCone] = SDPData;
        SDPData = NULL;
    }
    
    if ( nLpCols > 0 ) {
        HDSDP_CALL(HUserDataCreate(&LPData));
        HUserDataSetConeData(LPData, HDSDP_CONETYPE_LP, nConstrs, nLpCols, LpMatBeg, LpMatIdx, LpMatElem);
        HDSDP_CALL(HDSDPSetCone(hsolve, nCones, LPData));
    }
    
    HDSDPSetDualObjective(hsolve, rowRHS);
    HDSDP_CALL(HDSDPOptimize(hsolve, 1));
    
exit_cleanup:
    
    for ( int iBlk = 0; iBlk < nCones; ++iBlk ) {
        HUserDataDestroy(&SDPDatas[iBlk]);
    }
    HDSDP_FREE(SDPDatas);
    
    HUserDataDestroy(&SDPData);
    HUserDataDestroy(&LPData);
    
    HDSDPDestroy(&hsolve);
    
    HDSDP_FREE(BlkDims);
    HDSDP_FREE(rowRHS);
    
    for ( int iBlk = 0; iBlk < nCones; ++iBlk ) {
        HDSDP_FREE(coneMatBeg[iBlk]);
        HDSDP_FREE(coneMatIdx[iBlk]);
        HDSDP_FREE(coneMatElem[iBlk]);
    }
    
    HDSDP_FREE(coneMatBeg);
    HDSDP_FREE(coneMatIdx);
    HDSDP_FREE(coneMatElem);
    
    if ( nLpCols > 0 ) {
        HDSDP_FREE(LpMatBeg);
        HDSDP_FREE(LpMatIdx);
        HDSDP_FREE(LpMatElem);
    }
    
    return retcode;
}


#define DIMENSION ("dim.csv")
#define ABEG      ("Ap.csv")
#define AIDX      ("Ai.csv")
#define AELEM     ("Ax.csv")
#define NSLV      1000

//#define GPUTEST
#ifdef GPUTEST
#include "cuDSS.h"
#define CUDA_CALL_AND_CHECK(call, msg) \
    do { \
        cuda_error = call; \
        if (cuda_error != cudaSuccess) { \
            printf("Example FAILED: CUDA API returned error = %d, details: " #msg "\n", cuda_error); \
            CUDSS_EXAMPLE_FREE; \
            return -1; \
        } \
    } while(0);

#define CUDSS_CALL_AND_CHECK(call, status, msg) \
    do { \
        status = call; \
        if (status != CUDSS_STATUS_SUCCESS) { \
            printf("Example FAILED: CUDSS call ended unsuccessfully with status = %d, details: " #msg "\n", status); \
            CUDSS_EXAMPLE_FREE; \
            return -2; \
        } \
    } while(0);
#endif


int test_mat( char *path ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    FILE *file = NULL;
    char filename[100] = "";
    char line[128] = "";
    
    int nCol = 0;
    int nElem = 0;
    int *AMatBeg = NULL;
    int *AMatIdx = NULL;
    double *AMatElem = NULL;
    double *rhs = NULL;
    double *sol = NULL;
    
#ifdef GPUTEST
    int *AMatBegCuda = NULL;
    int *AMatIdxCuda = NULL;
    double *AMatElemCuda = N.ULL;
    double *rhsCuda = NULL;
    double *solCuda = NULL;
#endif
    
    hdsdp_linsys *cpuLinSolver = NULL;
    
    /* Get dimension of problem */
    strcpy(filename, path);
    strcat(filename, DIMENSION);
    
    file = fopen(filename, "r");
    
    if ( !file ) {
        printf("Failed to read %s \n", filename);
        retcode = HDSDP_RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    fgets(line, LINEBUFFER, file);
    sscanf(line, "%d", &nCol);
    
    /* Get data */
    HDSDP_CALL(getSparseMatData(nCol, nCol, path, ABEG, AIDX, AELEM, &AMatBeg, &AMatIdx, &AMatElem));
    nElem = AMatBeg[nCol];
    
    HDSDP_INIT(rhs, double, nCol);
    HDSDP_INIT(sol, double, nCol);
    
    for ( int iElem = 0; iElem < nCol; ++iElem ) {
        rhs[iElem] = 1.0 / (iElem + 1);
    }
    
    HDSDP_CALL(HFpLinsysCreate(&cpuLinSolver, nCol, HDSDP_LINSYS_SPARSE_DIRECT));
    HFpLinsysSetParam(cpuLinSolver, -1.0, -1.0, 1, -1, -1);
    HDSDP_CALL(HFpLinsysSymbolic(cpuLinSolver, AMatBeg, AMatIdx));
    
    double tStart = HUtilGetTimeStamp();
    double dCpuTime = 0.0;
    
    /* Test CPU solve */
    for ( int iTest = 0; iTest < NSLV; ++iTest ) {
        HDSDP_CALL(HFpLinsysNumeric(cpuLinSolver, AMatBeg, AMatIdx, AMatElem));
    }
    dCpuTime = HUtilGetTimeStamp() - tStart;
    
    HDSDP_CALL(HFpLinsysSolve(cpuLinSolver, 1, rhs, sol));
    
    printf("[CPU Test] n = %d | %d factors | %f seconds \n", nCol, NSLV, dCpuTime);
    
#ifdef GPUTEST
    
    cudaError_t cuda_error = cudaSuccess;
    cudssStatus_t status = CUDSS_STATUS_SUCCESS;
    /* Test GPU solve */
    /* Allocate device memory for A, x and b */
    CUDA_CALL_AND_CHECK(cudaMalloc(&AMatBegCuda, (nCol + 1) * sizeof(int)), "cudaMalloc for AMatBegCuda");
    CUDA_CALL_AND_CHECK(cudaMalloc(&AMatIdxCuda, nElem * sizeof(int)), "cudaMalloc for AMatIdxCuda");
    CUDA_CALL_AND_CHECK(cudaMalloc(&AMatElem, nElem * sizeof(double)), "cudaMalloc for AMatElemCuda");
    CUDA_CALL_AND_CHECK(cudaMalloc(&rhsCuda, nElem * sizeof(double)), "cudaMalloc for rhsCuda");
    CUDA_CALL_AND_CHECK(cudaMalloc(&solCuda, nElem * sizeof(double)), "cudaMalloc for solCuda");
    
    /* Copy host memory to device for A and b */
    CUDA_CALL_AND_CHECK(cudaMemcpy(AMatBegCuda, AMatBeg, (nCol + 1) * sizeof(int),
                            cudaMemcpyHostToDevice), "cudaMemcpy for AMatBeg");
    CUDA_CALL_AND_CHECK(cudaMemcpy(AMatIdxCuda, AMatIdx, nElem * sizeof(int),
                            cudaMemcpyHostToDevice), "cudaMemcpy for AMatIdx");
    CUDA_CALL_AND_CHECK(cudaMemcpy(AMatElemCuda, AMatElem, nElem * sizeof(double),
                            cudaMemcpyHostToDevice), "cudaMemcpy for AMatElem");
    
    /* Create a CUDA stream */
    cudaStream_t stream = NULL;
    CUDA_CALL_AND_CHECK(cudaStreamCreate(&stream), "cudaStreamCreate");

    /* Creating the cuDSS library handle */
    cudssHandle_t handle;
    CUDSS_CALL_AND_CHECK(cudssCreate(&handle), status, "cudssCreate");
    
    /* Creating cuDSS solver configuration and data objects */
    cudssConfig_t solverConfig;
    cudssData_t solverData;
    CUDSS_CALL_AND_CHECK(cudssConfigCreate(&solverConfig), status, "cudssConfigCreate");
    CUDSS_CALL_AND_CHECK(cudssDataCreate(handle, &solverData), status, "cudssDataCreate");
    
    /* Create matrix objects for the right-hand side b and solution x (as dense matrices). */
    cudssMatrix_t x, b;
    int64_t nrows = nCol, ncols = nCol;
    int ldb = ncols, ldx = nrows;
    CUDSS_CALL_AND_CHECK(cudssMatrixCreateDn(&b, ncols, 1, ldb, rhsCuda, CUDA_R_64F,
                        CUDSS_LAYOUT_COL_MAJOR), status, "cudssMatrixCreateDn for b");
    CUDSS_CALL_AND_CHECK(cudssMatrixCreateDn(&x, nrows, 1, ldx, solCuda, CUDA_R_64F,
                         CUDSS_LAYOUT_COL_MAJOR), status, "cudssMatrixCreateDn for x");
    
    /* Create a matrix object for the sparse input matrix. */
    cudssMatrix_t A;
    cudssMatrixType_t mtype     = CUDSS_MTYPE_SPD;
    cudssMatrixViewType_t mview = CUDSS_MVIEW_UPPER;
    cudssIndexBase_t base       = CUDSS_BASE_ZERO;
    
    CUDSS_CALL_AND_CHECK(cudssMatrixCreateCsr(&A, (int64_t) nCol, (int64_t) nCol, nElem, AMatBegCuda, NULL,
                         AMatIdxCuda, AMatElemCuda, CUDA_R_32I, CUDA_R_64F, mtype, mview,
                         base), status, "cudssMatrixCreateCsr");
    
    /* Symbolic factorization */
    CUDSS_CALL_AND_CHECK(cudssExecute(handle, CUDSS_PHASE_ANALYSIS, solverConfig, solverData,
                        A, x, b), status, "cudssExecute for analysis");
    
    /* Test factorization */
    tStart = HUtilGetTimeStamp();
    double dGpuTime = 0.0;
    for ( int iTest = 0; iTest < NSLV; ++iTest ) {
        CUDSS_CALL_AND_CHECK(cudssExecute(handle, CUDSS_PHASE_FACTORIZATION, solverConfig,
                                 solverData, A, x, b), status, "cudssExecute for factor");
    }
    dGpuTime = HUtilGetTimeStamp() - tStart;
    printf("[NVIDIA Test] n = %d | %d factors | %f seconds \n", nCol, NSLV, dCpuTime);
    
    CUDSS_CALL_AND_CHECK(cudssMatrixDestroy(A), status, "cudssMatrixDestroy for A");
    CUDSS_CALL_AND_CHECK(cudssMatrixDestroy(b), status, "cudssMatrixDestroy for b");
    CUDSS_CALL_AND_CHECK(cudssMatrixDestroy(x), status, "cudssMatrixDestroy for x");
    CUDSS_CALL_AND_CHECK(cudssDataDestroy(handle, solverData), status, "cudssDataDestroy");
    CUDSS_CALL_AND_CHECK(cudssConfigDestroy(solverConfig), status, "cudssConfigDestroy");
    CUDSS_CALL_AND_CHECK(cudssDestroy(handle), status, "cudssHandleDestroy");
    CUDA_CALL_AND_CHECK(cudaStreamSynchronize(stream), "cudaStreamSynchronize");
    
#endif
    
exit_cleanup:
    
    HDSDP_FREE(sol);
    HDSDP_FREE(rhs);
    HDSDP_FREE(AMatBeg);
    HDSDP_FREE(AMatIdx);
    HDSDP_FREE(AMatElem);
    
#ifdef GPUTEST
    cudaFree(solCuda);
    cudaFree(rhsCuda);
    cudaFree(AMatBegCuda);
    cudaFree(AMatIdxCuda);
    cudaFree(AMatElemCuda);
#endif
    
    HFpLinsysDestroy(&cpuLinSolver);
    return retcode;
}
