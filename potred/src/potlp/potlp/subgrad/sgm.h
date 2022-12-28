/** @file sgm.h
 *  @brief Implement an LP solver based on sharpness and subgradient
 *
 * */
#ifndef sgm_h
#define sgm_h

#include <stdio.h>

#define SGM_RETCODE_OK      0
#define SGM_RETCODE_FAILED  1

#define SGM_INFINITY 1e+30
extern int potSGMSolve( int nCol, int nEqRow, int nIneqRow, int *eqMatBeg, int *eqMatIdx, double *eqMatElem,
                        int *inEqMatBeg, int *inEqMatIdx, double *inEqMatElem, double *eqMatRhs, double *ineqMatRhs,
                        double *colObj, double *colBoundUp, double *colVal, double *rowDual, double *colSlack,
                        double *bdSlack, int maxIter, double relTol, double stepScal );

#endif /* sgm_h */
