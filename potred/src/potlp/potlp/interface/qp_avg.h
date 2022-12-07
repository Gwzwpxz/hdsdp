#ifndef qp_avg_h
#define qp_avg_h

#include "pot_def.h"
#include "pot_param.h"

/** @struct pot\_qpsolver
 * An internal QP driver using the interior point method
 *
 */
typedef struct {
    
    int nRowAll;
    int nRow;
    int nQuadCol;
    
    double *QMatElem;
    double *colMatElem;
    
    /* Internal iterator */
    double *s;
    double *ds;
    
    double *lbd;
    double *dlbd;
    
    double *alpha;
    double *dalpha;
    double *alphabest;
    
    double *M;
    double *SLV;
    double *d1;
    double *d2;
    double *Valpha;
    double *buffer;
    
    double pObjVal;
    double *xWindow;
    
} pot_qpsolver;

#ifdef __cplusplus
extern "C" {
#endif

extern pot_int potQPCreate( pot_qpsolver **ppotQP );
extern pot_int potQPInit( pot_qpsolver *potQP, int nQuadCol, int nRowAll, int nRow );
extern void potQPLoadProb( pot_qpsolver *potQP, int nIter, int nCol, int nRow,
                           int nCone, double *gradWindow, double *xVarWindow );
extern pot_int potQPSolveProb( pot_qpsolver *potQP, double targetObj, double *xAvg );
extern void potQPClear( pot_qpsolver *potQP );
extern void potQPDestroy( pot_qpsolver **ppotQP );

#ifdef __cplusplus
}
#endif

#endif /* qp_avg_h */
