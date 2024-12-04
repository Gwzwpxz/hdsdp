#ifndef lp_mps_h
#define lp_mps_h

#ifdef HEADERPATH
#include "interface/hdsdp.h"
#else
#include "hdsdp.h"
#endif

/* Implement an LP mps file reader */
int potLpMpsRead( char *fname, char *name, int *pnRow, int *pnEqRow, int *pnInEqRow, int *pnCol, int *pnElem,
                      int **peqMatBeg,  int **peqMatIdx, double **peqMatElem, int **pIneqMatBeg,
                      int **pIneqMatIdx, double **pIneqMatElem, double **prowRHS, double **pcolObj,
                      int *pnColUb, int **pcolUbIdx, double **pcolUbElem );

#endif /* lp_mps_h */
