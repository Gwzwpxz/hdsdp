#ifndef linsys_h
#define linsys_h

#include "pot_def.h"

typedef struct {
    
    pot_int nCol;
    void *solver;
    
    int  backUpLin;
    
    pot_int  (*LCreate)  ( void **, pot_int );
    void (*LDestroy) ( void ** );
    
    pot_int  (*LSFac)  ( void *, pot_int *, pot_int * );
    pot_int  (*LNFac)  ( void *, pot_int *, pot_int *, double * );
    pot_int  (*LNFacBackup)  ( void *, pot_int *, pot_int *, double * );
    pot_int  (*LSolve) ( void *, double * );

} pot_linsys;

extern pot_int potLinsysCreate( pot_linsys **ppotLinsys );
extern pot_int potLinsysInit( pot_linsys *potLinsys, pot_int nCol );
extern pot_int potLinsysSymFactorize( pot_linsys *potLinsys, pot_int *colMatBeg, pot_int *colMatIdx );
extern pot_int potLinsysNumFactorize( pot_linsys *potLinsys, pot_int *colMatBeg, pot_int *colMatIdx, double *colMatElem );
extern pot_int potLinsysSolve( pot_linsys *potLinsys, double *rhsVec, double *solVec );
extern void potLinsysSwitchToBackup( pot_linsys *potLinsys );
extern void potLinsysClear( pot_linsys *potLinsys );
extern void potLinsysDestroy( pot_linsys **ppotLinsys );

#endif /* linsys_h */
