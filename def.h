// def.h
#ifndef _DEF
#define _DEF

#include "myparameters.h"

#define NULL 0
#define FALSE 0
#define TRUE 1
typedef int BOOL;

#define NV  ((par.NVX)*par.NVY)
#define NNX ((par.NVX)+1)
#define NNY (par.NVY+1)
#define NN  (NNX*NNY)
#define NDOF (2*NN)


#define JCM ((par.NOSTICKJCM)*(par.VOXSIZE))  // cell1-medium
#define JCC ((par.NOSTICKJCC)*(par.VOXSIZE)) // cell1-cell1
#define JCM2 ((par.NOSTICKJCM2)*(par.VOXSIZE))  // cell2-medium
#define JCC2 ((par.NOSTICKJCC2)*(par.VOXSIZE)) // cell2-cell2
#define JCCb ((par.NOSTICKJCCb)*(par.VOXSIZE)) // cell1-cell2

#define SQ05 sqrt(0.5)
#endif
