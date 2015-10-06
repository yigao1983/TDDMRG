#ifndef parameters_h
#define parameters_h

#include "precision.h"
#include "matrix.h"
#include "vector.h"

extern int   Nspin, Length, Bridge;
extern dreal *Echain, *Uchain, *Tchain;

void parameters_Init(void);
void parameters_Quit(void);

#endif
