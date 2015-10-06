#ifndef davidson_h
#define davidson_h

#include "vector.h"
#include "conjgrad.h"

typedef void (*Precondt)(dreal rho, const vect_dreal *vect, vect_dreal *Pvect);

void davidson_Diag(int, dreal *, vect_dreal *, Operator, Precondt);
//void davidson_Check(void);

#endif
