#ifndef conjgrad_h
#define conjgrad_h

#include "vector.h"

typedef void (*Operator)(const vect_dreal *, vect_dreal *);

void conjgrad_Posdef(int, const vect_dreal *, vect_dreal *, Operator);
void conjgrad_Symmtr(int, const vect_dreal *, vect_dreal *, Operator);
void conjgrad_Check(void);

#endif
