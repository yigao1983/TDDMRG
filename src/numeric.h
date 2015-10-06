#ifndef numeric_h
#define numeric_h

#include "precision.h"

dreal  numeric_Maxval_dreal(int, const dreal *);
dreal  numeric_Minval_dreal(int, const dreal *);
dreal  numeric_Quadrat_dreal(int, const dreal *, const dreal *);
dcmplx numeric_Quadrat_dcmplx(int, const dreal *, const dcmplx *);
int    numeric_Ipow(int, int);
int    numeric_Ispow2(int);
void   numeric_Kronig(int, const dcmplx *, dcmplx *);
void   numeric_Convolut(int, dreal, int, const dcmplx *, int, const dcmplx *, dcmplx *);

#endif
