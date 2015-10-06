#ifndef fourier_h
#define fourier_h

#include "precision.h"

extern const int fourier_forward;
extern const int fourier_backward;

void fourier_dcmplx_1d(int, dcmplx *, dcmplx *, int);
void fourier_dcmplx_3d(int, int, int, dcmplx *, dcmplx *, int);

#endif
