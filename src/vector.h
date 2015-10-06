#ifndef vector_h
#define vector_h

#include "matrix.h"

typedef mat1d_dreal  vect_dreal;
typedef mat1d_dcmplx vect_dcmplx;

dreal  vect_Dotprod_dreal(const vect_dreal *, const vect_dreal *);
dcmplx vect_Dotprod_dcmplx(const vect_dcmplx *, const vect_dcmplx *);
dreal  vect_Norm_dreal(const vect_dreal *);
dreal  vect_Norm_dcmplx(const vect_dcmplx *);
void   vect_Renorm_dreal(vect_dreal *);
void   vect_Renorm_dcmplx(vect_dcmplx *);
void   vect_Orthog_dreal(int, vect_dreal * const * const, vect_dreal *);
void   vect_Orthog_dcmplx(int, vect_dcmplx * const * const, vect_dcmplx *);
void   vect_RandNorm_dreal(vect_dreal *);
void   vect_RandNorm_dcmplx(vect_dcmplx *);

#define vect_Reset_dreal  mat1d_Reset_dreal
#define vect_Reset_dcmplx mat1d_Reset_dcmplx
#define vect_New_dreal    mat1d_New_dreal
#define vect_New_dcmplx   mat1d_New_dcmplx
#define vect_Del_dreal    mat1d_Del_dreal
#define vect_Del_dcmplx   mat1d_Del_dcmplx
#define vect_Copy_dreal   mat1d_Copy_dreal
#define vect_Copy_dcmplx  mat1d_Copy_dcmplx

#endif
