#include <stdlib.h>
#include <math.h>
#include "dbg.h"
#include "randnum.h"
#include "vector.h"

dreal vect_Dotprod_dreal(const vect_dreal *vect_a, const vect_dreal *vect_b)
{
  int   ndim_a, ndim_b, idim;
  dreal prod;
  
  check_mem(vect_a, "vect_a");
  check_mem(vect_b, "vect_b");
  ndim_a = vect_a->ndim;
  ndim_b = vect_b->ndim;
  check(ndim_a == ndim_b, "Inequivalent ndim");
  
  for(idim = 0, prod = 0.0; idim < ndim_a; ++idim)
    prod += vect_a->addr[idim] * vect_b->addr[idim];
  
  return prod;
  
 error:
  abort();
}

dcmplx vect_Dotprod_dcmplx(const vect_dcmplx *vect_a, const vect_dcmplx *vect_b)
{
  int    ndim_a, ndim_b, idim;
  dcmplx prod;
 
  check_mem(vect_a, "vect_a");
  check_mem(vect_b, "vect_b");
  ndim_a = vect_a->ndim;
  ndim_b = vect_b->ndim;
  check(ndim_a == ndim_b, "Inequivalent ndim");
  
  for(idim = 0, prod = 0.0; idim < ndim_a; ++idim)
    prod += conj(vect_a->addr[idim]) * vect_b->addr[idim];
  
  return prod;
  
 error:
  abort();
}

dreal vect_Norm_dreal(const vect_dreal *vect)
{
  dreal prod;
  
  prod = vect_Dotprod_dreal(vect, vect);
  check(prod >= 0.0, "Invalid prod");
  
  return sqrt(prod);
  
 error:
  abort();
}

dreal vect_Norm_dcmplx(const vect_dcmplx *vect)
{
  dcmplx prod;
  
  prod = vect_Dotprod_dcmplx(vect, vect);
  check(creal(prod) >= 0.0, "Invalid prod");
  
  return sqrt(creal(prod));
  
 error:
  abort();
}

void vect_Renorm_dreal(vect_dreal *vect)
{
  int   ndim, idim;
  dreal norm;
  
  check_mem(vect, "vect");
  
  ndim = vect->ndim;
  norm = vect_Norm_dreal(vect);
  
  check(norm != 0.0, "Invalid norm");
  
  for(idim = 0; idim < ndim; ++idim)
    vect->addr[idim] /= norm;
  
  return;
  
 error:
  abort();
}

void vect_Renorm_dcmplx(vect_dcmplx *vect)
{
  int   ndim, idim;
  dreal norm;
  
  check_mem(vect, "vect");
  
  ndim = vect->ndim;
  norm = vect_Norm_dcmplx(vect);
  
  check(norm != 0.0, "Invalid norm");
  
  for(idim = 0; idim < ndim; ++idim)
    vect->addr[idim] /= norm;
  
  return;
  
 error:
  abort();
}

void vect_Orthog_dreal(int nsub, vect_dreal * const * const vect_arr, vect_dreal *vect)
{
  int        isub, ndim, idim;
  dreal      prod;
  vect_dreal *vect_sub = NULL, *vect_sav = NULL;
  
  check_mem(vect, "vect");
  ndim = vect->ndim;
  
  vect_sav = vect_New_dreal();
  check_mem(vect_sav, "vect_sav");
  vect_sav->Alloc(ndim, vect_sav);
  vect_Copy_dreal(vect, vect_sav);
  
  for(isub = 0; isub < nsub; ++isub) {
    vect_sub = vect_arr[isub];
    prod = vect_Dotprod_dreal(vect_sub, vect_sav);
    for(idim = 0; idim < ndim; ++idim)
      vect->addr[idim] -= prod * vect_sub->addr[idim];
  }
  
  vect_Del_dreal(vect_sav);
  
  return;
  
 error:
  abort();
}

void vect_Orthog_dcmplx(int nsub, vect_dcmplx * const * const vect_arr, vect_dcmplx *vect)
{
  int         isub, ndim, idim;
  dreal       prod;
  vect_dcmplx *vect_sub = NULL, *vect_sav = NULL;
  
  check_mem(vect, "vect");
  ndim = vect->ndim;
  
  vect_sav = vect_New_dcmplx();
  check_mem(vect_sav, "vect_sav");
  vect_sav->Alloc(ndim, vect_sav);
  vect_Copy_dcmplx(vect, vect_sav);
  
  for(isub = 0; isub < nsub; ++isub) {
    vect_sub = vect_arr[isub];
    prod = vect_Dotprod_dcmplx(vect_sub, vect_sav);
    for(idim = 0; idim < ndim; ++idim)
      vect->addr[idim] -= prod * vect_sub->addr[idim];
  }
  
  vect_Del_dcmplx(vect_sav);
  
  return;
  
 error:
  abort();
}

void vect_RandNorm_dreal(vect_dreal *vect)
{
  int ndim, idim;
  
  ndim = vect->ndim;
  
  for(idim = 0; idim < ndim; ++idim)
    vect->addr[idim] = 2.0 * randnum_Get() - 1.0;
  
  vect_Renorm_dreal(vect);
}

void vect_RandNorm_dcmplx(vect_dcmplx *vect)
{
  int   ndim, idim;
  dreal re, im;
  
  ndim = vect->ndim;
  
  for(idim = 0; idim < ndim; ++idim) {
    re = 2.0 * randnum_Get() - 1.0;
    im = 2.0 * randnum_Get() - 1.0;
    vect->addr[idim] = re + I * im;
  }
  
  vect_Renorm_dcmplx(vect);
}
