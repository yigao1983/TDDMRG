#include <stdlib.h>
#include <math.h>
#include "dbg.h"
#include "sys.h"
#include "constants.h"
#include "fourier.h"
#include "numeric.h"

dreal numeric_Maxval_dreal(int ng, const dreal *fval)
{
  int   ig;
  dreal maxval;
  
  maxval = fval[0];
  for(ig = 0; ig < ng; ++ig) {
    if(fval[ig] > maxval) maxval = fval[ig];
  }
  
  return maxval;
}

dreal numeric_Minval_dreal(int ng, const dreal *fval)
{
  int   ig;
  dreal minval;
  
  minval = fval[0];
  for(ig = 0; ig < ng; ++ig) {
    if(fval[ig] < minval) minval = fval[ig];
  }
  
  return minval;
}

dreal numeric_Quadrat_dreal(int ng, const dreal *xval, const dreal *fval)
{
  int   ig;
  dreal dxvl, sumf, quad;
  
  quad = 0.0;
  
  if(ng > 1) {
    for(ig = 0; ig < ng - 1; ++ig) {
      dxvl = xval[ig+1] - xval[ig];
      sumf = fval[ig+1] + fval[ig];
      quad+= 0.5 * dxvl * sumf;
    }
  }
  
  return quad;
}

dcmplx numeric_Quadrat_dcmplx(int ng, const dreal *xval, const dcmplx *fval)
{
  int    ig;
  dreal  dxvl;
  dcmplx sumf, quad;
  
  quad = 0.0;
  
  if(ng > 1) {
    for(ig = 0; ig < ng - 1; ++ig) {
      dxvl = xval[ig+1] - xval[ig];
      sumf = fval[ig+1] + fval[ig];
      quad+= 0.5 * dxvl * sumf;
    }
  }
  
  return quad;
}

int numeric_Ipow(int base, int exp)
{
  int result = 1;
  
  while(exp) {
    if(exp & 1) {
      result *= base;
    }
    exp >>= 1;
    base *= base;
  }
  
  return result;
}

int numeric_Ispow2(int x)
{
  if(x == 0) {
    return 1;
  } else if(x == 1) {
    return 0;
  } else {
    while(((x % 2) == 0) && x > 1) /* While x is even and > 1 */
      x /= 2;
    return (x == 1);
  }
}

void numeric_Kronig(int ngrd, const dcmplx *gg, dcmplx *kk)
{
  int    ig;
  dcmplx *fr = NULL, *ff = NULL;
  
  check(ngrd % 2 == 0, "Invalid ngrd");
  
  fr = (dcmplx *) calloc(ngrd, sizeof(dcmplx));
  check_mem(fr, "fr");
  ff = (dcmplx *) calloc(ngrd, sizeof(dcmplx));
  check_mem(ff, "ff");
  
  for(ig = 0; ig < ngrd/2; ++ig) {
    fr[ig]        = gg[ig+ngrd/2];
    fr[ig+ngrd/2] = gg[ig];
  }
  
  fourier_dcmplx_1d(ngrd, fr, ff, fourier_forward);
  
  ff[0] *= 0.5; // Removing plateau
  for(ig = 0; ig < ngrd / 2; ++ig) {
    ff[ig+ngrd/2] = 0.0;
  }
  
  fourier_dcmplx_1d(ngrd, ff, fr, fourier_backward);
  
  for(ig = 0; ig < ngrd/2; ++ig) {
    kk[ig]        = fr[ig+ngrd/2] / ngrd;
    kk[ig+ngrd/2] = fr[ig] / ngrd;
  }
  
  freeup(ff); freeup(fr);
  
  return;
  
 error:
  abort();
}

void numeric_Convolut(int ngrd, dreal dx, int isigna, const dcmplx *fa,
                      int isignb, const dcmplx *fb, dcmplx *fc)
{
  int    igrd;
  dcmplx *fra = NULL, *frb = NULL, *frc = NULL,
         *ffa = NULL, *ffb = NULL, *ffc = NULL;
  
  check(ngrd % 2 == 0, "Invalid ngrd: %5d", ngrd);
  
  fra = (dcmplx *) calloc(ngrd, sizeof(dcmplx));
  frb = (dcmplx *) calloc(ngrd, sizeof(dcmplx));
  frc = (dcmplx *) calloc(ngrd, sizeof(dcmplx));
  ffa = (dcmplx *) calloc(ngrd, sizeof(dcmplx));
  ffb = (dcmplx *) calloc(ngrd, sizeof(dcmplx));
  ffc = (dcmplx *) calloc(ngrd, sizeof(dcmplx));
  
  for(igrd = 0; igrd < ngrd/2; ++igrd) {
    fra[igrd]        = fa[igrd+ngrd/2];
    fra[igrd+ngrd/2] = fa[igrd];
    frb[igrd]        = fb[igrd+ngrd/2];
    frb[igrd+ngrd/2] = fb[igrd];
  }
  
  fourier_dcmplx_1d(ngrd, fra, ffa, isigna * fourier_forward);
  fourier_dcmplx_1d(ngrd, frb, ffb, isignb * fourier_forward);
  
  for(igrd = 0; igrd < ngrd; ++igrd) {
    ffc[igrd] = ffa[igrd] * ffb[igrd];
  }
  
  fourier_dcmplx_1d(ngrd, ffc, frc, fourier_backward);
  
  for(igrd = 0; igrd < ngrd/2; ++igrd) {
    fc[igrd]        = frc[igrd+ngrd/2] * dx / ngrd;
    fc[igrd+ngrd/2] = frc[igrd]        * dx / ngrd;
  }
  
  free(ffc); free(ffb); free(ffa);
  free(frc); free(frb); free(fra);
  
  return;
  
 error:
  abort();
}
