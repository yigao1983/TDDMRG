#include <stdlib.h>
#include "dbg.h"
#include "vector.h"
#include "conjgrad.h"

static const int   Nstep   = 16;
static const dreal Tolnorm = 1e-9;
static int         Ndim;
static Operator    CGopr, CGopr2;

static void conjgrad_Opr2(const vect_dreal *vect, vect_dreal *O2vect)
{
  vect_dreal *Ovect = NULL;
  
  Ovect = vect_New_dreal();
  check_mem(Ovect, "Ovect");
  Ovect->Alloc(Ndim, Ovect);
  
  CGopr(vect,  Ovect);
  CGopr(Ovect, O2vect);
  
  vect_Del_dreal(Ovect);
  
  return;
  
 error:
  abort();
}

static void conjgrad_Iter(const vect_dreal *Yvect, vect_dreal *Xvect, Operator Opr)
{
  int        istep, cnvg, idim;
  dreal      rnorm;
  dreal      alpha, beta, numer, denom;
  vect_dreal *Rvect = NULL, *Pvect = NULL, *Opvect = NULL;
  
  Rvect = vect_New_dreal();
  check_mem(Rvect, "Rvect");
  Rvect->Alloc(Ndim, Rvect);
  
  Pvect = vect_New_dreal();
  check_mem(Pvect, "Pvect");
  Pvect->Alloc(Ndim, Pvect);
  
  Opvect = vect_New_dreal();
  check_mem(Opvect, "Opvect");
  Opvect->Alloc(Ndim, Opvect);
  
  Xvect->Reset(Xvect);
  
  vect_Copy_dreal(Yvect, Rvect);
  vect_Copy_dreal(Rvect, Pvect);
  
  for(cnvg = 0, istep = 0; istep < Nstep && !cnvg; ++istep) {
    Opr(Pvect, Opvect);
    numer = vect_Dotprod_dreal(Rvect, Rvect);
    denom = vect_Dotprod_dreal(Pvect, Opvect);
    alpha = numer / denom;
    
    for(idim = 0; idim < Ndim; ++idim) {
      Xvect->addr[idim] += alpha*Pvect->addr[idim];
      Rvect->addr[idim] -= alpha*Opvect->addr[idim];
    }
    
    rnorm = vect_Norm_dreal(Rvect);
    if(rnorm < Tolnorm) {
      cnvg = 1;
    } else {
      denom = numer;
      numer = vect_Dotprod_dreal(Rvect, Rvect);
      beta  = numer / denom;
      for(idim = 0; idim < Ndim; ++idim)
        Pvect->addr[idim] = Rvect->addr[idim] + beta * Pvect->addr[idim];
    }
  }
  
  vect_Del_dreal(Opvect);
  vect_Del_dreal(Pvect);
  vect_Del_dreal(Rvect);
  
  return;
  
 error:
  abort();
}

void conjgrad_Posdef(int Ngrd, const vect_dreal *Yvect, vect_dreal *Xvect, Operator Opr)
{
  Ndim = Ngrd;
  
  CGopr = Opr;
  
  conjgrad_Iter(Yvect, Xvect, CGopr);
}

void conjgrad_Symmtr(int Ngrd, const vect_dreal *Yvect, vect_dreal *Xvect, Operator Opr)
{
  vect_dreal *Oyvect = NULL;
  
  Ndim = Ngrd;
  
  CGopr = Opr;
  
  CGopr2 = conjgrad_Opr2;
  
  Oyvect = vect_New_dreal();
  check_mem(Oyvect, "Oyvect");
  Oyvect->Alloc(Ndim, Oyvect);
  
  Opr(Yvect, Oyvect);
  
  conjgrad_Iter(Oyvect, Xvect, CGopr2);
  
  vect_Del_dreal(Oyvect);
  
  return;
  
 error:
  abort();
}
/*
void conjgrad_Opr(const vect_dreal *vect, vect_dreal *Ovect)
{
  int ndim, idim;
  
  ndim = vect->ndim;
  
  Ovect->Reset(Ovect);
  
  for(idim = 0; idim < ndim; ++idim) {
    Ovect->addr[idim] +=(-ndim+idim) * vect->addr[idim];
    if(idim+1 < ndim) Ovect->addr[idim] += 0.5 * vect->addr[idim+1];
    if(idim-1 > -1)   Ovect->addr[idim] += 0.5 * vect->addr[idim-1];
  }
}

void conjgrad_Check(void)
{
  int        Ngrd = 512;
  vect_dreal *Yv = NULL, *Xv = NULL;
  
  Yv = vect_New_dreal();
  Yv->Alloc(Ngrd, Yv);
  Xv = vect_New_dreal();
  Xv->Alloc(Ngrd, Xv);
  
  vect_RandNorm_dreal(Yv); Yv->Print("Yv", Yv);
  conjgrad_Symmtr(Ngrd, Yv, Xv, conjgrad_Opr);
  
  Yv->Reset(Yv); conjgrad_Opr(Xv, Yv);
  Yv->Print("Yvnew", Yv);
  
  vect_Del_dreal(Xv);
  vect_Del_dreal(Yv);
}
*/
