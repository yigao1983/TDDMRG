#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "dbg.h"
#include "sys.h"
#include "randnum.h"
#include "matrix.h"
#include "vector.h"
#include "lapack.h"
#include "numeric.h"
#include "davidson.h"

static const int   Nsubmax = 150;
static const dreal Tolnorm = 1e-10;
static const char  *davlabel = "davidson:       ";

static int         Ndim;
static dreal       Rho;
static mat2d_dreal *Matsub = NULL;
static vect_dreal  *Xvect = NULL, *Rvect = NULL, *Cvect = NULL;
static vect_dreal  **Basis = NULL, **Obasis = NULL;
static Operator    Davopr;
static Precondt    Davpcd;

static void davidson_Init(void)
{
  int isub;
  
  Xvect = vect_New_dreal();
  check_mem(Xvect, "Xvect");
  Xvect->Alloc(Ndim, Xvect);
  
  Rvect = vect_New_dreal();
  check_mem(Rvect, "Rvect");
  Rvect->Alloc(Ndim, Rvect);
  
  Cvect = vect_New_dreal();
  check_mem(Cvect, "Cvect");
  Cvect->Alloc(Ndim, Cvect);
  
  Basis = (vect_dreal **) calloc(Nsubmax, sizeof(vect_dreal *));
  for(isub = 0; isub < Nsubmax; ++isub) {
    Basis[isub] = vect_New_dreal();
    check_mem(Basis[isub], "Basis[isub]");
    Basis[isub]->Alloc(Ndim, Basis[isub]);
  }
  
  Obasis = (vect_dreal **) calloc(Nsubmax, sizeof(vect_dreal *));
  for(isub = 0; isub < Nsubmax; ++isub) {
    Obasis[isub] = vect_New_dreal();
    check_mem(Obasis[isub], "Obasis[isub]");
    Obasis[isub]->Alloc(Ndim, Obasis[isub]);
  }
  
  Matsub = mat2d_New_dreal();
  check_mem(Matsub, "Matsub");
  Matsub->Alloc(Nsubmax, Nsubmax, Matsub);
  
  return;
  
 error:
  abort();
}

static void davidson_Quit(void)
{
  int isub;
  
  mat2d_Del_dreal(Matsub);
  
  for(isub = Nsubmax-1; isub > -1; --isub)
    vect_Del_dreal(Obasis[isub]);
  freeup(Obasis);
  
  for(isub = Nsubmax-1; isub > -1; --isub)
    vect_Del_dreal(Basis[isub]);
  freeup(Basis);
  
  vect_Del_dreal(Cvect);
  vect_Del_dreal(Rvect);
  vect_Del_dreal(Xvect);
}

static void davidson_Expand(int msub)
{
  int        isub;
  dreal      vOprv;
  vect_dreal *Cbasis = NULL, *Ocbasis = NULL;
  
  Cbasis = vect_New_dreal();
  check_mem(Cbasis, "Cbasis");
  Cbasis->Alloc(Ndim, Cbasis);
  
  Ocbasis = vect_New_dreal();
  check_mem(Ocbasis, "Ocbasis");
  Ocbasis->Alloc(Ndim, Ocbasis);
  
  vect_Copy_dreal(Cvect, Cbasis);
  
  if(msub > 1) vect_Orthog_dreal(msub-1, Basis, Cbasis);
  
  vect_Renorm_dreal(Cbasis);
  
  Davopr(Cbasis, Ocbasis);
  
  vect_Copy_dreal(Cbasis,  Basis[msub-1]);
  vect_Copy_dreal(Ocbasis, Obasis[msub-1]);
  
  for(isub = 0; isub < msub; ++isub) {
    vOprv = vect_Dotprod_dreal(Basis[isub], Ocbasis);
    Matsub->ptr[msub-1][isub] = vOprv;
    Matsub->ptr[isub][msub-1] = vOprv;
  }
  
  vect_Del_dreal(Ocbasis);
  vect_Del_dreal(Cbasis);
  
  return;
  
 error:
  abort();
}

static void davidson_Hmat(int msub, dreal *Hmat)
{
  int isub, jsub;
  
  for(isub = 0; isub < msub; ++isub) {
    for(jsub = 0; jsub < msub; ++jsub) {
      Hmat[isub*msub+jsub] = Matsub->ptr[isub][jsub];
    }
  }
}

static void davidson_Xvect(int msub, const dreal *Hvect)
{
  int idim, isub;
  
  Xvect->Reset(Xvect);
  
  for(idim = 0; idim < Ndim; ++idim) {
    for(isub = 0; isub < msub; ++isub) {
      Xvect->addr[idim] += Hvect[isub] * Basis[isub]->addr[idim]; 
    }
  }
}

static void davidson_Rvect(int msub, const dreal *Hvect)
{
  int        idim, isub;
  vect_dreal *Oxvect = NULL;
  
  Rvect->Reset(Rvect);
  
  Oxvect = vect_New_dreal();
  check_mem(Oxvect, "Oxvect");
  Oxvect->Alloc(Ndim, Oxvect);
  
  for(idim = 0; idim < Ndim; ++idim) {
    for(isub = 0; isub < msub; ++isub) {
      Oxvect->addr[idim] += Hvect[isub] * Obasis[isub]->addr[idim];
    }
    Rvect->addr[idim] = Oxvect->addr[idim] - Rho * Xvect->addr[idim];
  }
  
  vect_Del_dreal(Oxvect);
  
  return;
  
 error:
  abort();
}

static void davidson_Cvect_Init(vect_dreal *Evc)
{
  dreal cnorm;
  
  vect_Copy_dreal(Evc, Cvect);
  
  cnorm = vect_Norm_dreal(Cvect);
  
  if(cnorm < Tolnorm) vect_RandNorm_dreal(Cvect);
}

static void davidson_Cvect(void)
{
  Davpcd(Rho, Rvect, Cvect);
}

void davidson_Diag(int Ngrd, dreal *Eig, vect_dreal *Evc, Operator Opr, Precondt Pcd)
{
  int   mmax, msub, cnvg;
  dreal rnorm;
  dreal *Earr = NULL;
  dreal *Hmat = NULL;
  
  Ndim = Ngrd;
  
  Davopr = Opr;
  Davpcd = Pcd;
  
  davidson_Init();
  
  randnum_Set();
  
  davidson_Cvect_Init(Evc);
  
  mmax = (Ndim > Nsubmax) ? Nsubmax : Ndim;
  
  if(Rank == Root) {
    printf("%s\n", davlabel);
    printf("%s ==== Davidson ====\n", davlabel);
    printf("%s  msub       rnorm \n", davlabel);
    printf("%s ------------------\n", davlabel);
  }
  for(msub = 1, cnvg = 0; msub <= mmax && !cnvg; ++msub) {
    davidson_Expand(msub);
    Earr = (dreal *) calloc(msub, sizeof(dreal));
    Hmat = (dreal *) calloc(msub*msub, sizeof(dreal));
    davidson_Hmat(msub, Hmat);
    lapack_dsyev(msub, Hmat, Earr);
    Rho = Earr[0];
    davidson_Xvect(msub, Hmat);
    davidson_Rvect(msub, Hmat);
    rnorm = vect_Norm_dreal(Rvect);
    if(Rank == Root && numeric_Ispow2(msub)) {
      printf("%s %5d %10.5e\n", davlabel, msub, rnorm);
      fflush(stdout);
    }
    if(rnorm < Tolnorm) {
      cnvg = 1;
    } else {
      davidson_Cvect();
    }
    freeup(Hmat); freeup(Earr);
  }
  if(Rank == Root) {
    printf("%s ==================\n", davlabel);
    if(!cnvg) printf("%s %s\n", davlabel, "Convergence not reached!");
    printf("%s\n", davlabel);
    fflush(stdout);
  }
  
  *Eig = Rho; vect_Copy_dreal(Xvect, Evc);
  
  davidson_Quit();
}
