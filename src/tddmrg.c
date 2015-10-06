#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include "dbg.h"
#include "sys.h"
#include "constants.h"
#include "vector.h"
#include "matrix.h"
#include "lapack.h"
#include "numeric.h"
#include "parameters.h"
#include "dmrg.h"

static int          Nstep;
static dreal        dT, Alpha, Time, Coef;
static dcmplx       JL, JR;
static dreal        *Bchain = NULL, *Fchain = NULL;
static vect_dcmplx  *TDvc = NULL;
static mat2d_dreal  *HblockL = NULL, *HblockR = NULL, *HsiteL = NULL, *HsiteR = NULL;
static mat2d_dreal  *VblocksiteL = NULL, *VblocksiteR = NULL, *Vsitesite = NULL;
static mat2d_dcmplx *Dmat_L2R = NULL, *Dmat_R2L = NULL;
static mat2d_dreal  **SZblockL = NULL, **SZblockR = NULL, **SPblockL = NULL, **SPblockR = NULL;
static mat2d_dreal  **SZsiteL = NULL, **SZsiteR = NULL, **SPsiteL = NULL, **SPsiteR = NULL;
static mat2d_dcmplx **zWblock_L2R = NULL, **zWblock_R2L = NULL;

static const char *tdchainin   = "tdchain.in";
static const char *tddmrglabel = "tddmrg:         ";
static const char *tdcurrent   = "tdCurrent.ascii";

static void tddmrg_dvect2zvect(const vect_dreal *dvect, vect_dcmplx *zvect)
{
  int ndim, idim;
  
  ndim = dvect->ndim;
  
  for(idim = 0; idim < ndim; ++idim)
    zvect->addr[idim] = dvect->addr[idim];
}

static void tddmrg_dmat2d2zmat2d(const mat2d_dreal *dmat, mat2d_dcmplx *zmat)
{
  int nrow, ncol, irow, icol;
  
  nrow = dmat->nrow;
  ncol = dmat->ncol;
  
  for(icol = 0; icol < ncol; ++icol)
    for(irow = 0; irow < nrow; ++irow)
      zmat->ptr[icol][irow] = dmat->ptr[icol][irow];
}

static void tddmrg_Jmat(dreal tt, int ndim1, int ndim2,
                        mat2d_dreal **SPmat1, mat2d_dreal **SPmat2,
                        mat2d_dcmplx *Jmat)
{
  int   iflav, idim1, jdim1, idim2, jdim2, idim, jdim;
  dreal xp1, xm1, xp2, xm2;
  
  Jmat->Reset(Jmat);
  
  for(iflav = 0; iflav < Nflavor; ++iflav) {
    for(idim2 = 0; idim2 < ndim2; ++idim2) {
      for(idim1 = 0; idim1 < ndim1; ++idim1) {
        idim = idim2*ndim1 + idim1;
        for(jdim2 = 0; jdim2 < ndim2; ++jdim2) {
          for(jdim1 = 0; jdim1 < ndim1; ++jdim1) {
            jdim = jdim2*ndim1 + jdim1;
            xp1 = SPmat1[iflav]->ptr[idim1][jdim1];
            xm1 = SPmat1[iflav]->ptr[jdim1][idim1];
            xp2 = SPmat2[iflav]->ptr[idim2][jdim2];
            xm2 = SPmat2[iflav]->ptr[jdim2][idim2];
            Jmat->ptr[idim][jdim] += I * tt*(xp1*xm2 - xm1*xp2);
          }
        }
      }
    }
  }
}

static void tddmrg_Input(void)
{
  FILE *fp = NULL;
  int  iblock, ndimBB;
  int  *ndimW = NULL, *mdimW = NULL;
  
  if(Rank == Root) {
    printf("%s Reading input for tddmrg...\n", tddmrglabel);
    
    printf("%s Reading wavefunction...\n", tddmrglabel);
    fp = fopen(evcgsbin, "r");
    fread(&ndimBB, sizeof(int), 1, fp);
    EvcGS->Deall(EvcGS);
    EvcGS->Alloc(ndimBB, EvcGS);
    fread(EvcGS->addr, sizeof(dreal), ndimBB, fp);
    fclose(fp);
    
    printf("%s ndimBB = %5d\n%s\n", tddmrglabel, ndimBB, tddmrglabel);
    
    printf("%s Reading basis...\n", tddmrglabel);
    fp = fopen(wblockbin, "r");
    ndimW = (int *) calloc(Nblock, sizeof(int));
    mdimW = (int *) calloc(Nblock, sizeof(int));
    for(iblock = 0; iblock < Nblock; ++iblock) {
      fread(ndimW+iblock, sizeof(int), 1, fp);
      fread(mdimW+iblock, sizeof(int), 1, fp);
      Wblock[iblock]->Deall(Wblock[iblock]);
      Wblock[iblock]->Alloc(ndimW[iblock], mdimW[iblock], Wblock[iblock]);
      fread(Wblock[iblock]->addr, sizeof(dreal), ndimW[iblock]*mdimW[iblock], fp);
    }
    fclose(fp);
    
    printf("%s iblock, ndimW, mdimW =\n", tddmrglabel);
    for(iblock = 0; iblock < Nblock; ++iblock) {
      printf("%s %5d %5d %5d\n", tddmrglabel, iblock, ndimW[iblock], mdimW[iblock]);
    }
    printf("%s\n", tddmrglabel);
    
    printf("%s Reading successfully completed\n%s\n", tddmrglabel, tddmrglabel);
  }
  
  sys_Bcast_int(&ndimBB, 1, Root);
  sys_Bcast_dreal(EvcGS->addr, ndimBB, Root); 
  
  TDvc->Deall(TDvc);
  TDvc->Alloc(ndimBB, TDvc);
  tddmrg_dvect2zvect(EvcGS, TDvc);
  
  if(Rank != Root) {
    ndimW = (int *) calloc(Nblock, sizeof(int));
    mdimW = (int *) calloc(Nblock, sizeof(int));
  }
  
  sys_Bcast_int(ndimW, Nblock, Root);
  sys_Bcast_int(mdimW, Nblock, Root);
  
  if(Rank != Root) {
    EvcGS->Deall(EvcGS);
    EvcGS->Alloc(ndimBB, EvcGS);
    for(iblock = 0; iblock < Nblock; ++iblock) {
      Wblock[iblock]->Deall(Wblock[iblock]);
      Wblock[iblock]->Alloc(ndimW[iblock], mdimW[iblock], Wblock[iblock]);
    }
  }
  
  for(iblock = 0; iblock < Nblock; ++iblock)
    sys_Bcast_dreal(Wblock[iblock]->addr, ndimW[iblock]*mdimW[iblock], Root);
  
  for(iblock = 0; iblock < Nblock; ++iblock) {
    zWblock_L2R[iblock]->Deall(zWblock_L2R[iblock]);
    zWblock_L2R[iblock]->Alloc(ndimW[iblock], mdimW[iblock], zWblock_L2R[iblock]);
    tddmrg_dmat2d2zmat2d(Wblock[iblock], zWblock_L2R[iblock]);
  }
  
  for(iblock = 0; iblock < Nblock; ++iblock) {
    zWblock_R2L[iblock]->Deall(zWblock_R2L[iblock]);
    zWblock_R2L[iblock]->Alloc(ndimW[iblock], mdimW[iblock], zWblock_R2L[iblock]);
    tddmrg_dmat2d2zmat2d(Wblock[iblock], zWblock_R2L[iblock]);
  }
  
  freeup(mdimW);
  freeup(ndimW);
}

static void tddmrg_Print(void)
{
  int len;
  
  if(Rank == Root) {
    printf("%s\n", tddmrglabel);
    printf("%s TDDMRG parameters:\n", tddmrglabel);
    printf("%s\n", tddmrglabel);
    printf("%s dT     = %15.5e\n", tddmrglabel, dT);
    printf("%s Nstep  = %5d\n",    tddmrglabel, Nstep);
    printf("%s Alpha  = %15.5e\n", tddmrglabel, Alpha);
    printf("%s B (Ha) =\n", tddmrglabel);
    for(len = 0; len < Length; ++len)
      printf("%s %5d %15.5e\n", tddmrglabel, len, Bchain[len]);
    printf("%s\n", tddmrglabel);
  }
}

static void tddmrg_Read(void)
{
  FILE *fp = NULL;
  int  len;
  char line[SHRT_MAX];
  
  if(Rank == Root) {
    fp = fopen(tdchainin, "r");
    if(fgets(line, sizeof(line), fp)) sscanf(line, "%lf %d", &dT, &Nstep);
    if(fgets(line, sizeof(line), fp)) sscanf(line, "%lf" , &Alpha);
    for(len = 0; len < Length; ++len)
      if(fgets(line, sizeof(line), fp)) sscanf(line, "%lf" , Bchain+len);
    fclose(fp);
  }
  
  sys_Bcast_dreal(Bchain, Length, Root);
  
  for(len = 0; len < Length; ++len)
    Bchain[len] /= EV_HA;
  
  tddmrg_Print();
}

static void tddmrg_Echain(void)
{
  int len;
  
  for(len = 0; len < Length; ++len) {
    Echain[len] = Fchain[len] + Coef * Bchain[len];
  }
}

static void tddmrg_TDvc_Save(const vect_dcmplx *TDvc_src, vect_dcmplx *TDvc_des)
{
  int ndimBB;
  
  ndimBB = TDvc_src->ndim;
  
  TDvc_des->Deall(TDvc_des);
  TDvc_des->Alloc(ndimBB, TDvc_des);
  
  vect_Copy_dcmplx(TDvc_src, TDvc_des);
}

void tddmrg_Init(void)
{
  int iflav, iblock;
  
  Bchain = (dreal *) calloc(Length, sizeof(dreal));
  Fchain = (dreal *) calloc(Length, sizeof(dreal));
  
  memcpy(Fchain, Echain, Length*sizeof(dreal));
  
  tddmrg_Read();
  
  HblockL = mat2d_New_dreal(); HblockR = mat2d_New_dreal();
  HsiteL  = mat2d_New_dreal(); HsiteR  = mat2d_New_dreal();
  
  SZblockL = (mat2d_dreal **) calloc(Nflavor, sizeof(mat2d_dreal *));
  SZblockR = (mat2d_dreal **) calloc(Nflavor, sizeof(mat2d_dreal *));
  SPblockL = (mat2d_dreal **) calloc(Nflavor, sizeof(mat2d_dreal *));
  SPblockR = (mat2d_dreal **) calloc(Nflavor, sizeof(mat2d_dreal *));
  
  SZsiteL = (mat2d_dreal **) calloc(Nflavor, sizeof(mat2d_dreal *));
  SZsiteR = (mat2d_dreal **) calloc(Nflavor, sizeof(mat2d_dreal *));
  SPsiteL = (mat2d_dreal **) calloc(Nflavor, sizeof(mat2d_dreal *));
  SPsiteR = (mat2d_dreal **) calloc(Nflavor, sizeof(mat2d_dreal *));
  for(iflav = 0; iflav < Nflavor; ++iflav) {
    SZblockL[iflav] = mat2d_New_dreal(); SZblockR[iflav] = mat2d_New_dreal();
    SPblockL[iflav] = mat2d_New_dreal(); SPblockR[iflav] = mat2d_New_dreal();
    
    SZsiteL[iflav] = mat2d_New_dreal(); SZsiteR[iflav] = mat2d_New_dreal();
    SPsiteL[iflav] = mat2d_New_dreal(); SPsiteR[iflav] = mat2d_New_dreal();
  }
  
  VblocksiteL = mat2d_New_dreal();
  VblocksiteR = mat2d_New_dreal();
  Vsitesite   = mat2d_New_dreal();
  
  TDvc = vect_New_dcmplx();
  
  zWblock_L2R = (mat2d_dcmplx **) calloc(Nblock, sizeof(mat2d_dcmplx *));
  zWblock_R2L = (mat2d_dcmplx **) calloc(Nblock, sizeof(mat2d_dcmplx *));
  
  for(iblock = 0; iblock < Nblock; ++iblock) {
    zWblock_L2R[iblock] = mat2d_New_dcmplx();
    zWblock_R2L[iblock] = mat2d_New_dcmplx();
  }
  
  Dmat_L2R = mat2d_New_dcmplx();
  Dmat_R2L = mat2d_New_dcmplx();
  
  tddmrg_Input();
}

void tddmrg_Quit(void)
{
  int iblock, iflav;
  
  mat2d_Del_dcmplx(Dmat_R2L);
  mat2d_Del_dcmplx(Dmat_L2R);
  
  for(iblock = Nblock-1; iblock > -1; --iblock) {
    mat2d_Del_dcmplx(zWblock_R2L[iblock]);
    mat2d_Del_dcmplx(zWblock_L2R[iblock]);
  }
  
  freeup(zWblock_R2L);
  freeup(zWblock_L2R);
  
  vect_Del_dcmplx(TDvc);
  
  mat2d_Del_dreal(Vsitesite);
  mat2d_Del_dreal(VblocksiteR);
  mat2d_Del_dreal(VblocksiteL);
  
  if(SPsiteR) {
    for(iflav = Nflavor-1; iflav > -1; --iflav)
      mat2d_Del_dreal(SPsiteR[iflav]);
  } freeup(SPsiteR);
  if(SPsiteL) {
    for(iflav = Nflavor-1; iflav > -1; --iflav)
      mat2d_Del_dreal(SPsiteL[iflav]);
  } freeup(SPsiteL);
  
  if(SZsiteR) {
    for(iflav = Nflavor-1; iflav > -1; --iflav)
      mat2d_Del_dreal(SZsiteR[iflav]);
  } freeup(SZsiteR);
  if(SZsiteL) {
    for(iflav = Nflavor-1; iflav > -1; --iflav)
      mat2d_Del_dreal(SZsiteL[iflav]);
  } freeup(SZsiteL);
  
  if(SPblockR) {
    for(iflav = Nflavor-1; iflav > -1; --iflav)
      mat2d_Del_dreal(SPblockR[iflav]);
  } freeup(SPblockR);
  if(SPblockL) {
    for(iflav = Nflavor-1; iflav > -1; --iflav)
      mat2d_Del_dreal(SPblockL[iflav]);
  } freeup(SPblockL);
  
  if(SZblockR) {
    for(iflav = Nflavor-1; iflav > -1; --iflav)
      mat2d_Del_dreal(SZblockR[iflav]);
  } freeup(SZblockR);
  if(SZblockL) {
    for(iflav = Nflavor-1; iflav > -1; --iflav)
      mat2d_Del_dreal(SZblockL[iflav]);
  } freeup(SZblockL);
  
  mat2d_Del_dreal(HsiteR);  mat2d_Del_dreal(HsiteL);
  mat2d_Del_dreal(HblockR); mat2d_Del_dreal(HblockL);
  
  freeup(Fchain);
  freeup(Bchain);
}

static void tddmrg_expMat(dreal tau, const mat2d_dreal *Mat, mat2d_dcmplx *expMat)
{
  int    ndim, idim, jdim;
  dreal  *mat = NULL, *eig = NULL;
  dcmplx *zmat = NULL, *expdiag = NULL, *tmpmat = NULL;
  
  ndim = Mat->nrow;
  
  expMat->Reset(expMat);
  
  mat = (dreal *) calloc(ndim*ndim, sizeof(dreal));
  eig = (dreal *) calloc(ndim, sizeof(dreal));
  
  zmat    = (dcmplx *) calloc(ndim*ndim, sizeof(dcmplx));
  expdiag = (dcmplx *) calloc(ndim*ndim, sizeof(dcmplx));
  tmpmat  = (dcmplx *) calloc(ndim*ndim, sizeof(dcmplx));
  
  memcpy(mat, Mat->addr, ndim*ndim*sizeof(dreal));
  
  lapack_dsyev(ndim, mat, eig);
  
  for(idim = 0; idim < ndim; ++idim)
    for(jdim = 0; jdim < ndim; ++jdim)
      zmat[idim*ndim+jdim] = (dcmplx) mat[idim*ndim+jdim];
  
  for(idim = 0; idim < ndim; ++idim)
    expdiag[idim*ndim+idim] = cexp(-I * eig[idim] * tau);
  
  lapack_zgemm(ndim, ndim, ndim, 'N', 'C', 1.0, expdiag, zmat, 0.0, tmpmat); 
  lapack_zgemm(ndim, ndim, ndim, 'N', 'N', 1.0, zmat, tmpmat, 0.0, expMat->addr); 
  
  freeup(tmpmat);
  freeup(expdiag);
  freeup(zmat);
  freeup(eig);
  freeup(mat);
}
// Block operator block-site vector multiplication, NO value added
static void tddmrg_OprBvectBS(int ndimB, int ndimS, int mdimB, const mat2d_dcmplx *OprB,
                              const vect_dcmplx *vectBS, vect_dcmplx *OprBvectBS)
{
  int         idimS;
  vect_dcmplx *vectB = NULL, *OprBvectB = NULL;
  
  vectB = vect_New_dcmplx();
  check_mem(vectB, "vectB");
  vectB->Alloc(mdimB, vectB);
  
  OprBvectB = vect_New_dcmplx();
  check_mem(OprBvectB, "OprBvectB");
  OprBvectB->Alloc(ndimB, OprBvectB);
  
  for(idimS = 0; idimS < ndimS; ++idimS) {
    memcpy(vectB->addr, vectBS->addr+idimS*mdimB, sizeof(dcmplx)*mdimB);
    lapack_zgemm(ndimB, 1, mdimB, 'N', 'N', 1.0, OprB->addr, vectB->addr, 0.0, OprBvectB->addr);
    memcpy(OprBvectBS->addr+idimS*ndimB, OprBvectB->addr, sizeof(dcmplx)*ndimB);
  }
  
  vect_Del_dcmplx(OprBvectB);
  vect_Del_dcmplx(vectB); 
  
  return;
  
 error:
  abort();
}
// Site operator block-site vector multiplication, NO value added
static void tddmrg_OprSvectBS(int ndimB, int ndimS, int mdimS, const mat2d_dcmplx *OprS,
                              const vect_dcmplx *vectBS, vect_dcmplx *OprSvectBS)
{
  int         idimS, idimB, idimBS;
  vect_dcmplx *vectS = NULL, *OprSvectS = NULL;
  
  vectS = vect_New_dcmplx();
  check_mem(vectS, "vectS");
  vectS->Alloc(mdimS, vectS);
  
  OprSvectS = vect_New_dcmplx();
  check_mem(OprSvectS, "OprSvectS");
  OprSvectS->Alloc(ndimS, OprSvectS);
  
  for(idimB = 0; idimB < ndimB; ++idimB) {
    for(idimS = 0; idimS < mdimS; ++idimS) {
      idimBS = idimS*ndimB + idimB;
      vectS->addr[idimS] = vectBS->addr[idimBS];
    }
    lapack_zgemm(ndimS, 1, mdimS, 'N', 'N', 1.0, OprS->addr, vectS->addr, 0.0, OprSvectS->addr);
    for(idimS = 0; idimS < ndimS; ++idimS) {
      idimBS = idimS*ndimB + idimB;
      OprSvectBS->addr[idimBS] = OprSvectS->addr[idimS];
    }
  }
  
  vect_Del_dcmplx(OprSvectS);
  vect_Del_dcmplx(vectS); 
  
  return;
  
 error:
  abort();
}
// Block-site operator block-site vector multiplication, NO value added
static void tddmrg_OprBSvectBS(int ndimBS, int mdimBS, const mat2d_dcmplx *OprBS,
                               const vect_dcmplx *vectBS, vect_dcmplx *OprBSvectBS)
{
  lapack_zgemm(ndimBS, 1, mdimBS, 'N', 'N', 1.0, OprBS->addr, vectBS->addr, 0.0, OprBSvectBS->addr);
}
// BlockL operator superblock vector multiplication, NO value added
static void tddmrg_OprBLvect(int ndimBL, int ndimSL, int ndimBR, int ndimSR, int mdimBL,
                             const mat2d_dcmplx *OprBL, const vect_dcmplx *vect, vect_dcmplx *OprBLvect)
{
  int         ndimL, mdimL, ndimBB, idimBR, idimSR, idimR;
  vect_dcmplx *vectL = NULL, *OprBLvectL = NULL;
  
  ndimBB = ndimBL*ndimSL*ndimBR*ndimSL;
  ndimL  = ndimBL*ndimSL;
  mdimL  = mdimBL*ndimSL;
  
  vectL = vect_New_dcmplx();
  check_mem(vectL, "vectL");
  vectL->Alloc(mdimL, vectL);
  
  OprBLvectL = vect_New_dcmplx();
  check_mem(OprBLvectL, "OprBLvectL");
  OprBLvectL->Alloc(ndimL, OprBLvectL);
  
  for(idimSR = 0; idimSR < ndimSR; ++idimSR) {
    for(idimBR = 0; idimBR < ndimBR; ++idimBR) {
      idimR = idimSR*ndimBR + idimBR;
      memcpy(vectL->addr, vect->addr+idimR*mdimL, sizeof(dcmplx)*mdimL);
      tddmrg_OprBvectBS(ndimBL, ndimSL, mdimBL, OprBL, vectL, OprBLvectL);
      memcpy(OprBLvect->addr+idimR*ndimL, OprBLvectL->addr, sizeof(dcmplx)*ndimL);
    }
  }
  
  vect_Del_dcmplx(OprBLvectL);
  vect_Del_dcmplx(vectL);
  
  return;
  
 error:
  abort();
}
// SiteL operator superblock vector multiplication, NO value added
static void tddmrg_OprSLvect(int ndimBL, int ndimSL, int ndimBR, int ndimSR, int mdimSL,
                             const mat2d_dcmplx *OprSL, const vect_dcmplx *vect, vect_dcmplx *OprSLvect)
{
  int         ndimL, mdimL, ndimBB, idimBR, idimSR, idimR;
  vect_dcmplx *vectL = NULL, *OprSLvectL = NULL;
  
  ndimBB = ndimBL*ndimSL*ndimBR*ndimSL;
  ndimL  = ndimBL*ndimSL;
  mdimL  = ndimBL*mdimSL;
  
  vectL = vect_New_dcmplx();
  check_mem(vectL, "vectL");
  vectL->Alloc(mdimL, vectL);
  
  OprSLvectL = vect_New_dcmplx();
  check_mem(OprSLvectL, "OprSLvectL");
  OprSLvectL->Alloc(ndimL, OprSLvectL);
  
  for(idimSR = 0; idimSR < ndimSR; ++idimSR) {
    for(idimBR = 0; idimBR < ndimBR; ++idimBR) {
      idimR = idimSR*ndimBR + idimBR;
      memcpy(vectL->addr, vect->addr+idimR*mdimL, sizeof(dcmplx)*mdimL);
      tddmrg_OprSvectBS(ndimBL, ndimSL, mdimSL, OprSL, vectL, OprSLvectL);
      memcpy(OprSLvect->addr+idimR*ndimL, OprSLvectL->addr, sizeof(dcmplx)*ndimL);
    }
  }
  
  vect_Del_dcmplx(OprSLvectL);
  vect_Del_dcmplx(vectL);
  
  return;
  
 error:
  abort();
}
// BlockR operator superblock operator multiplication, NO value added
static void tddmrg_OprBRvect(int ndimBL, int ndimSL, int ndimBR, int ndimSR, int mdimBR,
                             const mat2d_dcmplx *OprBR, const vect_dcmplx *vect, vect_dcmplx *OprBRvect)
{
  int         ndimR, mdimR, ndimBB, idimBL, idimSL, idimBR, idimSR, idimR, idimBB;
  vect_dcmplx *vectR = NULL, *OprBRvectR = NULL;
  
  ndimBB = ndimBL*ndimSL*ndimBR*ndimSL;
  ndimR  = ndimBR*ndimSR;
  mdimR  = mdimBR*ndimSR;
  
  vectR = vect_New_dcmplx();
  check_mem(vectR, "vectR");
  vectR->Alloc(mdimR, vectR);
  
  OprBRvectR = vect_New_dcmplx();
  check_mem(OprBRvectR, "OprBRvectR");
  OprBRvectR->Alloc(ndimR, OprBRvectR);
  
  for(idimSL = 0; idimSL < ndimSL; ++idimSL) {
    for(idimBL = 0; idimBL < ndimBL; ++idimBL) {
      for(idimSR = 0; idimSR < ndimSR; ++idimSR) {
        for(idimBR = 0; idimBR < mdimBR; ++idimBR) {
          idimR  = idimSR*mdimBR + idimBR;
          idimBB = dmrg_indxBB(ndimBL, ndimSL, mdimBR, ndimSR, idimBL, idimSL, idimBR, idimSR);
          vectR->addr[idimR] = vect->addr[idimBB];
        }
      }
      tddmrg_OprBvectBS(ndimBR, ndimSR, mdimBR, OprBR, vectR, OprBRvectR);
      for(idimSR = 0; idimSR < ndimSR; ++idimSR) {
        for(idimBR = 0; idimBR < ndimBR; ++idimBR) {
          idimR  = idimSR*ndimBR + idimBR;
          idimBB = dmrg_indxBB(ndimBL, ndimSL, ndimBR, ndimSR, idimBL, idimSL, idimBR, idimSR);
          OprBRvect->addr[idimBB] = OprBRvectR->addr[idimR];
        }
      }
    }
  }
  
  vect_Del_dcmplx(OprBRvectR);
  vect_Del_dcmplx(vectR);
  
  return;
  
 error:
  abort();
}
// SiteR operator superblock vector multiplication, NO value added
static void tddmrg_OprSRvect(int ndimBL, int ndimSL, int ndimBR, int ndimSR, int mdimSR,
                             const mat2d_dcmplx *OprSR, const vect_dcmplx *vect, vect_dcmplx *OprSRvect)
{
  int         ndimR, mdimR, ndimBB, idimBL, idimSL, idimBR, idimSR, idimR, idimBB;
  vect_dcmplx *vectR = NULL, *OprSRvectR = NULL;
  
  ndimBB = ndimBL*ndimSL*ndimBR*ndimSL;
  ndimR  = ndimBR*ndimSR;
  mdimR  = ndimBR*mdimSR;
  
  vectR = vect_New_dcmplx();
  check_mem(vectR, "vectR");
  vectR->Alloc(mdimR, vectR);
  
  OprSRvectR = vect_New_dcmplx();
  check_mem(OprSRvectR, "OprSRvectR");
  OprSRvectR->Alloc(ndimR, OprSRvectR);
  
  for(idimSL = 0; idimSL < ndimSL; ++idimSL) {
    for(idimBL = 0; idimBL < ndimBL; ++idimBL) {
      for(idimSR = 0; idimSR < mdimSR; ++idimSR) {
        for(idimBR = 0; idimBR < ndimBR; ++idimBR) {
          idimR  = idimSR*ndimBR + idimBR;
          idimBB = dmrg_indxBB(ndimBL, ndimSL, ndimBR, mdimSR, idimBL, idimSL, idimBR, idimSR);
          vectR->addr[idimR] = vect->addr[idimBB];
        }
      }
      tddmrg_OprSvectBS(ndimBR, ndimSR, mdimSR, OprSR, vectR, OprSRvectR);
      for(idimSR = 0; idimSR < ndimSR; ++idimSR) {
        for(idimBR = 0; idimBR < ndimBR; ++idimBR) {
          idimR  = idimSR*ndimBR + idimBR;
          idimBB = dmrg_indxBB(ndimBL, ndimSL, ndimBR, ndimSR, idimBL, idimSL, idimBR, idimSR);
          OprSRvect->addr[idimBB] = OprSRvectR->addr[idimR];
        }
      }
    }
  }
  
  vect_Del_dcmplx(OprSRvectR);
  vect_Del_dcmplx(vectR);
  
  return;
  
 error:
  abort();
}
// L operator superblock vector multiplication, NO value added
static void tddmrg_OprLvect(int ndimBL, int ndimSL, int ndimBR, int ndimSR, int mdimBL, int mdimSL,
                            const mat2d_dcmplx *OprL, const vect_dcmplx *vect, vect_dcmplx *OprLvect)
{
  int        ndimBB, ndimL, mdimL, idimBR, idimSR, idimR;
  vect_dcmplx *vectL = NULL, *OprLvectL = NULL;
  
  ndimBB = ndimBL*ndimSL*ndimBR*ndimSR;
  ndimL  = ndimBL*ndimSL;
  mdimL  = mdimBL*mdimSL;
  
  vectL = vect_New_dcmplx();
  check_mem(vectL, "vectL");
  vectL->Alloc(mdimL, vectL);
  
  OprLvectL = vect_New_dcmplx();
  check_mem(OprLvectL, "OprLvectL");
  OprLvectL->Alloc(ndimL, OprLvectL);
  
  for(idimSR = 0; idimSR < ndimSR; ++idimSR) {
    for(idimBR = 0; idimBR < ndimBR; ++idimBR) {
      idimR = idimSR*ndimBR + idimBR;
      memcpy(vectL->addr, vect->addr+idimR*mdimL, sizeof(dcmplx)*mdimL);
      tddmrg_OprBSvectBS(ndimL, mdimL, OprL, vectL, OprLvectL);
      memcpy(OprLvect->addr+idimR*ndimL, OprLvectL->addr, sizeof(dcmplx)*ndimL);
    }
  }
  
  vect_Del_dcmplx(OprLvectL);
  vect_Del_dcmplx(vectL);
  
  return;
  
 error:
  abort();
}
// R operator superblock vector multiplication, NO value added
static void tddmrg_OprRvect(int ndimBL, int ndimSL, int ndimBR, int ndimSR, int mdimBR, int mdimSR,
                            const mat2d_dcmplx *OprR, const vect_dcmplx *vect, vect_dcmplx *OprRvect)
{
  int        ndimR, mdimR, ndimBB, idimBL, idimSL, idimBR, idimSR, idimR, idimBB;
  vect_dcmplx *vectR = NULL, *OprRvectR = NULL;
  
  ndimBB = ndimBL*ndimSL*ndimBR*ndimSL;
  ndimR  = ndimBR*ndimSR;
  mdimR  = mdimBR*mdimSR;
  
  vectR = vect_New_dcmplx();
  check_mem(vectR, "vectR");
  vectR->Alloc(mdimR, vectR);
  
  OprRvectR = vect_New_dcmplx();
  check_mem(OprRvectR, "OprRvectR");
  OprRvectR->Alloc(ndimR, OprRvectR);
  
  for(idimSL = 0; idimSL < ndimSL; ++idimSL) {
    for(idimBL = 0; idimBL < ndimBL; ++idimBL) {
      for(idimSR = 0; idimSR < mdimSR; ++idimSR) {
        for(idimBR = 0; idimBR < mdimBR; ++idimBR) {
          idimR  = idimSR*mdimBR + idimBR;
          idimBB = dmrg_indxBB(ndimBL, ndimSL, mdimBR, mdimSR, idimBL, idimSL, idimBR, idimSR);
          vectR->addr[idimR] = vect->addr[idimBB];
        }
      }
      tddmrg_OprBSvectBS(ndimR, mdimR, OprR, vectR, OprRvectR);
      for(idimSR = 0; idimSR < ndimSR; ++idimSR) {
        for(idimBR = 0; idimBR < ndimBR; ++idimBR) {
          idimR  = idimSR*ndimBR + idimBR;
          idimBB = dmrg_indxBB(ndimBL, ndimSL, ndimBR, ndimSR, idimBL, idimSL, idimBR, idimSR);
          OprRvect->addr[idimBB] = OprRvectR->addr[idimR];
        }
      }
    }
  }
  
  vect_Del_dcmplx(OprRvectR);
  vect_Del_dcmplx(vectR);
  
  return;
  
 error:
  abort();
}
// Site-site operator superblock vector multiplication, NO value added
static void tddmrg_OprSSvect(int ndimBL, int ndimSL, int ndimBR, int ndimSR, int mdimSL, int mdimSR,
                             const mat2d_dcmplx *OprSS, const vect_dcmplx *vect, vect_dcmplx *OprSSvect)
{
  int         ndimSS, mdimSS, idimBL, idimBR, idimSL, idimSR, idimSS, idimBB;
  vect_dcmplx *vectSS = NULL, *OprSSvectSS = NULL;
  
  ndimSS = ndimSL*ndimSR;
  mdimSS = mdimSL*mdimSR;
  
  vectSS = vect_New_dcmplx();
  check_mem(vectSS, "vectSS");
  vectSS->Alloc(mdimSS, vectSS);
  
  OprSSvectSS = vect_New_dcmplx();
  check_mem(OprSSvectSS, "OprSSvectSS");
  OprSSvectSS->Alloc(ndimSS, OprSSvectSS);
  
  for(idimBR = 0; idimBR < ndimBR; ++idimBR) {
    for(idimBL = 0; idimBL < ndimBL; ++idimBL) {
      for(idimSR = 0; idimSR < mdimSR; ++idimSR) {
        for(idimSL = 0; idimSL < mdimSL; ++idimSL) {
          idimBB = dmrg_indxBB(ndimBL, mdimSL, ndimBR, mdimSR, idimBL, idimSL, idimBR, idimSR);
          idimSS = idimSR*mdimSL + idimSL;
          vectSS->addr[idimSS] = vect->addr[idimBB];
        }
      }
      lapack_zgemm(ndimSS, 1, mdimSS, 'N', 'N', 1.0, OprSS->addr, vectSS->addr, 0.0, OprSSvectSS->addr); 
      for(idimSR = 0; idimSR < ndimSR; ++idimSR) {
        for(idimSL = 0; idimSL < ndimSL; ++idimSL) {
          idimBB = dmrg_indxBB(ndimBL, ndimSL, ndimBR, ndimSR, idimBL, idimSL, idimBR, idimSR);
          idimSS = idimSR*ndimSL + idimSL;
          OprSSvect->addr[idimBB] = OprSSvectSS->addr[idimSS];
        }
      }
    }
  }
  
  vect_Del_dcmplx(OprSSvectSS);
  vect_Del_dcmplx(vectSS);
  
  return;
  
 error:
  abort();
}

static void tddmrg_Superblock(int indxBL, int indxSL, int indxBR, int indxSR)
{
  int   iflav;
  dreal EE, UU, tL, tR, tC;
  
  NdimBL = Mretain[indxBL];
  NdimBR = Mretain[indxBR];
  NdimSL = Nsite;
  NdimSR = Nsite;
  NdimBB = NdimBL * NdimSL * NdimBR * NdimSR;
  NdimL  = NdimBL * NdimSL;
  NdimR  = NdimBR * NdimSR;
  NdimSS = NdimSL * NdimSR;
  
  if(HsiteL) HsiteL->Deall(HsiteL);
  if(HsiteR) HsiteR->Deall(HsiteR);
  if(Vsitesite) Vsitesite->Deall(Vsitesite);
  
  HsiteL->Alloc(NdimSL, NdimSL, HsiteL);
  HsiteR->Alloc(NdimSR, NdimSR, HsiteR);
  Vsitesite->Alloc(NdimSS, NdimSS, Vsitesite);
  
  for(iflav = 0; iflav < Nflavor; ++iflav) {
    if(SZsiteL[iflav]) SZsiteL[iflav]->Deall(SZsiteL[iflav]);
    if(SZsiteR[iflav]) SZsiteR[iflav]->Deall(SZsiteR[iflav]);
    if(SPsiteL[iflav]) SPsiteL[iflav]->Deall(SPsiteL[iflav]);
    if(SPsiteR[iflav]) SPsiteR[iflav]->Deall(SPsiteR[iflav]);
    
    SZsiteL[iflav]->Alloc(NdimSL, NdimSL, SZsiteL[iflav]);
    SZsiteR[iflav]->Alloc(NdimSR, NdimSR, SZsiteR[iflav]);
    SPsiteL[iflav]->Alloc(NdimSL, NdimSL, SPsiteL[iflav]);
    SPsiteR[iflav]->Alloc(NdimSR, NdimSR, SPsiteR[iflav]);
  }
  
  dmrg_Ssite(SZsiteL, dmrg_SigmaZ); dmrg_Ssite(SZsiteR, dmrg_SigmaZ);
  dmrg_Ssite(SPsiteL, dmrg_SigmaP); dmrg_Ssite(SPsiteR, dmrg_SigmaP);
  
  EE = Echain[indxSL]; UU = Uchain[indxSL];
  dmrg_Hmat(EE, UU, SZsiteL, HsiteL);
  
  EE = Echain[indxSR]; UU = Uchain[indxSR];
  dmrg_Hmat(EE, UU, SZsiteR, HsiteR);
  
  tC = Tchain[indxSL];
  dmrg_Vmat(tC, NdimSL, NdimSR, SPsiteL, SPsiteR, Vsitesite);
  
  if(indxBL == 0) {
    if(HblockL) HblockL->Deall(HblockL);
    if(VblocksiteL) VblocksiteL->Deall(VblocksiteL);
    
    HblockL->Alloc(NdimBL, NdimBL, HblockL);
    VblocksiteL->Alloc(NdimL, NdimL, VblocksiteL);
    
    for(iflav = 0; iflav < Nflavor; ++iflav) {
      if(SZblockL[iflav]) SZblockL[iflav]->Deall(SZblockL[iflav]);
      if(SPblockL[iflav]) SPblockL[iflav]->Deall(SPblockL[iflav]);
      SZblockL[iflav]->Alloc(NdimBL, NdimBL, SZblockL[iflav]);
      SPblockL[iflav]->Alloc(NdimBL, NdimBL, SPblockL[iflav]);
    }
    dmrg_Ssite(SZblockL, dmrg_SigmaZ);
    dmrg_Ssite(SPblockL, dmrg_SigmaP); 
    
    EE = Echain[indxBL]; UU = Uchain[indxBL];
    dmrg_Hmat(EE, UU, SZblockL, HblockL);
    
    tL = Tchain[indxBL];
    dmrg_Vmat(tL, NdimBL, NdimSL, SPblockL, SPsiteL, VblocksiteL);
  }
  
  if(indxBR == 0) {
    if(HblockR) HblockR->Deall(HblockR);
    if(VblocksiteR) VblocksiteR->Deall(VblocksiteR);
    
    HblockR->Alloc(NdimBR, NdimBR, HblockR);
    VblocksiteR->Alloc(NdimR, NdimR, VblocksiteR);
    
    for(iflav = 0; iflav < Nflavor; ++iflav) {
      if(SZblockR[iflav]) SZblockR[iflav]->Deall(SZblockR[iflav]);
      if(SPblockR[iflav]) SPblockR[iflav]->Deall(SPblockR[iflav]);
      SZblockR[iflav]->Alloc(NdimBR, NdimBR, SZblockR[iflav]);
      SPblockR[iflav]->Alloc(NdimBR, NdimBR, SPblockR[iflav]);
    }
    dmrg_Ssite(SZblockR, dmrg_SigmaZ);
    dmrg_Ssite(SPblockR, dmrg_SigmaP); 
    
    EE = Echain[Length-indxBR-1]; UU = Uchain[Length-indxBR-1]; // NOTE: crucial!
    dmrg_Hmat(EE, UU, SZblockR, HblockR);
    
    tR = Tchain[indxBR];
    dmrg_Vmat(tR, NdimBR, NdimSR, SPblockR, SPsiteR, VblocksiteR);
  }
  
  if(TDvc->ndim != NdimBB) {
    TDvc->Deall(TDvc);
    TDvc->Alloc(NdimBB, TDvc);
  }
  
  Dmat_L2R->Deall(Dmat_L2R);
  Dmat_L2R->Alloc(NdimL, NdimL, Dmat_L2R);
  
  Dmat_R2L->Deall(Dmat_R2L);
  Dmat_R2L->Alloc(NdimR, NdimR, Dmat_R2L);
}

static void tddmrg_TDvc_Predict_L2R(int indxBL, int indxBR, const vect_dcmplx *EvcGS_last, vect_dcmplx *EvcGS_next)
{
  int         ndimBL_next, ndimSL_next, ndimBR_next, ndimSR_next, ndimL_next, ndimR_next,
              ndimBL_last, ndimSL_last, ndimBR_last, ndimSR_last, ndimL_last, ndimR_last,
              ndimBB_contr;
  int         idimBL_next, idimSL_next, idimBR_next, idimSR_next, idimR_next, idimBB_next,
              idimBR_last, idimSR_last, idimR_last, idimBB_contr;
  dcmplx      *Oprptr = NULL, *vectptr = NULL, *Oprvectptr = NULL;
  vect_dcmplx *vectBR_last = NULL, *vectR_next = NULL, *vect_contr = NULL;
  
  check(indxBL > 0, "Invalid indxBL");
  check(indxBR < Nblock-1, "Invalid indxBL");
  
  ndimBL_next = Mretain[indxBL];
  ndimBR_next = Mretain[indxBR];
  ndimSL_next = Nsite;
  ndimSR_next = Nsite;
  ndimL_next  = ndimBL_next*ndimSL_next;
  ndimR_next  = ndimBR_next*ndimSR_next;
  
  ndimBL_last = Mretain[indxBL-1];
  ndimBR_last = Mretain[indxBR+1];
  ndimSL_last = Nsite;
  ndimSR_last = Nsite;
  ndimL_last  = ndimBL_last*ndimSL_last;
  ndimR_last  = ndimBR_last*ndimSR_last;
  
  ndimBB_contr= ndimBL_next*ndimBR_last*ndimSR_last;
  
  vectBR_last = vect_New_dcmplx();
  check_mem(vectBR_last, "vectBR_last");
  vectBR_last->Alloc(ndimBR_last, vectBR_last);
  
  vectR_next = vect_New_dcmplx();
  check_mem(vectR_next, "vectR_next");
  vectR_next->Alloc(ndimR_next, vectR_next);
  
  vect_contr = vect_New_dcmplx();
  check_mem(vect_contr, "vect_contr");
  vect_contr->Alloc(ndimBB_contr, vect_contr);
  
  Oprptr = zWblock_L2R[indxBL-1]->addr;
  for(idimSR_last = 0; idimSR_last < ndimSR_last; ++idimSR_last) {
    for(idimBR_last = 0; idimBR_last < ndimBR_last; ++idimBR_last) {
      idimR_last = idimSR_last*ndimBR_last + idimBR_last;
      vectptr    = EvcGS_last->addr + idimR_last*ndimL_last;
      Oprvectptr = vect_contr->addr + idimR_last*ndimBL_next;
      lapack_zgemm(ndimBL_next, 1, ndimL_last, 'C', 'N', 1.0, Oprptr, vectptr, 0.0, Oprvectptr);
    }
  }
  
  Oprptr = zWblock_R2L[indxBR]->addr;
  for(idimSR_last = 0; idimSR_last < ndimSR_last; ++idimSR_last) {
    idimSL_next = idimSR_last;
    for(idimBL_next = 0; idimBL_next < ndimBL_next; ++idimBL_next) {
      for(idimBR_last = 0; idimBR_last < ndimBR_last; ++idimBR_last) {
        idimBB_contr = dmrg_indxBB(ndimBL_next, 1, ndimBR_last, ndimSR_last,
                                   idimBL_next, 0, idimBR_last, idimSR_last);
        vectBR_last->addr[idimBR_last] = vect_contr->addr[idimBB_contr];
      }
      vectptr    = vectBR_last->addr;
      Oprvectptr = vectR_next->addr;
      lapack_zgemm(ndimR_next, 1, ndimBR_last, 'N', 'N', 1.0, Oprptr, vectptr, 0.0, Oprvectptr);
      for(idimSR_next = 0; idimSR_next < ndimSR_next; ++idimSR_next) {
        for(idimBR_next = 0; idimBR_next < ndimBR_next; ++idimBR_next) {
          idimR_next  = idimSR_next*ndimBR_next + idimBR_next;
          idimBB_next = dmrg_indxBB(ndimBL_next, ndimSL_next, ndimBR_next, ndimSR_next,
                                    idimBL_next, idimSL_next, idimBR_next, idimSR_next);
          EvcGS_next->addr[idimBB_next] = vectR_next->addr[idimR_next];
        }
      }
    }
  }
  
  vect_Del_dcmplx(vect_contr);
  vect_Del_dcmplx(vectR_next);
  vect_Del_dcmplx(vectBR_last);
  
  return;
  
 error:
  abort();
}

static void tddmrg_TDvc_Predict_R2L(int indxBL, int indxBR, const vect_dcmplx *TDvc_last, vect_dcmplx *TDvc_next)
{
  int         ndimBL_next, ndimSL_next, ndimBR_next, ndimSR_next, ndimL_next, ndimR_next,
              ndimBL_last, ndimSL_last, ndimBR_last, ndimSR_last, ndimL_last, ndimR_last,
              ndimBB_contr;
  int         idimBR_next, idimSR_next, idimR_next,
              idimBL_last, idimSL_last, idimBR_last, idimSR_last, idimR_last, idimBB_last, idimBB_contr;
  dcmplx      *Oprptr = NULL, *vectptr = NULL, *Oprvectptr = NULL;
  vect_dcmplx *vectR_last = NULL, *vectBR_next = NULL, *vect_contr = NULL;
  
  ndimBL_next = Mretain[indxBL];
  ndimBR_next = Mretain[indxBR];
  ndimSL_next = Nsite;
  ndimSR_next = Nsite;
  ndimL_next  = ndimBL_next*ndimSL_next;
  ndimR_next  = ndimBR_next*ndimSR_next;
  
  ndimBL_last = Mretain[indxBL+1];
  ndimBR_last = Mretain[indxBR-1];
  ndimSL_last = Nsite;
  ndimSR_last = Nsite;
  ndimL_last  = ndimBL_last*ndimSL_last;
  ndimR_last  = ndimBR_last*ndimSR_last;
  
  ndimBB_contr = ndimBL_last*ndimBR_next*ndimSR_next;
  
  vectR_last = vect_New_dcmplx();
  check_mem(vectR_last, "vectR_last");
  vectR_last->Alloc(ndimR_last, vectR_last);
  
  vectBR_next = vect_New_dcmplx();
  check_mem(vectBR_next, "vectR_next");
  vectBR_next->Alloc(ndimBR_next, vectBR_next);
  
  vect_contr = vect_New_dcmplx();
  check_mem(vect_contr, "vect_contr");
  vect_contr->Alloc(ndimBB_contr, vect_contr);
  
  Oprptr = zWblock_R2L[indxBR-1]->addr;
  for(idimSL_last = 0; idimSL_last < ndimSL_last; ++idimSL_last) {
    idimSR_next = idimSL_last;
    for(idimBL_last = 0; idimBL_last < ndimBL_last; ++idimBL_last) {
      for(idimSR_last = 0; idimSR_last < ndimSR_last; ++idimSR_last) {
        for(idimBR_last = 0; idimBR_last < ndimBR_last; ++idimBR_last) {
          idimR_last = idimSR_last*ndimBR_last + idimBR_last;
          idimBB_last = dmrg_indxBB(ndimBL_last, ndimSL_last, ndimBR_last, ndimSR_last,
                                    idimBL_last, idimSL_last, idimBR_last, idimSR_last);
          vectR_last->addr[idimR_last] = TDvc_last->addr[idimBB_last];
        }
      }
      lapack_zgemm(ndimBR_next, 1, ndimR_last, 'C', 'N', 1.0, Oprptr, vectR_last->addr, 0.0, vectBR_next->addr);
      for(idimBR_next = 0; idimBR_next < ndimBR_next; ++idimBR_next) {
        idimBB_contr = dmrg_indxBB(ndimBL_last, 1, ndimBR_next, ndimSR_next,
                                   idimBL_last, 0, idimBR_next, idimSR_next);
        vect_contr->addr[idimBB_contr] = vectBR_next->addr[idimBR_next];
      }
    }
  }
  
  Oprptr = zWblock_L2R[indxBL]->addr;
  for(idimSR_next = 0; idimSR_next < ndimSR_next; ++idimSR_next) {
    for(idimBR_next = 0; idimBR_next < ndimBR_next; ++idimBR_next) {
      idimR_next = idimSR_next*ndimBR_next + idimBR_next;
      vectptr = vect_contr->addr + idimR_next*ndimBL_last;
      Oprvectptr = TDvc_next->addr + idimR_next*ndimL_next;
      lapack_zgemm(ndimL_next, 1, ndimBL_last, 'N', 'N', 1.0, Oprptr, vectptr, 0.0, Oprvectptr);
    }
  }
  
  vect_Del_dcmplx(vect_contr);
  vect_Del_dcmplx(vectBR_next);
  vect_Del_dcmplx(vectR_last);
  
  return;
  
 error:
  abort();
}

static void tddmrg_Evolut_L2R(int indxBL, int indxBR, vect_dcmplx *TDvc)
{
  int          ndimBB;
  vect_dcmplx  *OprTDvc = NULL;
  mat2d_dcmplx *expHblockL = NULL, *expHblockR = NULL, *expHsiteL = NULL, *expHsiteR = NULL,
               *expVblocksiteL = NULL, *expVblocksiteR = NULL, *expVsitesite = NULL;
  
  ndimBB = TDvc->ndim;
  
  OprTDvc = vect_New_dcmplx();
  check_mem(OprTDvc, "OprTDvc");
  OprTDvc->Alloc(ndimBB, OprTDvc);
  
  if(indxBL == 0) {
    expHblockL = mat2d_New_dcmplx();
    check_mem(expHblockL, "expHblockL");
    expHblockL->Alloc(NdimBL, NdimBL, expHblockL);
    
    expVblocksiteL = mat2d_New_dcmplx();
    check_mem(expVblocksiteL, "expVblocksiteL");
    expVblocksiteL->Alloc(NdimL, NdimL, expVblocksiteL);
    
    tddmrg_expMat(0.5*dT, HblockL, expHblockL);
    tddmrg_expMat(0.5*dT, VblocksiteL, expVblocksiteL);
    
    tddmrg_OprBLvect(NdimBL, NdimSL, NdimBR, NdimSR, NdimBL, expHblockL, TDvc, OprTDvc);
    vect_Copy_dcmplx(OprTDvc, TDvc);
    tddmrg_OprLvect(NdimBL, NdimSL, NdimBR, NdimSR, NdimBL, NdimSL, expVblocksiteL, TDvc, OprTDvc);
    vect_Copy_dcmplx(OprTDvc, TDvc);
    
    mat2d_Del_dcmplx(expVblocksiteL);
    mat2d_Del_dcmplx(expHblockL);
  }
  
  expHsiteL = mat2d_New_dcmplx();
  check_mem(expHsiteL, "expHsiteL");
  expHsiteL->Alloc(NdimSL, NdimSL, expHsiteL);
  
  expVsitesite = mat2d_New_dcmplx();
  check_mem(expVsitesite, "expVsitesite");
  expVsitesite->Alloc(NdimSS, NdimSS, expVsitesite);
  
  tddmrg_expMat(0.5*dT, HsiteL, expHsiteL);
  tddmrg_expMat(0.5*dT, Vsitesite, expVsitesite);
  
  tddmrg_OprSLvect(NdimBL, NdimSL, NdimBR, NdimSR, NdimSL, expHsiteL, TDvc, OprTDvc);
  vect_Copy_dcmplx(OprTDvc, TDvc);
  tddmrg_OprSSvect(NdimBL, NdimSL, NdimBR, NdimSR, NdimSL, NdimSR, expVsitesite, TDvc, OprTDvc);
  vect_Copy_dcmplx(OprTDvc, TDvc);
  
  mat2d_Del_dcmplx(expVsitesite);
  mat2d_Del_dcmplx(expHsiteL);
  
  if(indxBR == 0) {
    expHsiteR = mat2d_New_dcmplx();
    check_mem(expHsiteR, "expHsiteR");
    expHsiteR->Alloc(NdimSR, NdimSR, expHsiteR);
    
    expVblocksiteR = mat2d_New_dcmplx();
    check_mem(expVblocksiteR, "expVblocksiteR");
    expVblocksiteR->Alloc(NdimR, NdimR, expVblocksiteR);
    
    expHblockR = mat2d_New_dcmplx();
    check_mem(expHblockR, "expHblockR");
    expHblockR->Alloc(NdimBR, NdimBR, expHblockR);
    
    tddmrg_expMat(0.5*dT, HsiteR, expHsiteR);
    tddmrg_expMat(0.5*dT, VblocksiteR, expVblocksiteR);
    tddmrg_expMat(0.5*dT, HblockR, expHblockR);
    
    tddmrg_OprSRvect(NdimBL, NdimSL, NdimBR, NdimSR, NdimSR, expHsiteR, TDvc, OprTDvc);
    vect_Copy_dcmplx(OprTDvc, TDvc);
    tddmrg_OprRvect(NdimBL, NdimSL, NdimBR, NdimSR, NdimBR, NdimSR, expVblocksiteR, TDvc, OprTDvc);
    vect_Copy_dcmplx(OprTDvc, TDvc);
    tddmrg_OprBRvect(NdimBL, NdimSL, NdimBR, NdimSR, NdimBR, expHblockR, TDvc, OprTDvc);
    vect_Copy_dcmplx(OprTDvc, TDvc);
    
    mat2d_Del_dcmplx(expHblockR);
    mat2d_Del_dcmplx(expVblocksiteR);
    mat2d_Del_dcmplx(expHsiteR);
  }
  
  vect_Del_dcmplx(OprTDvc);
  
  return;
  
 error:
  abort();
}

static void tddmrg_Evolut_R2L(int indxBL, int indxBR, vect_dcmplx *TDvc)
{
  int          ndimBB;
  vect_dcmplx  *OprTDvc = NULL;
  mat2d_dcmplx *expHblockL = NULL, *expHblockR = NULL, *expHsiteL = NULL, *expHsiteR = NULL,
               *expVblocksiteL = NULL, *expVblocksiteR = NULL, *expVsitesite = NULL;
  
  ndimBB = TDvc->ndim;
  
  OprTDvc = vect_New_dcmplx();
  check_mem(OprTDvc, "OprTDvc");
  OprTDvc->Alloc(ndimBB, OprTDvc);
  
  if(indxBR == 0) {
    expHblockR = mat2d_New_dcmplx();
    check_mem(expHblockR, "expHblockR");
    expHblockR->Alloc(NdimBR, NdimBR, expHblockR);
    
    expVblocksiteR = mat2d_New_dcmplx();
    check_mem(expVblocksiteR, "expVblocksiteR");
    expVblocksiteR->Alloc(NdimR, NdimR, expVblocksiteR);
    
    expHsiteR = mat2d_New_dcmplx();
    check_mem(expHsiteR, "expHsiteR");
    expHsiteR->Alloc(NdimSR, NdimSR, expHsiteR);
    
    tddmrg_expMat(0.5*dT, HblockR, expHblockR);
    tddmrg_expMat(0.5*dT, VblocksiteR, expVblocksiteR);
    tddmrg_expMat(0.5*dT, HsiteR, expHsiteR);
    
    tddmrg_OprBRvect(NdimBL, NdimSL, NdimBR, NdimSR, NdimBR, expHblockR, TDvc, OprTDvc);
    vect_Copy_dcmplx(OprTDvc, TDvc);
    tddmrg_OprRvect(NdimBL, NdimSL, NdimBR, NdimSR, NdimBR, NdimSR, expVblocksiteR, TDvc, OprTDvc);
    vect_Copy_dcmplx(OprTDvc, TDvc);
    tddmrg_OprSRvect(NdimBL, NdimSL, NdimBR, NdimSR, NdimSR, expHsiteR, TDvc, OprTDvc);
    vect_Copy_dcmplx(OprTDvc, TDvc);
    
    mat2d_Del_dcmplx(expHsiteR);
    mat2d_Del_dcmplx(expVblocksiteR);
    mat2d_Del_dcmplx(expHblockR);
  }
  
  expVsitesite = mat2d_New_dcmplx();
  check_mem(expVsitesite, "expVsitesite");
  expVsitesite->Alloc(NdimSS, NdimSS, expVsitesite);
  
  expHsiteL = mat2d_New_dcmplx();
  check_mem(expHsiteL, "expHsiteL");
  expHsiteL->Alloc(NdimSL, NdimSL, expHsiteL);
  
  tddmrg_expMat(0.5*dT, Vsitesite, expVsitesite);
  tddmrg_expMat(0.5*dT, HsiteL, expHsiteL);
  
  tddmrg_OprSSvect(NdimBL, NdimSL, NdimBR, NdimSR, NdimSL, NdimSR, expVsitesite, TDvc, OprTDvc);
  vect_Copy_dcmplx(OprTDvc, TDvc);
  tddmrg_OprSLvect(NdimBL, NdimSL, NdimBR, NdimSR, NdimSL, expHsiteL, TDvc, OprTDvc);
  vect_Copy_dcmplx(OprTDvc, TDvc);
  
  mat2d_Del_dcmplx(expHsiteL);
  mat2d_Del_dcmplx(expVsitesite);
  
  if(indxBL == 0) {
    expVblocksiteL = mat2d_New_dcmplx();
    check_mem(expVblocksiteL, "expVblocksiteL");
    expVblocksiteL->Alloc(NdimL, NdimL, expVblocksiteL);
    
    expHblockL = mat2d_New_dcmplx();
    check_mem(expHblockL, "expHblockL");
    expHblockL->Alloc(NdimBL, NdimBL, expHblockL);
    
    tddmrg_expMat(0.5*dT, VblocksiteL, expVblocksiteL);
    tddmrg_expMat(0.5*dT, HblockL, expHblockL);
    
    tddmrg_OprLvect(NdimBL, NdimSL, NdimBR, NdimSR, NdimBL, NdimSL, expVblocksiteL, TDvc, OprTDvc);
    vect_Copy_dcmplx(OprTDvc, TDvc);
    tddmrg_OprBLvect(NdimBL, NdimSL, NdimBR, NdimSR, NdimBL, expHblockL, TDvc, OprTDvc);
    vect_Copy_dcmplx(OprTDvc, TDvc);
    
    mat2d_Del_dcmplx(expHblockL);
    mat2d_Del_dcmplx(expVblocksiteL);
  }
  
  vect_Del_dcmplx(OprTDvc);
  
  return;
  
 error:
  abort();
}

static void tddmrg_Dmat_L2R(void)
{
  int idimBL, jdimBL, idimSL, jdimSL, idimL, jdimL, idimBB, jdimBB, idimSR, idimBR;
  
  Dmat_L2R->Reset(Dmat_L2R);
  
  for(idimSL = 0; idimSL < NdimSL; ++idimSL) {
    for(idimBL = 0; idimBL < NdimBL; ++idimBL) {
      idimL = idimSL*NdimBL + idimBL;
      
      for(jdimSL = 0; jdimSL < NdimSL; ++jdimSL) {
        for(jdimBL = 0; jdimBL < NdimBL; ++jdimBL) {
          jdimL = jdimSL*NdimBL + jdimBL;
          
          for(idimSR = 0; idimSR < NdimSR; ++idimSR) {
            for(idimBR = 0; idimBR < NdimBR; ++idimBR) {
              idimBB = dmrg_indxBB(NdimBL, NdimSL, NdimBR, NdimSR, idimBL, idimSL, idimBR, idimSR);
              jdimBB = dmrg_indxBB(NdimBL, NdimSL, NdimBR, NdimSR, jdimBL, jdimSL, idimBR, idimSR);
              Dmat_L2R->ptr[idimL][jdimL] += TDvc->addr[jdimBB] * conj(TDvc->addr[idimBB]);
            }
          }
        }
      }
    }
  }
}

static void tddmrg_Dmat_R2L(void)
{
  int idimBR, jdimBR, idimSR, jdimSR, idimR, jdimR, idimBB, jdimBB, idimSL, idimBL;
  
  Dmat_R2L->Reset(Dmat_R2L);
  
  for(idimSR = 0; idimSR < NdimSR; ++idimSR) {
    for(idimBR = 0; idimBR < NdimBR; ++idimBR) {
      idimR = idimSR*NdimBR + idimBR;
      
      for(jdimSR = 0; jdimSR < NdimSR; ++jdimSR) {
        for(jdimBR = 0; jdimBR < NdimBR; ++jdimBR) {
          jdimR = jdimSR*NdimBR + jdimBR;
          
          for(idimSL = 0; idimSL < NdimSL; ++idimSL) {
            for(idimBL = 0; idimBL < NdimBL; ++idimBL) {
              idimBB = dmrg_indxBB(NdimBL, NdimSL, NdimBR, NdimSR, idimBL, idimSL, idimBR, idimSR);
              jdimBB = dmrg_indxBB(NdimBL, NdimSL, NdimBR, NdimSR, idimBL, idimSL, jdimBR, jdimSR);
              Dmat_R2L->ptr[idimR][jdimR] += TDvc->addr[jdimBB] * conj(TDvc->addr[idimBB]);
            }
          }
        }
      }
    }
  }
}

static void tddmrg_zWblock_L2R(mat2d_dcmplx *Wvc)
{
  int    ndim, mdim, idim, idimL;
  dreal  error;
  size_t nlen;
  dreal  *eig = NULL;
  
  ndim = Wvc->nrow;
  mdim = Wvc->ncol;
  
  nlen = ndim*sizeof(dcmplx);
  
  eig = (dreal *) calloc(NdimL, sizeof(dreal));
  
  lapack_zheev(NdimL, Dmat_L2R->addr, eig);
  
  for(idimL = 0, error = 0.0; idimL < NdimL-mdim; ++idimL)
    error += eig[idimL];
  
  //if(Rank == Root) {
  //  printf("%s\n", tddmrglabel);
  //  printf("%s error = %15.5e\n", tddmrglabel, error);
  //}
  
  for(idim = 0; idim < mdim; ++idim) {
    idimL = NdimL-idim-1;
    memcpy(Wvc->ptr[idim], Dmat_L2R->ptr[idimL], nlen);
  }
  
  freeup(eig);
}

static void tddmrg_zWblock_R2L(mat2d_dcmplx *Wvc)
{
  int    ndim, mdim, idim, idimR;
  dreal  error;
  size_t nlen;
  dreal  *eig = NULL;
  
  ndim = Wvc->nrow;
  mdim = Wvc->ncol;
  
  nlen = ndim*sizeof(dcmplx);
  
  eig = (dreal *) calloc(NdimR, sizeof(dreal));
  
  lapack_zheev(NdimR, Dmat_R2L->addr, eig);
  
  for(idimR = 0, error = 0.0; idimR < NdimR-mdim; ++idimR)
    error += eig[idimR];
  
  //if(Rank == Root) {
  //  printf("%s\n", tddmrglabel);
  //  printf("%s error = %15.5e\n", tddmrglabel, error);
  //}
  
  for(idim = 0; idim < mdim; ++idim) {
    idimR = NdimR-idim-1;
    memcpy(Wvc->ptr[idim], Dmat_R2L->ptr[idimR], nlen);
  }
  
  freeup(eig);
}

static void tddmrg_Left2Right(void)
{
  int         iblock, indxBL, indxBR, indxSL, indxSR;
  vect_dcmplx *TDvc_old = NULL;
  
  TDvc_old = vect_New_dcmplx();
  
  for(iblock = 0; iblock < Nblock; ++iblock) {
    indxBL = iblock;   indxBR = Length-iblock-4;
    indxSL = iblock+1; indxSR = iblock+2;
    tddmrg_Superblock(indxBL, indxSL, indxBR, indxSR);
    if(indxBL > 0) tddmrg_TDvc_Predict_L2R(indxBL, indxBR, TDvc_old, TDvc);
    tddmrg_Evolut_L2R(indxBL, indxBR, TDvc);
    tddmrg_TDvc_Save(TDvc, TDvc_old);
    tddmrg_Dmat_L2R();
    tddmrg_zWblock_L2R(zWblock_L2R[indxBL]);
  }
  
  vect_Del_dcmplx(TDvc_old);
}

static void tddmrg_Right2Left(void)
{
  int         iblock, indxBL, indxBR, indxSL, indxSR;
  vect_dcmplx *TDvc_old = NULL;
  
  TDvc_old = vect_New_dcmplx();
  
  for(iblock = 0; iblock < Nblock; ++iblock) {
    indxBL = Length-iblock-4; indxBR = iblock;
    indxSL = Length-iblock-3; indxSR = Length-iblock-2;
    tddmrg_Superblock(indxBL, indxSL, indxBR, indxSR);
    if(indxBR > 0) tddmrg_TDvc_Predict_R2L(indxBL, indxBR, TDvc_old, TDvc);
    tddmrg_Evolut_R2L(indxBL, indxBR, TDvc);
    tddmrg_TDvc_Save(TDvc, TDvc_old);
    tddmrg_Dmat_R2L();
    tddmrg_zWblock_R2L(zWblock_R2L[indxBR]);
  }
  
  vect_Del_dcmplx(TDvc_old);
}

static void tddmrg_CurrentL(void)
{
  int          ndimBL, ndimSL, ndimBR, ndimSR, ndimSS, ndimBB, indxBL, indxSL, indxBR, indxSR,
               NblockHalf, iblock, iflav;
  dreal        tC;
  vect_dcmplx  *TDvc_old = NULL, *TDvc_new = NULL;
  mat2d_dcmplx *JLmat = NULL;
  mat2d_dreal  **spsiteL = NULL, **spsiteR = NULL;
  
  ndimBB = TDvc->ndim;
  
  TDvc_old = vect_New_dcmplx();
  TDvc_new = vect_New_dcmplx();
  
  TDvc_old->Alloc(ndimBB, TDvc_old);
  vect_Copy_dcmplx(TDvc, TDvc_old);
  
  NblockHalf = (Length - Length%2) / 2 - 1;
  
  //printf("CurrentL: NblockHalf = %d\n", NblockHalf);
  
  for(iblock = 0; iblock < NblockHalf; ++iblock) {
    indxBL = iblock;   indxBR = Length-iblock-4;
    indxSL = iblock+1; indxSR = iblock+2;
    
    ndimBL = Mretain[indxBL];
    ndimBR = Mretain[indxBR];
    ndimSL = Nsite;
    ndimSR = Nsite;
    ndimBB = ndimBL * ndimSL * ndimBR * ndimSR;
    
    TDvc_new->Deall(TDvc_new);
    TDvc_new->Alloc(ndimBB, TDvc_new);
    
    if(indxBL == 0) {
      vect_Copy_dcmplx(TDvc_old, TDvc_new);
    } else {
      tddmrg_TDvc_Predict_L2R(indxBL, indxBR, TDvc_old, TDvc_new);
    }
    
    tddmrg_TDvc_Save(TDvc_new, TDvc_old);
  }
  
  //printf("CurrentL: indxBL, indxSL, indxBR, indxSR = %d %d %d %d\n", indxBL, indxSL, indxBR, indxSR);
  indxBL = NblockHalf-1; indxBR = Length-NblockHalf-3;
  indxSL = NblockHalf;   indxSR = NblockHalf+1;
  //printf("CurrentL: indxBL, indxSL, indxBR, indxSR = %d %d %d %d\n", indxBL, indxSL, indxBR, indxSR);
  
  ndimBL = Mretain[indxBL];
  ndimBR = Mretain[indxBR];
  ndimSL = Nsite;
  ndimSR = Nsite;
  ndimSS = ndimSL*ndimSR;
  
  JLmat = mat2d_New_dcmplx();
  check_mem(JLmat, "JLmat");
  JLmat->Alloc(ndimSS, ndimSS, JLmat);
  
  tC = Tchain[indxSL];
  
  spsiteL = (mat2d_dreal **) calloc(Nflavor, sizeof(mat2d_dreal *));
  spsiteR = (mat2d_dreal **) calloc(Nflavor, sizeof(mat2d_dreal *));
  for(iflav = 0; iflav < Nflavor; ++iflav) {
    spsiteL[iflav] = mat2d_New_dreal();
    spsiteR[iflav] = mat2d_New_dreal();
    spsiteL[iflav]->Alloc(Nsite, Nsite, spsiteL[iflav]);
    spsiteR[iflav]->Alloc(Nsite, Nsite, spsiteR[iflav]);
  }
  
  dmrg_Ssite(spsiteL, dmrg_SigmaP); dmrg_Ssite(spsiteR, dmrg_SigmaP);
  
  tddmrg_Jmat(tC, ndimSL, ndimSR, spsiteL, spsiteR, JLmat);
  
  tddmrg_OprSSvect(ndimBL, ndimSL, ndimBR, ndimSR, ndimSL, ndimSR, JLmat, TDvc_new, TDvc_old);
  
  JL = vect_Dotprod_dcmplx(TDvc_new, TDvc_old);
  
  for(iflav = Nflavor-1; iflav > -1; --iflav) {
    mat2d_Del_dreal(spsiteR[iflav]);
    mat2d_Del_dreal(spsiteL[iflav]);
  }
  freeup(spsiteR); freeup(spsiteL);
  
  mat2d_Del_dcmplx(JLmat);
  
  vect_Del_dcmplx(TDvc_new);
  vect_Del_dcmplx(TDvc_old);
  
  return;
  
 error:
  abort();
}

static void tddmrg_CurrentR(void)
{
  int          ndimBL, ndimSL, ndimBR, ndimSR, ndimSS, ndimBB, indxBL, indxSL, indxBR, indxSR,
               NblockHalf, iblock, iflav;
  dreal        tC;
  vect_dcmplx  *TDvc_old = NULL, *TDvc_new = NULL;
  mat2d_dcmplx *JRmat = NULL;
  mat2d_dreal  **spsiteL = NULL, **spsiteR = NULL;
  
  ndimBB = TDvc->ndim;
  
  TDvc_old = vect_New_dcmplx();
  TDvc_new = vect_New_dcmplx();
  
  TDvc_old->Alloc(ndimBB, TDvc_old);
  vect_Copy_dcmplx(TDvc, TDvc_old);
  
  NblockHalf = Length - ( (Length - Length%2) / 2 - 1 ) - 2;
  
  for(iblock = 0; iblock < NblockHalf; ++iblock) {
    indxBL = iblock;   indxBR = Length-iblock-4;
    indxSL = iblock+1; indxSR = iblock+2;
    
    ndimBL = Mretain[indxBL];
    ndimBR = Mretain[indxBR];
    ndimSL = Nsite;
    ndimSR = Nsite;
    ndimBB = ndimBL * ndimSL * ndimBR * ndimSR;
    
    TDvc_new->Deall(TDvc_new);
    TDvc_new->Alloc(ndimBB, TDvc_new);
    
    if(indxBL == 0) {
      vect_Copy_dcmplx(TDvc_old, TDvc_new);
    } else {
      tddmrg_TDvc_Predict_L2R(indxBL, indxBR, TDvc_old, TDvc_new);
    }
    
    tddmrg_TDvc_Save(TDvc_new, TDvc_old);
  }
  
  indxBL = NblockHalf-1; indxBR = Length-NblockHalf-3;
  indxSL = NblockHalf;   indxSR = NblockHalf+1;
  
  ndimBL = Mretain[indxBL];
  ndimBR = Mretain[indxBR];
  ndimSL = Nsite;
  ndimSR = Nsite;
  ndimSS = ndimSL*ndimSR;
  
  JRmat = mat2d_New_dcmplx();
  check_mem(JRmat, "JRmat");
  JRmat->Alloc(ndimSS, ndimSS, JRmat);
  
  tC = Tchain[indxSL];
  
  spsiteL = (mat2d_dreal **) calloc(Nflavor, sizeof(mat2d_dreal *));
  spsiteR = (mat2d_dreal **) calloc(Nflavor, sizeof(mat2d_dreal *));
  for(iflav = 0; iflav < Nflavor; ++iflav) {
    spsiteL[iflav] = mat2d_New_dreal();
    spsiteR[iflav] = mat2d_New_dreal();
    spsiteL[iflav]->Alloc(Nsite, Nsite, spsiteL[iflav]);
    spsiteR[iflav]->Alloc(Nsite, Nsite, spsiteR[iflav]);
  }
  
  dmrg_Ssite(spsiteL, dmrg_SigmaP); dmrg_Ssite(spsiteR, dmrg_SigmaP);
  
  tddmrg_Jmat(tC, ndimSL, ndimSR, spsiteL, spsiteR, JRmat);
  
  tddmrg_OprSSvect(ndimBL, ndimSL, ndimBR, ndimSR, ndimSL, ndimSR, JRmat, TDvc_new, TDvc_old);
  
  JR = vect_Dotprod_dcmplx(TDvc_new, TDvc_old);
  
  for(iflav = Nflavor-1; iflav > -1; --iflav) {
    mat2d_Del_dreal(spsiteR[iflav]);
    mat2d_Del_dreal(spsiteL[iflav]);
  }
  freeup(spsiteR); freeup(spsiteL);
  
  mat2d_Del_dcmplx(JRmat);
  
  vect_Del_dcmplx(TDvc_new);
  vect_Del_dcmplx(TDvc_old);
  
  return;
  
 error:
  abort();
}

static inline void tddmrg_Info(int istp, dreal time)
{
  FILE *fp = NULL;
  char mode[2];
  
  if(Rank == Root) {
    if(numeric_Ispow2(istp)) {
      printf("%s istep, Time = %5d %10.5f\n", tddmrglabel, istp, time);
      fflush(stdout);
    }
    if(istp == 0) strcpy(mode, "w");
    else strcpy(mode, "a");
    fp = fopen(tdcurrent, mode);
    fprintf(fp, "%10.5f %15.5e %15.5e\n", time, creal(JL), creal(JR));
    fclose(fp);
  }
}

void tddmrg_Propagate(void)
{
  int istep;
  
  istep = 0; Time = 0.0;
  tddmrg_CurrentL();
  tddmrg_CurrentR();
  tddmrg_Info(istep, Time);
  
  for(istep = 1; istep < Nstep; ++istep) {
    Time = dT * istep;
    Coef = erf(Alpha * (Time-0.5*dT));
    tddmrg_Echain();
    tddmrg_Left2Right();
    tddmrg_Right2Left();
    tddmrg_CurrentL();
    tddmrg_CurrentR();
    tddmrg_Info(istep, Time);
  }
}
