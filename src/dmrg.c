#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include "dbg.h"
#include "sys.h"
#include "matrix.h"
#include "vector.h"
#include "numeric.h"
#include "lapack.h"
#include "constants.h"
#include "parameters.h"
#include "davidson.h"
#include "dmrg.h"

// Density Matrix Renormalization Group Parameter list
//
// NSP:         Number of spin, always 2
// Mret:        Number of states retained in a block, if number of states exceed Mret
// Nsweep:      Number of sweeps in finite DMRG calculations
// Nflavor:     Number of flavors at each site, 1 or 2 to mimic nonspin/spin fermionic chain
// Nsite:       Number of states per site
// Nblock:      Number of blocks that will be calculated in finite DMRG
// NblockInf:   Number of blocks that will be calculated in infinite DMRG
// NdimBL:      Dimension of left block
// NdimSL:      Dimension of left site, essentially Nsite
// NdimBR:      Dimension of right block
// NdimSR:      Dimension of right site, essentially Nsite
// NdimBB:      Dimension of superblock, NdimBL*NdimSL*NdimSR*NdimBR
// NdimSS:      Dimension of central part, NdimSL*NdimSR
// NdimL:       Dimension of left part, NdimBL*NdimSL
// NdimR:       Dimension of right part, NdimBR*NdimSR
// EigGS:       Ground state energy of superblock
// EigD:        Eigenvalues of density matrix
// EvcGS:       Ground state wavefunction of superblock
// Dmat:        Density matrix
// HblockL:     Hamilnotian of left block
// HblockR:     Hamilnotian of right block
// HsiteL:      Hamiltonian of left site
// HsiteR:      Hamiltonian of right site
// VblocksiteL: Tunneling Hamiltonian for left block and left site
// VblocksiteR: Tunneling Hamiltonian for right block and left site
// Vsitesite:   Tunneling Hamiltonian for right block and left site
// SZblockL:    SigmaZ operator in left block
// SZblockR:    SigmaZ operator in right block
// SZsiteL:     SigmaZ operator in left site 
// SZsiteR:     SigmaZ operator in right site 
// SPblockL:    SigmaPlus operator in left block
// SPblockR:    SigmaPlus operator in right block
// SPsiteL:     SigmaPlus operator in left site 
// SPsiteR:     SigmaPlus operator in right site 
// Hblock:      Hamiltonian of each block
// Wblock:      Basis wavefunctions of each block
// SZblock:     SigmaZ operators of each block
// SPblock:     SigmaPlus operators of each block

const int   NSP = 2;
int         isGS, isTD;
int         Nblock, Mret, Nsweep, Nflavor, Nsite, NblockInf;
int         NdimBL, NdimSL, NdimBR, NdimSR, NdimBB, NdimSS, NdimL, NdimR;
int         *Mretain = NULL;
vect_dreal  *DenGS = NULL, *EvcGS = NULL;
mat2d_dreal **Wblock = NULL;

static dreal       EigGS;
static dreal       *EigD = NULL;
static mat2d_dreal *Dmat = NULL;
static mat2d_dreal *HblockL = NULL, *HblockR = NULL, *HsiteL = NULL, *HsiteR = NULL;
static mat2d_dreal *VblocksiteL = NULL, *VblocksiteR = NULL, *Vsitesite = NULL;
static mat2d_dreal **SZblockL = NULL, **SZblockR = NULL, **SZsiteL = NULL, **SZsiteR = NULL;
static mat2d_dreal **SPblockL = NULL, **SPblockR = NULL, **SPsiteL = NULL, **SPsiteR = NULL;
static mat2d_dreal **Hblock = NULL;
static mat2d_dreal ***SZblock = NULL, ***SPblock = NULL;

const char *evcgsbin  = "dmrgEvcGS.bin";
const char *wblockbin = "dmrgWblock.bin";
//const char *hblockbin  = "dmrgHblock.bin";
//const char *szblockbin = "dmrgSZblock.bin";
//const char *spblockbin = "dmrgSPblock.bin";

static const char *dmrgin    = "dmrg.in";
static const char *dmrglabel = "dmrg:           ";
// Read in DMRG control parameters
static void dmrg_Read(void)
{
  FILE *fp = NULL;
  char line[SHRT_MAX];
  
  if(Rank == Root) {
    fp = fopen(dmrgin, "r");
    if(fgets(line, sizeof(line), fp)) sscanf(line, "%d", &Mret);
    if(fgets(line, sizeof(line), fp)) sscanf(line, "%d", &Nsweep);
    if(fgets(line, sizeof(line), fp)) sscanf(line, "%d", &isGS);
    if(fgets(line, sizeof(line), fp)) sscanf(line, "%d", &isTD);
    fclose(fp);
  }
  sys_Bcast_int(&Mret,   1, Root);
  sys_Bcast_int(&Nsweep, 1, Root);
  sys_Bcast_int(&isGS,   1, Root);
  sys_Bcast_int(&isTD,   1, Root);
  
  check(Mret > Nsite, "Invalid Mret");
  check(Nsweep > 0, "Invalid Nsweep");
  
  return;
  
 error:
  abort();
}

static void dmrg_Print(void)
{
  if(Rank == Root) {
    printf("%s\n", dmrglabel);
    printf("%s DMRG parameters:\n", dmrglabel);
    printf("%s\n", dmrglabel);
    printf("%s Nflavor = %5d\n", dmrglabel, Nflavor);
    printf("%s Nsite   = %5d\n", dmrglabel, Nsite);
    printf("%s Nblock  = %5d\n", dmrglabel, Nblock);
    printf("%s Mret    = %5d\n", dmrglabel, Mret);
    printf("%s Nsweep  = %5d\n", dmrglabel, Nsweep);
    printf("%s\n", dmrglabel);
  }
}
// DMRG initialization
void dmrg_Init(void)
{
  int iblock, block, blockmax, ndim, mdim, iflav;
  
  Nflavor   = Nspin;
  Nsite     = numeric_Ipow(NSP, Nflavor);
  Nblock    = Length - 3;
  NblockInf = (Length + Length%2) / 2 - 1;
  
  dmrg_Read();
  dmrg_Print();
  
  blockmax = log(Mret) / log(Nsite);
  
  Mretain = (int *) calloc(Nblock, sizeof(int));
  
  Hblock  = (mat2d_dreal **)  calloc(Nblock, sizeof(mat2d_dreal *));
  Wblock  = (mat2d_dreal **)  calloc(Nblock, sizeof(mat2d_dreal *));
  SZblock = (mat2d_dreal ***) calloc(Nblock, sizeof(mat2d_dreal **));
  SPblock = (mat2d_dreal ***) calloc(Nblock, sizeof(mat2d_dreal **));
  for(iblock = 0; iblock < Nblock; ++iblock) {
    block = iblock+1;
    ndim = (block > blockmax) ? Mret : numeric_Ipow(Nsite, block);
    Mretain[iblock] = ndim;
    Hblock[iblock] = mat2d_New_dreal();
    check_mem(Hblock[iblock], "Hblock[iblock]");
    Hblock[iblock]->Alloc(ndim, ndim, Hblock[iblock]);
    
    block = iblock+2;
    mdim = (block > blockmax) ? Mret : numeric_Ipow(Nsite, block);
    Wblock[iblock] = mat2d_New_dreal();
    check_mem(Wblock[iblock], "Wblock[iblock]");
    Wblock[iblock]->Alloc(ndim*Nsite, mdim, Wblock[iblock]);
    
    SZblock[iblock] = (mat2d_dreal **) calloc(Nflavor, sizeof(mat2d_dreal *));
    SPblock[iblock] = (mat2d_dreal **) calloc(Nflavor, sizeof(mat2d_dreal *));
    for(iflav = 0; iflav < Nflavor; ++iflav) {
      SZblock[iblock][iflav] = mat2d_New_dreal();
      SZblock[iblock][iflav]->Alloc(ndim, ndim, SZblock[iblock][iflav]);
      SPblock[iblock][iflav] = mat2d_New_dreal();
      SPblock[iblock][iflav]->Alloc(ndim, ndim, SPblock[iblock][iflav]);
    }
  }
  
  HblockL = mat2d_New_dreal(); HblockR = mat2d_New_dreal();
  HsiteL  = mat2d_New_dreal(); HsiteR  = mat2d_New_dreal();
  
  VblocksiteL = mat2d_New_dreal();
  VblocksiteR = mat2d_New_dreal();
  Vsitesite   = mat2d_New_dreal();
  
  SZblockL = (mat2d_dreal **) calloc(Nflavor, sizeof(mat2d_dreal *));
  SZblockR = (mat2d_dreal **) calloc(Nflavor, sizeof(mat2d_dreal *));
  SZsiteL  = (mat2d_dreal **) calloc(Nflavor, sizeof(mat2d_dreal *));
  SZsiteR  = (mat2d_dreal **) calloc(Nflavor, sizeof(mat2d_dreal *));
  SPblockL = (mat2d_dreal **) calloc(Nflavor, sizeof(mat2d_dreal *));
  SPblockR = (mat2d_dreal **) calloc(Nflavor, sizeof(mat2d_dreal *));
  SPsiteL  = (mat2d_dreal **) calloc(Nflavor, sizeof(mat2d_dreal *));
  SPsiteR  = (mat2d_dreal **) calloc(Nflavor, sizeof(mat2d_dreal *));
  for(iflav = 0; iflav < Nflavor; ++iflav) {
    SZblockL[iflav] = mat2d_New_dreal(); SZblockR[iflav] = mat2d_New_dreal();
    SZsiteL[iflav]  = mat2d_New_dreal(); SZsiteR[iflav]  = mat2d_New_dreal();
    SPblockL[iflav] = mat2d_New_dreal(); SPblockR[iflav] = mat2d_New_dreal();
    SPsiteL[iflav]  = mat2d_New_dreal(); SPsiteR[iflav]  = mat2d_New_dreal();
  }
  
  DenGS = vect_New_dreal();
  DenGS->Alloc(Length, DenGS);
  
  EvcGS = vect_New_dreal();
  
  Dmat = mat2d_New_dreal();
  
  return;
  
 error:
  abort();
}
// DMRG exit
void dmrg_Quit(void)
{
  int iblock, iflav;
  
  mat2d_Del_dreal(Dmat);
  vect_Del_dreal(EvcGS);
  vect_Del_dreal(DenGS);
  
  if(SPsiteR) {
    for(iflav = Nflavor-1; iflav > -1; --iflav)
      mat2d_Del_dreal(SPsiteR[iflav]);
  } freeup(SPsiteR);
  if(SPsiteL) {
    for(iflav = Nflavor-1; iflav > -1; --iflav)
      mat2d_Del_dreal(SPsiteL[iflav]);
  } freeup(SPsiteL);
  if(SPblockR) {
    for(iflav = Nflavor-1; iflav > -1; --iflav)
      mat2d_Del_dreal(SPblockR[iflav]);
  } freeup(SPblockR);
  if(SPblockL) {
    for(iflav = Nflavor-1; iflav > -1; --iflav)
      mat2d_Del_dreal(SPblockL[iflav]);
  } freeup(SPblockL);
  
  if(SZsiteR) {
    for(iflav = Nflavor-1; iflav > -1; --iflav)
      mat2d_Del_dreal(SZsiteR[iflav]);
  } freeup(SZsiteR);
  if(SZsiteL) {
    for(iflav = Nflavor-1; iflav > -1; --iflav)
      mat2d_Del_dreal(SZsiteL[iflav]);
  } freeup(SZsiteL);
  if(SZblockR) {
    for(iflav = Nflavor-1; iflav > -1; --iflav)
      mat2d_Del_dreal(SZblockR[iflav]);
  } freeup(SZblockR);
  if(SZblockL) {
    for(iflav = Nflavor-1; iflav > -1; --iflav)
      mat2d_Del_dreal(SZblockL[iflav]);
  } freeup(SZblockL);
  
  mat2d_Del_dreal(Vsitesite);
  mat2d_Del_dreal(VblocksiteR);
  mat2d_Del_dreal(VblocksiteL);
  
  mat2d_Del_dreal(HsiteR);  mat2d_Del_dreal(HsiteL);
  mat2d_Del_dreal(HblockR); mat2d_Del_dreal(HblockL);
  
  if(SPblock) {
    for(iblock = Nblock-1; iblock > -1; --iblock) {
      for(iflav = Nflavor-1; iflav > -1; --iflav)
        mat2d_Del_dreal(SPblock[iblock][iflav]);
      freeup(SPblock[iblock]);
    }
    freeup(SPblock);
  }
  if(SZblock) {
    for(iblock = Nblock-1; iblock > -1; --iblock) {
      for(iflav = Nflavor-1; iflav > -1; --iflav)
        mat2d_Del_dreal(SZblock[iblock][iflav]);
      freeup(SZblock[iblock]);
    }
    freeup(SZblock);
  }
  if(Wblock) {
    for(iblock = Nblock-1; iblock > -1; --iblock)
      mat2d_Del_dreal(Wblock[iblock]);
    freeup(Wblock);
  }
  if(Hblock) {
    for(iblock = Nblock-1; iblock > -1; --iblock)
      mat2d_Del_dreal(Hblock[iblock]);
    freeup(Hblock);
  }
  
  freeup(Mretain);
  freeup(EigD);
}
// Index in a superblock
inline int dmrg_indxBB(int ndimBL, int ndimSL, int ndimBR, int ndimSR,
                       int indxBL, int indxSL, int indxBR, int indxSR)
{
  int indxBB;
  
  indxBB = indxSR * ndimBL * ndimSL * ndimBR;
  indxBB+= indxBR * ndimBL * ndimSL;
  indxBB+= indxSL * ndimBL;
  indxBB+= indxBL;
  
  return indxBB;
}
// Pauli matrix SigmaZ
void dmrg_SigmaZ(mat2d_dreal *Sigma)
{
  Sigma->Reset(Sigma);
  Sigma->ptr[0][0] = 1.0;
  Sigma->ptr[1][1] =-1.0;
}
// Pauli matrix SigmaPlus
void dmrg_SigmaP(mat2d_dreal *Sigma)
{
  Sigma->Reset(Sigma);
  Sigma->ptr[1][0] = 1.0;
}
// Sigma operator on a site
void dmrg_Ssite(mat2d_dreal **Ssite, Paulimat dmrg_Sigma)
{
  int         isp_1, isp_2, jsp_1, jsp_2, ifock, jfock;
  mat2d_dreal *Sigma = NULL;
  
  Sigma = mat2d_New_dreal();
  Sigma->Alloc(NSP, NSP, Sigma);
  
  dmrg_Sigma(Sigma);
  
  switch(Nflavor) {
  case 1:
    mat2d_Copy_dreal(Sigma, Ssite[0]);
    break;
  case 2:
    for(isp_2 = 0; isp_2 < NSP; ++isp_2) {
      for(isp_1 = 0; isp_1 < NSP; ++isp_1) {
        ifock = isp_2*NSP + isp_1;
        for(jsp_1 = 0; jsp_1 < NSP; ++jsp_1) {
          jfock = isp_2*NSP + jsp_1;
          Ssite[0]->ptr[ifock][jfock] = Sigma->ptr[isp_1][jsp_1];
        }
      }
    }
    for(isp_1 = 0; isp_1 < NSP; ++isp_1) {
      for(isp_2 = 0; isp_2 < NSP; ++isp_2) {
        ifock = isp_2*NSP + isp_1;
        for(jsp_2 = 0; jsp_2 < NSP; ++jsp_2) {
          jfock = jsp_2*NSP + isp_1;
          Ssite[1]->ptr[ifock][jfock] = Sigma->ptr[isp_2][jsp_2];
        }
      }
    }
    break;
  default:
    abort();
  }
  
  mat2d_Del_dreal(Sigma);
}
// Occupation matrix on a site or block, with known SZmat
void dmrg_Occmat(mat2d_dreal **SZmat, mat2d_dreal **Occmat)
{
  int iflav, ndim, idim, jdim;
  
  ndim = SZmat[0]->nrow;
  
  for(iflav = 0; iflav < Nflavor; ++iflav) {
    Occmat[iflav]->Reset(Occmat[iflav]);
    for(idim = 0; idim < ndim; ++idim) {
      for(jdim = 0; jdim < ndim; ++jdim) {
        Occmat[iflav]->ptr[idim][jdim] += 0.5 * SZmat[iflav]->ptr[idim][jdim];
      }
      Occmat[iflav]->ptr[idim][idim] += 0.5;
    }
  }
}
// Hamiltonian matrix on a site or block, with known SimgaZ matrix
void dmrg_Hmat(dreal EE, dreal UU, mat2d_dreal **SZmat, mat2d_dreal *Hmat)
{
  int         ndim, idim, jdim, iflav;
  mat2d_dreal **Occmat = NULL;
  
  ndim = Hmat->nrow;
  
  Occmat = (mat2d_dreal **) calloc(Nflavor, sizeof(mat2d_dreal *));
  for(iflav = 0; iflav < Nflavor; ++iflav) {
    Occmat[iflav] = mat2d_New_dreal();
    Occmat[iflav]->Alloc(ndim, ndim, Occmat[iflav]);
  }
  
  dmrg_Occmat(SZmat, Occmat);
  
  Hmat->Reset(Hmat);
  
  for(iflav = 0; iflav < Nflavor; ++iflav) {
    for(idim = 0; idim < ndim; ++idim) {
      for(jdim = 0; jdim < ndim; ++jdim) {
        Hmat->ptr[idim][jdim] += EE * Occmat[iflav]->ptr[idim][jdim];
      }
    }
  }
  
  if(Nflavor == 2)
    lapack_dgemm(ndim, ndim, ndim, 'N', 'N', UU, Occmat[0]->addr, Occmat[1]->addr, 1.0, Hmat->addr);
  
  for(iflav = Nflavor-1; iflav > -1; --iflav)
    mat2d_Del_dreal(Occmat[iflav]);
  freeup(Occmat);
}
// Tunneling matrix between two sites or blocks
void dmrg_Vmat(dreal tt, int ndim1, int ndim2,
               mat2d_dreal **SPmat1, mat2d_dreal **SPmat2,
               mat2d_dreal *Vmat)
{
  int   iflav, idim1, jdim1, idim2, jdim2, idim, jdim;
  dreal xp1, xm1, xp2, xm2;
  
  Vmat->Reset(Vmat);
  
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
            Vmat->ptr[idim][jdim] += tt*(xp1*xm2 + xm1*xp2);
          }
        }
      }
    }
  }
}
// Superblock initialization
static void dmrg_Superblock(int indxBL, int indxSL, int indxBR, int indxSR)
{
  int   iflav;
  dreal EE, UU, tL, tR, tC;
  
  NdimBL = Hblock[indxBL]->nrow;
  NdimBR = Hblock[indxBR]->nrow;
  NdimSL = Nsite;
  NdimSR = Nsite;
  NdimBB = NdimBL * NdimSL * NdimBR * NdimSR;
  NdimSS = NdimSL * NdimSR;
  NdimL  = NdimBL * NdimSL;
  NdimR  = NdimBR * NdimSR;
  
  if(HblockL) HblockL->Deall(HblockL);
  if(HblockR) HblockR->Deall(HblockR);
  if(HsiteL)  HsiteL->Deall(HsiteL);
  if(HsiteR)  HsiteR->Deall(HsiteR);
  
  HblockL->Alloc(NdimBL, NdimBL, HblockL);
  HblockR->Alloc(NdimBR, NdimBR, HblockR);
  HsiteL->Alloc(NdimSL, NdimSL, HsiteL);
  HsiteR->Alloc(NdimSR, NdimSR, HsiteR);
  
  for(iflav = 0; iflav < Nflavor; ++iflav) {
    if(SZblockL[iflav]) SZblockL[iflav]->Deall(SZblockL[iflav]);
    if(SZblockR[iflav]) SZblockR[iflav]->Deall(SZblockR[iflav]);
    if(SZsiteL[iflav])  SZsiteL[iflav]->Deall(SZsiteL[iflav]);
    if(SZsiteR[iflav])  SZsiteR[iflav]->Deall(SZsiteR[iflav]);
    if(SPblockL[iflav]) SPblockL[iflav]->Deall(SPblockL[iflav]);
    if(SPblockR[iflav]) SPblockR[iflav]->Deall(SPblockR[iflav]);
    if(SPsiteL[iflav])  SPsiteL[iflav]->Deall(SPsiteL[iflav]);
    if(SPsiteR[iflav])  SPsiteR[iflav]->Deall(SPsiteR[iflav]);
    
    SZblockL[iflav]->Alloc(NdimBL, NdimBL, SZblockL[iflav]);
    SZblockR[iflav]->Alloc(NdimBR, NdimBR, SZblockR[iflav]);
    SZsiteL[iflav]->Alloc(NdimSL, NdimSL, SZsiteL[iflav]);
    SZsiteR[iflav]->Alloc(NdimSR, NdimSR, SZsiteR[iflav]);
    SPblockL[iflav]->Alloc(NdimBL, NdimBL, SPblockL[iflav]);
    SPblockR[iflav]->Alloc(NdimBR, NdimBR, SPblockR[iflav]);
    SPsiteL[iflav]->Alloc(NdimSL, NdimSL, SPsiteL[iflav]);
    SPsiteR[iflav]->Alloc(NdimSR, NdimSR, SPsiteR[iflav]);
  }
  
  mat2d_Copy_dreal(Hblock[indxBL], HblockL);
  mat2d_Copy_dreal(Hblock[indxBR], HblockR);
  
  for(iflav = 0; iflav < Nflavor; ++iflav) {
    mat2d_Copy_dreal(SZblock[indxBL][iflav], SZblockL[iflav]);
    mat2d_Copy_dreal(SZblock[indxBR][iflav], SZblockR[iflav]);
    mat2d_Copy_dreal(SPblock[indxBL][iflav], SPblockL[iflav]);
    mat2d_Copy_dreal(SPblock[indxBR][iflav], SPblockR[iflav]);
  }
  
  dmrg_Ssite(SZsiteL, dmrg_SigmaZ); dmrg_Ssite(SZsiteR, dmrg_SigmaZ);
  dmrg_Ssite(SPsiteL, dmrg_SigmaP); dmrg_Ssite(SPsiteR, dmrg_SigmaP);
  
  EE = Echain[indxSL]; UU = Uchain[indxSL];
  dmrg_Hmat(EE, UU, SZsiteL, HsiteL);
  
  EE = Echain[indxSR]; UU = Uchain[indxSR];
  dmrg_Hmat(EE, UU, SZsiteR, HsiteR);
  
  tL = Tchain[indxBL];
  tR = Tchain[indxBR];
  tC = Tchain[indxSL];
  
  if(VblocksiteL) VblocksiteL->Deall(VblocksiteL);
  if(VblocksiteR) VblocksiteR->Deall(VblocksiteR);
  if(Vsitesite)   Vsitesite->Deall(Vsitesite);
  
  VblocksiteL->Alloc(NdimL, NdimL, VblocksiteL);
  VblocksiteR->Alloc(NdimR, NdimR, VblocksiteR);
  Vsitesite->Alloc(NdimSS, NdimSS, Vsitesite);
  
  dmrg_Vmat(tL, NdimBL, NdimSL, SPblockL, SPsiteL, VblocksiteL);
  dmrg_Vmat(tR, NdimBR, NdimSR, SPblockR, SPsiteR, VblocksiteR);
  dmrg_Vmat(tC, NdimSL, NdimSR, SPsiteL,  SPsiteR, Vsitesite);
  
  if(EvcGS->ndim != NdimBB) {
    EvcGS->Deall(EvcGS);
    EvcGS->Alloc(NdimBB, EvcGS);
  }
  
  Dmat->Deall(Dmat);
  Dmat->Alloc(NdimL, NdimL, Dmat);
}
// Ground state wavefunction intialization for infinite DMRG calculations
static void dmrg_EvcGS_Init_Inf(void)
{
  int        idimBL, idimSL, idimL, idimBR, idimSR, idimR, idimBB;
  vect_dreal *vect_rand = NULL;
  
  check(NdimBL == NdimBR, "Inequivalent NdimBL and NdimBR");
  check(NdimSL == NdimSR, "Inequivalent NdimSL and NdimSR");
  
  vect_rand = vect_New_dreal();
  check_mem(vect_rand, "vect_rand");
  vect_rand->Alloc(NdimL, vect_rand);
  
  vect_RandNorm_dreal(vect_rand);
  
  for(idimSL = 0; idimSL < NdimSL; ++idimSL) {
    for(idimBL = 0; idimBL < NdimBL; ++idimBL) {
      idimL = idimSL*NdimBL + idimBL;
      for(idimSR = 0; idimSR < NdimSR; ++idimSR) {
        for(idimBR = 0; idimBR < NdimBR; ++idimBR) {
          idimR  = idimSR*NdimBR + idimBR;
          idimBB = dmrg_indxBB(NdimBL, NdimSL, NdimBR, NdimSR, idimBL, idimSL, idimBR, idimSR);
          EvcGS->addr[idimBB] = vect_rand->addr[idimL] * vect_rand->addr[idimR];
        }
      }
    }
  }
  
  vect_Del_dreal(vect_rand);
  
  return;
  
 error:
  abort();
}
// Ground state wavefunction initialization for remnant calculations
static void dmrg_EvcGS_Init_Rem(void)
{
  int        idimBL, idimSL, idimL, idimBR, idimSR, idimR, idimBB;
  vect_dreal *vectL_rand = NULL, *vectR_rand = NULL;
  
  vectL_rand = vect_New_dreal();
  check_mem(vectL_rand, "vectL_rand");
  vectL_rand->Alloc(NdimL, vectL_rand);
  
  vectR_rand = vect_New_dreal();
  check_mem(vectR_rand, "vectR_rand");
  vectR_rand->Alloc(NdimR, vectR_rand);
  
  vect_RandNorm_dreal(vectL_rand);
  vect_RandNorm_dreal(vectR_rand);
  
  for(idimSL = 0; idimSL < NdimSL; ++idimSL) {
    for(idimBL = 0; idimBL < NdimBL; ++idimBL) {
      idimL = idimSL*NdimBL + idimBL;
      for(idimSR = 0; idimSR < NdimSR; ++idimSR) {
        for(idimBR = 0; idimBR < NdimBR; ++idimBR) {
          idimR  = idimSR*NdimBR + idimBR;
          idimBB = dmrg_indxBB(NdimBL, NdimSL, NdimBR, NdimSR, idimBL, idimSL, idimBR, idimSR);
          EvcGS->addr[idimBB] = vectL_rand->addr[idimL] * vectR_rand->addr[idimR];
        }
      }
    }
  }
  
  vect_Del_dreal(vectR_rand);
  vect_Del_dreal(vectL_rand);
  
  return;
  
 error:
  abort();
}
// Transform from a mirror to a normal wavefunction
static void dmrg_Evc_Mirror(int ndimBL, int ndimSL, int ndimBR, int ndimSR,
                            const vect_dreal *Evc_R2L, vect_dreal *Evc_L2R)
{
  int idimBL, idimSL, idimBR, idimSR, idimBB_R2L, idimBB_L2R;
  
  for(idimSR = 0; idimSR < ndimSR; ++idimSR) {
    for(idimBR = 0; idimBR < ndimBR; ++idimBR) {
      for(idimSL = 0; idimSL < ndimSL; ++idimSL) {
        for(idimBL = 0; idimBL < ndimBL; ++idimBL) {
          idimBB_R2L = dmrg_indxBB(ndimBR, ndimSR, ndimBL, ndimSL, idimBR, idimSR, idimBL, idimSL);
          idimBB_L2R = dmrg_indxBB(ndimBL, ndimSL, ndimBR, ndimSR, idimBL, idimSL, idimBR, idimSR);
          Evc_L2R->addr[idimBB_L2R] = Evc_R2L->addr[idimBB_R2L];
        }
      }
    }
  }
}
// Ground state wavefunction prediction in finite DMRG calculations
static void dmrg_EvcGS_Predict(int indxBL, int indxBR, const vect_dreal *EvcGS_last, vect_dreal *EvcGS_next)
{
  int        ndimBL_next, ndimSL_next, ndimBR_next, ndimSR_next, ndimL_next, ndimR_next,
             ndimBL_last, ndimSL_last, ndimBR_last, ndimSR_last, ndimL_last, ndimR_last,
             ndimBB_contr;
  int        idimBL_next, idimSL_next, idimBR_next, idimSR_next, idimR_next, idimBB_next,
             idimBR_last, idimSR_last, idimR_last, idimBB_contr;
  dreal      *Oprptr = NULL, *vectptr = NULL, *Oprvectptr = NULL;
  vect_dreal *vectBR_last = NULL, *vectR_next = NULL, *vect_contr = NULL;
  
  check(indxBL > 0, "Invalid indxBL");
  check(indxBR < Nblock-1, "Invalid indxBL");
  
  ndimBL_next = Hblock[indxBL]->nrow;
  ndimBR_next = Hblock[indxBR]->nrow;
  ndimSL_next = Nsite;
  ndimSR_next = Nsite;
  ndimL_next  = ndimBL_next*ndimSL_next;
  ndimR_next  = ndimBR_next*ndimSR_next;
  
  ndimBL_last = Hblock[indxBL-1]->nrow;
  ndimBR_last = Hblock[indxBR+1]->nrow;
  ndimSL_last = Nsite;
  ndimSR_last = Nsite;
  ndimL_last  = ndimBL_last*ndimSL_last;
  ndimR_last  = ndimBR_last*ndimSR_last;
  
  ndimBB_contr= ndimBL_next*ndimBR_last*ndimSR_last;
  
  vectBR_last = vect_New_dreal();
  check_mem(vectBR_last, "vectBR_last");
  vectBR_last->Alloc(ndimBR_last, vectBR_last);
  
  vectR_next = vect_New_dreal();
  check_mem(vectR_next, "vectR_next");
  vectR_next->Alloc(ndimR_next, vectR_next);
  
  vect_contr = vect_New_dreal();
  check_mem(vect_contr, "vect_contr");
  vect_contr->Alloc(ndimBB_contr, vect_contr);
  
  Oprptr = Wblock[indxBL-1]->addr;
  for(idimSR_last = 0; idimSR_last < ndimSR_last; ++idimSR_last) {
    for(idimBR_last = 0; idimBR_last < ndimBR_last; ++idimBR_last) {
      idimR_last = idimSR_last*ndimBR_last + idimBR_last;
      vectptr    = EvcGS_last->addr + idimR_last*ndimL_last;
      Oprvectptr = vect_contr->addr + idimR_last*ndimBL_next;
      lapack_dgemm(ndimBL_next, 1, ndimL_last, 'T', 'N', 1.0, Oprptr, vectptr, 0.0, Oprvectptr);
    }
  }
  
  Oprptr = Wblock[indxBR]->addr;
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
      lapack_dgemm(ndimR_next, 1, ndimBR_last, 'N', 'N', 1.0, Oprptr, vectptr, 0.0, Oprvectptr);
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
  
  vect_Del_dreal(vect_contr);
  vect_Del_dreal(vectR_next);
  vect_Del_dreal(vectBR_last);
  
  return;
  
 error:
  abort();
}
// Preconditioner for Davidson diagonalization
void dmrg_PBBvect(dreal rho, const vect_dreal *vect, vect_dreal *Pvect)
{
  const dreal offset = 1e-1;
  int         idimBL, idimBR, idimSL, idimSR, idimBB;
  dreal       diag;
  
  for(idimSR = 0; idimSR < NdimSR; ++idimSR) {
    for(idimBR = 0; idimBR < NdimBR; ++idimBR) {
      for(idimSL = 0; idimSL < NdimSL; ++idimSL) {
        for(idimBL = 0; idimBL < NdimBL; ++idimBL) {
          idimBB = dmrg_indxBB(NdimBL, NdimSL, NdimBR, NdimSR, idimBL, idimSL, idimBR, idimSR);
          diag =-HblockL->ptr[idimBL][idimBL] - HsiteL->ptr[idimSL][idimSL]
               - HblockR->ptr[idimBR][idimBR] - HsiteR->ptr[idimSR][idimSR];
          diag += rho;
          Pvect->addr[idimBB] = 1.0/(diag - offset) * vect->addr[idimBB];
        }
      }
    }
  }
}
// Block operator block-site vector multiplication, NO value added
static void dmrg_OprBvectBS(int ndimB, int ndimS, int mdimB, const mat2d_dreal *OprB,
                            const vect_dreal *vectBS, vect_dreal *OprBvectBS)
{
  int        idimS;
  vect_dreal *vectB = NULL, *OprBvectB = NULL;
  
  vectB = vect_New_dreal();
  check_mem(vectB, "vectB");
  vectB->Alloc(mdimB, vectB);
  
  OprBvectB = vect_New_dreal();
  check_mem(OprBvectB, "OprBvectB");
  OprBvectB->Alloc(ndimB, OprBvectB);
  
  for(idimS = 0; idimS < ndimS; ++idimS) {
    memcpy(vectB->addr, vectBS->addr+idimS*mdimB, sizeof(dreal)*mdimB);
    lapack_dgemm(ndimB, 1, mdimB, 'N', 'N', 1.0, OprB->addr, vectB->addr, 0.0, OprBvectB->addr);
    memcpy(OprBvectBS->addr+idimS*ndimB, OprBvectB->addr, sizeof(dreal)*ndimB);
  }
  
  vect_Del_dreal(OprBvectB);
  vect_Del_dreal(vectB); 
  
  return;
  
 error:
  abort();
}
// Site operator block-site vector multiplication, NO value added
static void dmrg_OprSvectBS(int ndimB, int ndimS, int mdimS, const mat2d_dreal *OprS,
                            const vect_dreal *vectBS, vect_dreal *OprSvectBS)
{
  int        idimS, idimB, idimBS;
  vect_dreal *vectS = NULL, *OprSvectS = NULL;
  
  vectS = vect_New_dreal();
  check_mem(vectS, "vectS");
  vectS->Alloc(mdimS, vectS);
  
  OprSvectS = vect_New_dreal();
  check_mem(OprSvectS, "OprSvectS");
  OprSvectS->Alloc(ndimS, OprSvectS);
  
  for(idimB = 0; idimB < ndimB; ++idimB) {
    for(idimS = 0; idimS < mdimS; ++idimS) {
      idimBS = idimS*ndimB + idimB;
      vectS->addr[idimS] = vectBS->addr[idimBS];
    }
    lapack_dgemm(ndimS, 1, mdimS, 'N', 'N', 1.0, OprS->addr, vectS->addr, 0.0, OprSvectS->addr);
    for(idimS = 0; idimS < ndimS; ++idimS) {
      idimBS = idimS*ndimB + idimB;
      OprSvectBS->addr[idimBS] = OprSvectS->addr[idimS];
    }
  }
  
  vect_Del_dreal(OprSvectS);
  vect_Del_dreal(vectS); 
  
  return;
  
 error:
  abort();
}
// Block-site operator block-site vector multiplication, NO value added
static void dmrg_OprBSvectBS(int ndimBS, int mdimBS, const mat2d_dreal *OprBS,
                             const vect_dreal *vectBS, vect_dreal *OprBSvectBS)
{
  lapack_dgemm(ndimBS, 1, mdimBS, 'N', 'N', 1.0, OprBS->addr, vectBS->addr, 0.0, OprBSvectBS->addr);
}
// BlockL operator superblock vector multiplication, value added
static void dmrg_OprBLvect(int ndimBL, int ndimSL, int ndimBR, int ndimSR, int mdimBL,
                           const mat2d_dreal *OprBL, const vect_dreal *vect, vect_dreal *OprBLvect)
{
  int        ndimL, mdimL, ndimBB, idimBR, idimSR, idimR, idimBB;
  vect_dreal *vectL = NULL, *OprBLvectL = NULL, *OprBLvect_save = NULL;
  
  ndimBB = ndimBL*ndimSL*ndimBR*ndimSL;
  ndimL  = ndimBL*ndimSL;
  mdimL  = mdimBL*ndimSL;
  
  vectL = vect_New_dreal();
  check_mem(vectL, "vectL");
  vectL->Alloc(mdimL, vectL);
  
  OprBLvectL = vect_New_dreal();
  check_mem(OprBLvectL, "OprBLvectL");
  OprBLvectL->Alloc(ndimL, OprBLvectL);
  
  OprBLvect_save = vect_New_dreal();
  check_mem(OprBLvect_save, "OprBLvect_save");
  OprBLvect_save->Alloc(ndimBB, OprBLvect_save);
  
  vect_Copy_dreal(OprBLvect, OprBLvect_save);
  
  for(idimSR = 0; idimSR < ndimSR; ++idimSR) {
    for(idimBR = 0; idimBR < ndimBR; ++idimBR) {
      idimR = idimSR*ndimBR + idimBR;
      memcpy(vectL->addr, vect->addr+idimR*mdimL, sizeof(dreal)*mdimL);
      dmrg_OprBvectBS(ndimBL, ndimSL, mdimBL, OprBL, vectL, OprBLvectL);
      memcpy(OprBLvect->addr+idimR*ndimL, OprBLvectL->addr, sizeof(dreal)*ndimL);
    }
  }
  
  for(idimBB = 0; idimBB < ndimBB; ++idimBB)
    OprBLvect->addr[idimBB] += OprBLvect_save->addr[idimBB];
  
  vect_Del_dreal(OprBLvect_save);
  vect_Del_dreal(OprBLvectL);
  vect_Del_dreal(vectL);
  
  return;
  
 error:
  abort();
}
// SiteL operator superblock vector multiplication, value added
static void dmrg_OprSLvect(int ndimBL, int ndimSL, int ndimBR, int ndimSR, int mdimSL,
                           const mat2d_dreal *OprSL, const vect_dreal *vect, vect_dreal *OprSLvect)
{
  int        ndimL, mdimL, ndimBB, idimBR, idimSR, idimR, idimBB;
  vect_dreal *vectL = NULL, *OprSLvectL = NULL,
             *OprSLvect_save = NULL;
  
  ndimBB = ndimBL*ndimSL*ndimBR*ndimSL;
  ndimL  = ndimBL*ndimSL;
  mdimL  = ndimBL*mdimSL;
  
  vectL = vect_New_dreal();
  check_mem(vectL, "vectL");
  vectL->Alloc(mdimL, vectL);
  
  OprSLvectL = vect_New_dreal();
  check_mem(OprSLvectL, "OprSLvectL");
  OprSLvectL->Alloc(ndimL, OprSLvectL);
  
  OprSLvect_save = vect_New_dreal();
  check_mem(OprSLvect_save, "OprSLvect_save");
  OprSLvect_save->Alloc(ndimBB, OprSLvect_save);
  
  vect_Copy_dreal(OprSLvect, OprSLvect_save);
  
  for(idimSR = 0; idimSR < ndimSR; ++idimSR) {
    for(idimBR = 0; idimBR < ndimBR; ++idimBR) {
      idimR = idimSR*ndimBR + idimBR;
      memcpy(vectL->addr, vect->addr+idimR*mdimL, sizeof(dreal)*mdimL);
      dmrg_OprSvectBS(ndimBL, ndimSL, mdimSL, OprSL, vectL, OprSLvectL);
      memcpy(OprSLvect->addr+idimR*ndimL, OprSLvectL->addr, sizeof(dreal)*ndimL);
    }
  }
  
  for(idimBB = 0; idimBB < ndimBB; ++idimBB)
    OprSLvect->addr[idimBB] += OprSLvect_save->addr[idimBB];
  
  vect_Del_dreal(OprSLvect_save);
  vect_Del_dreal(OprSLvectL);
  vect_Del_dreal(vectL);
  
  return;
  
 error:
  abort();
}
// BlockR operator superblock operator multiplication, value added
static void dmrg_OprBRvect(int ndimBL, int ndimSL, int ndimBR, int ndimSR, int mdimBR,
                           const mat2d_dreal *OprBR, const vect_dreal *vect, vect_dreal *OprBRvect)
{
  int        ndimR, mdimR, ndimBB, idimBL, idimSL, idimBR, idimSR, idimR, idimBB;
  vect_dreal *vectR = NULL, *OprBRvectR = NULL;
  
  ndimBB = ndimBL*ndimSL*ndimBR*ndimSL;
  ndimR  = ndimBR*ndimSR;
  mdimR  = mdimBR*ndimSR;
  
  vectR = vect_New_dreal();
  check_mem(vectR, "vectR");
  vectR->Alloc(mdimR, vectR);
  
  OprBRvectR = vect_New_dreal();
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
      dmrg_OprBvectBS(ndimBR, ndimSR, mdimBR, OprBR, vectR, OprBRvectR);
      for(idimSR = 0; idimSR < ndimSR; ++idimSR) {
        for(idimBR = 0; idimBR < ndimBR; ++idimBR) {
          idimR  = idimSR*ndimBR + idimBR;
          idimBB = dmrg_indxBB(ndimBL, ndimSL, ndimBR, ndimSR, idimBL, idimSL, idimBR, idimSR);
          OprBRvect->addr[idimBB] += OprBRvectR->addr[idimR];
        }
      }
    }
  }
  
  vect_Del_dreal(OprBRvectR);
  vect_Del_dreal(vectR);
  
  return;
  
 error:
  abort();
}
// SiteR operator superblock vector multiplication, value added
static void dmrg_OprSRvect(int ndimBL, int ndimSL, int ndimBR, int ndimSR, int mdimSR,
                           const mat2d_dreal *OprSR, const vect_dreal *vect, vect_dreal *OprSRvect)
{
  int        ndimR, mdimR, ndimBB, idimBL, idimSL, idimBR, idimSR, idimR, idimBB;
  vect_dreal *vectR = NULL, *OprSRvectR = NULL;
  
  ndimBB = ndimBL*ndimSL*ndimBR*ndimSL;
  ndimR  = ndimBR*ndimSR;
  mdimR  = ndimBR*mdimSR;
  
  vectR = vect_New_dreal();
  check_mem(vectR, "vectR");
  vectR->Alloc(mdimR, vectR);
  
  OprSRvectR = vect_New_dreal();
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
      dmrg_OprSvectBS(ndimBR, ndimSR, mdimSR, OprSR, vectR, OprSRvectR);
      for(idimSR = 0; idimSR < ndimSR; ++idimSR) {
        for(idimBR = 0; idimBR < ndimBR; ++idimBR) {
          idimR  = idimSR*ndimBR + idimBR;
          idimBB = dmrg_indxBB(ndimBL, ndimSL, ndimBR, ndimSR, idimBL, idimSL, idimBR, idimSR);
          OprSRvect->addr[idimBB] += OprSRvectR->addr[idimR];
        }
      }
    }
  }
  
  vect_Del_dreal(OprSRvectR);
  vect_Del_dreal(vectR);
  
  return;
  
 error:
  abort();
}

static void dmrg_OprLvect(int ndimBL, int ndimSL, int ndimBR, int ndimSR, int mdimBL, int mdimSL,
                          const mat2d_dreal *OprL, const vect_dreal *vect, vect_dreal *OprLvect)
{
  int        ndimBB, ndimL, mdimL, idimBR, idimSR, idimBB, idimR;
  vect_dreal *vectL = NULL, *OprLvectL = NULL, *OprLvect_save = NULL;
  
  ndimBB = ndimBL*ndimSL*ndimBR*ndimSR;
  ndimL  = ndimBL*ndimSL;
  mdimL  = mdimBL*mdimSL;
  
  vectL = vect_New_dreal();
  check_mem(vectL, "vectL");
  vectL->Alloc(mdimL, vectL);
  
  OprLvectL = vect_New_dreal();
  check_mem(OprLvectL, "OprLvectL");
  OprLvectL->Alloc(ndimL, OprLvectL);
  
  OprLvect_save = vect_New_dreal();
  check_mem(OprLvect_save, "OprLvect_save");
  OprLvect_save->Alloc(ndimBB, OprLvect_save);
  
  vect_Copy_dreal(OprLvect, OprLvect_save);
  
  for(idimSR = 0; idimSR < ndimSR; ++idimSR) {
    for(idimBR = 0; idimBR < ndimBR; ++idimBR) {
      idimR = idimSR*ndimBR + idimBR;
      memcpy(vectL->addr, vect->addr+idimR*mdimL, sizeof(dreal)*mdimL);
      dmrg_OprBSvectBS(ndimL, mdimL, OprL, vectL, OprLvectL);
      memcpy(OprLvect->addr+idimR*ndimL, OprLvectL->addr, sizeof(dreal)*ndimL);
    }
  }
  
  for(idimBB = 0; idimBB < ndimBB; ++idimBB)
    OprLvect->addr[idimBB] += OprLvect_save->addr[idimBB];
  
  vect_Del_dreal(OprLvect_save);
  vect_Del_dreal(OprLvectL);
  vect_Del_dreal(vectL);
  
  return;
  
 error:
  abort();
}

static void dmrg_OprRvect(int ndimBL, int ndimSL, int ndimBR, int ndimSR, int mdimBR, int mdimSR,
                          const mat2d_dreal *OprR, const vect_dreal *vect, vect_dreal *OprRvect)
{
  int        ndimR, mdimR, ndimBB, idimBL, idimSL, idimBR, idimSR, idimR, idimBB;
  vect_dreal *vectR = NULL, *OprRvectR = NULL;
  
  ndimBB = ndimBL*ndimSL*ndimBR*ndimSL;
  ndimR  = ndimBR*ndimSR;
  mdimR  = mdimBR*mdimSR;
  
  vectR = vect_New_dreal();
  check_mem(vectR, "vectR");
  vectR->Alloc(mdimR, vectR);
  
  OprRvectR = vect_New_dreal();
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
      dmrg_OprBSvectBS(ndimR, mdimR, OprR, vectR, OprRvectR);
      for(idimSR = 0; idimSR < ndimSR; ++idimSR) {
        for(idimBR = 0; idimBR < ndimBR; ++idimBR) {
          idimR  = idimSR*ndimBR + idimBR;
          idimBB = dmrg_indxBB(ndimBL, ndimSL, ndimBR, ndimSR, idimBL, idimSL, idimBR, idimSR);
          OprRvect->addr[idimBB] += OprRvectR->addr[idimR];
        }
      }
    }
  }
  
  vect_Del_dreal(OprRvectR);
  vect_Del_dreal(vectR);
  
  return;
  
 error:
  abort();
}
// Site-site operator superblock vector multiplication, value added
static void dmrg_OprSSvect(int ndimBL, int ndimSL, int ndimBR, int ndimSR, int mdimSL, int mdimSR,
                           const mat2d_dreal *OprSS, const vect_dreal *vect, vect_dreal *OprSSvect)
{
  int        ndimSS, mdimSS, idimBL, idimBR, idimSL, idimSR, idimSS, idimBB;
  vect_dreal *vectSS = NULL, *OprSSvectSS = NULL;
  
  ndimSS = ndimSL*ndimSR;
  mdimSS = mdimSL*mdimSR;
  
  vectSS = vect_New_dreal();
  check_mem(vectSS, "vectSS");
  vectSS->Alloc(mdimSS, vectSS);
  
  OprSSvectSS = vect_New_dreal();
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
      lapack_dgemm(ndimSS, 1, mdimSS, 'N', 'N', 1.0, OprSS->addr, vectSS->addr, 0.0, OprSSvectSS->addr); 
      for(idimSR = 0; idimSR < ndimSR; ++idimSR) {
        for(idimSL = 0; idimSL < ndimSL; ++idimSL) {
          idimBB = dmrg_indxBB(ndimBL, ndimSL, ndimBR, ndimSR, idimBL, idimSL, idimBR, idimSR);
          idimSS = idimSR*ndimSL + idimSL;
          OprSSvect->addr[idimBB] += OprSSvectSS->addr[idimSS];
        }
      }
    }
  }
  
  vect_Del_dreal(OprSSvectSS);
  vect_Del_dreal(vectSS);
  
  return;
  
 error:
  abort();
}
// Superblock Hamiltonian superblock vector multiplication
void dmrg_HBBvect(const vect_dreal *vect, vect_dreal *HBBvect)
{
  HBBvect->Reset(HBBvect);
  
  dmrg_OprBLvect(NdimBL, NdimSL, NdimBR, NdimSR, NdimBL, HblockL, vect, HBBvect);
  dmrg_OprSLvect(NdimBL, NdimSL, NdimBR, NdimSR, NdimSL, HsiteL,  vect, HBBvect);
  dmrg_OprBRvect(NdimBL, NdimSL, NdimBR, NdimSR, NdimBR, HblockR, vect, HBBvect);
  dmrg_OprSRvect(NdimBL, NdimSL, NdimBR, NdimSR, NdimSR, HsiteR,  vect, HBBvect);
  
  dmrg_OprLvect(NdimBL, NdimSL, NdimBR, NdimSR, NdimBL, NdimSL, VblocksiteL, vect, HBBvect);
  dmrg_OprRvect(NdimBL, NdimSL, NdimBR, NdimSR, NdimBR, NdimSR, VblocksiteR, vect, HBBvect);
  
  dmrg_OprSSvect(NdimBL, NdimSL, NdimBR, NdimSR, NdimSL, NdimSR, Vsitesite, vect, HBBvect);
}
// Density matrix update
static void dmrg_Dmat(void)
{
  int idimBL, jdimBL, idimSL, jdimSL, idimL, jdimL, idimBB, jdimBB, idimSR, idimBR;
  
  Dmat->Reset(Dmat);
  
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
              Dmat->ptr[idimL][jdimL] += EvcGS->addr[jdimBB] * EvcGS->addr[idimBB];
            }
          }
        }
      }
    }
  }
}

static void dmrg_Dmat_fixed(mat2d_dreal *Wvc)
{
  int ndim, idim, idimL;
  
  Dmat->Reset(Dmat);
  
  ndim = Wvc->ncol;
  
  for(idim = 0; idim < ndim; ++idim) {
    idimL = NdimL-idim-1;
    memcpy(Dmat->addr+idimL*NdimL, Wvc->addr+idim*NdimL, sizeof(dreal)*NdimL);
  }
}

static void dmrg_Matsite(dreal EE, dreal UU, mat2d_dreal **SZsitemat, mat2d_dreal **SPsitemat, mat2d_dreal *Hsitemat)
{
  dmrg_Ssite(SZsitemat, dmrg_SigmaZ);
  dmrg_Ssite(SPsitemat, dmrg_SigmaP);
  dmrg_Hmat(EE, UU, SZsitemat, Hsitemat);
}
                         
static void dmrg_Matblock(mat2d_dreal **SZblockmat, mat2d_dreal **SPblockmat, mat2d_dreal *Hblockmat)
{
  int        Ndim, idim, jdim, idimL, jdimL, iflav;
  dreal      vOprv;
  vect_dreal *vectLi = NULL, *vectLj = NULL,
             *OprBLvectL = NULL, *OprSLvectL = NULL, *OprLvectL = NULL;
  
  Ndim = Hblockmat->nrow;
  
  vectLi = vect_New_dreal();
  check_mem(vectLi, "vectLi");
  vectLi->Alloc(NdimL, vectLi);
  
  vectLj = vect_New_dreal();
  check_mem(vectLj, "vectLj");
  vectLj->Alloc(NdimL, vectLj);
  
  OprBLvectL = vect_New_dreal();
  check_mem(OprBLvectL, "OprBLvectL");
  OprBLvectL->Alloc(NdimL, OprBLvectL);
  
  OprSLvectL = vect_New_dreal();
  check_mem(OprSLvectL, "OprSLvectL");
  OprSLvectL->Alloc(NdimL, OprSLvectL);
  
  OprLvectL = vect_New_dreal();
  check_mem(OprLvectL, "OprLvectL");
  OprLvectL->Alloc(NdimL, OprLvectL);
  // Update Hblock, note that Dmat now contains basis
  for(idim = 0; idim < Ndim; ++idim) {
    idimL = NdimL-idim-1;
    memcpy(vectLi->addr, Dmat->ptr[idimL], NdimL*sizeof(dreal));
    dmrg_OprBvectBS(NdimBL, NdimSL, NdimBL, HblockL, vectLi, OprBLvectL);
    dmrg_OprSvectBS(NdimBL, NdimSL, NdimSL, HsiteL, vectLi, OprSLvectL);
    dmrg_OprBSvectBS(NdimL, NdimL, VblocksiteL, vectLi, OprLvectL);
    for(idimL = 0; idimL < NdimL; ++idimL)
      OprLvectL->addr[idimL] += OprBLvectL->addr[idimL] + OprSLvectL->addr[idimL];
    for(jdim = 0; jdim < Ndim; ++jdim) {
      jdimL = NdimL-jdim-1;
      memcpy(vectLj->addr, Dmat->ptr[jdimL], NdimL*sizeof(dreal));
      vOprv = vect_Dotprod_dreal(vectLj, OprLvectL);
      Hblockmat->ptr[idim][jdim] = vOprv;
    }
  }
  
  for(iflav = 0; iflav < Nflavor; ++iflav) {
    // Update SZblock
    for(idim = 0; idim < Ndim; ++idim) {
      idimL = NdimL-idim-1;
      memcpy(vectLi->addr, Dmat->ptr[idimL], NdimL*sizeof(dreal));
      dmrg_OprSvectBS(NdimBL, NdimSL, NdimSL, SZsiteL[iflav], vectLi, OprSLvectL);
      for(jdim = 0; jdim < Ndim; ++jdim) {
        jdimL = NdimL-jdim-1;
        memcpy(vectLj->addr, Dmat->ptr[jdimL], NdimL*sizeof(dreal));
        vOprv = vect_Dotprod_dreal(vectLj, OprSLvectL);
        SZblockmat[iflav]->ptr[idim][jdim] = vOprv;
      }
    }
    // Update SPsiteL
    for(idim = 0; idim < Ndim; ++idim) {
      idimL = NdimL-idim-1;
      memcpy(vectLi->addr, Dmat->ptr[idimL], NdimL*sizeof(dreal));
      dmrg_OprSvectBS(NdimBL, NdimSL, NdimSL, SPsiteL[iflav], vectLi, OprSLvectL);
      for(jdim = 0; jdim < Ndim; ++jdim) {
        jdimL = NdimL-jdim-1;
        memcpy(vectLj->addr, Dmat->ptr[jdimL], NdimL*sizeof(dreal));
        vOprv = vect_Dotprod_dreal(vectLj, OprSLvectL);
        SPblockmat[iflav]->ptr[idim][jdim] = vOprv;
      }
    }
  }
  
  vect_Del_dreal(OprLvectL);
  vect_Del_dreal(OprSLvectL);
  vect_Del_dreal(OprBLvectL);
  vect_Del_dreal(vectLj);
  vect_Del_dreal(vectLi);
  
  return;
  
 error:
  abort();
}
// Basis wavefunction update
static void dmrg_Wblock(mat2d_dreal *Wvc)
{
  int    ndim, mdim, idim, idimL;
  dreal  error;
  size_t nlen;
  
  ndim = Wvc->nrow;
  mdim = Wvc->ncol;
  
  nlen = ndim*sizeof(dreal);
  
  if(EigD) freeup(EigD);
  
  EigD = (dreal *) calloc(NdimL, sizeof(dreal));
  
  lapack_dsyev(NdimL, Dmat->addr, EigD);
  
  for(idimL = 0, error = 0.0; idimL < NdimL-mdim; ++idimL)
    error += EigD[idimL];
  
  if(Rank == Root) {
    printf("%s\n", dmrglabel);
    printf("%s error = %15.5e\n", dmrglabel, error);
  }
  
  for(idim = 0; idim < mdim; ++idim) {
    idimL = NdimL-idim-1;
    memcpy(Wvc->ptr[idim], Dmat->ptr[idimL], nlen);
  }
}

static void dmrg_EvcGS_Save(const vect_dreal *EvcGS_src, vect_dreal *EvcGS_des)
{
  int ndimBB;
  
  ndimBB = EvcGS_src->ndim;
  
  EvcGS_des->Deall(EvcGS_des);
  EvcGS_des->Alloc(ndimBB, EvcGS_des);
  
  vect_Copy_dreal(EvcGS_src, EvcGS_des);
}
// Calculate DenGS
static void dmrg_DenGS(const vect_dreal *EvcGS_save)
{
  int         ndimBL, ndimSL, ndimBR, ndimSR, ndimBB, ndimL, ndimR,
              indxBL, indxSL, indxBR, indxSR, iblock, iflav;
  vect_dreal  *EvcGS_old = NULL, *EvcGS_new = NULL, *OprEvcGS = NULL;
  mat2d_dreal **szsite = NULL, **occsite = NULL;
  
  DenGS->Reset(DenGS);
  
  ndimBL = Mretain[0];
  ndimBR = Mretain[Nblock-1];
  ndimSL = Nsite;
  ndimSR = Nsite;
  ndimBB = ndimBL*ndimSL*ndimBR*ndimSR;
  
  EvcGS_old = vect_New_dreal();
  EvcGS_old->Alloc(ndimBB, EvcGS_old);
  
  vect_Copy_dreal(EvcGS_save, EvcGS_old);
  
  EvcGS_new = vect_New_dreal();
  OprEvcGS = vect_New_dreal();
  
  szsite  = (mat2d_dreal **) calloc(Nflavor, sizeof(mat2d_dreal *));
  occsite = (mat2d_dreal **) calloc(Nflavor, sizeof(mat2d_dreal *));
  for(iflav = 0; iflav < Nflavor; ++iflav) {
    szsite[iflav]  = mat2d_New_dreal(); szsite[iflav]->Alloc(Nsite, Nsite, szsite[iflav]);
    occsite[iflav] = mat2d_New_dreal(); occsite[iflav]->Alloc(Nsite, Nsite, occsite[iflav]);
  }
  dmrg_Ssite(szsite, dmrg_SigmaZ);
  dmrg_Occmat(szsite, occsite);
  
  for(iblock = 0; iblock < Nblock; ++iblock) {
    indxBL = iblock;   indxBR = Length-iblock-4;
    indxSL = iblock+1; indxSR = iblock+2;
    
    ndimBL = Mretain[indxBL];
    ndimBR = Mretain[indxBR];
    ndimSL = Nsite;
    ndimSR = Nsite;
    ndimBB = ndimBL*ndimSL*ndimBR*ndimSR;
    ndimL  = ndimBL*ndimSL;
    ndimR  = ndimBR*ndimSR;
    
    EvcGS_new->Deall(EvcGS_new);
    EvcGS_new->Alloc(ndimBB, EvcGS_new);
    
    OprEvcGS->Deall(OprEvcGS);
    OprEvcGS->Alloc(ndimBB, OprEvcGS);
    
    if(indxBL == 0) {
      vect_Copy_dreal(EvcGS_old, EvcGS_new);
    } else {
      dmrg_EvcGS_Predict(indxBL, indxBR, EvcGS_old, EvcGS_new);
    }
    
    if(indxBL == 0) {
      for(iflav = 0; iflav < Nflavor; ++iflav) {
        OprEvcGS->Reset(OprEvcGS);
        dmrg_OprBLvect(ndimBL, ndimSL, ndimBR, ndimSR, ndimBL, occsite[iflav], EvcGS_new, OprEvcGS);
        DenGS->addr[indxBL] += vect_Dotprod_dreal(EvcGS_new, OprEvcGS);
      }
    }
    
    for(iflav = 0; iflav < Nflavor; ++iflav) {
      OprEvcGS->Reset(OprEvcGS);
      dmrg_OprSLvect(ndimBL, ndimSL, ndimBR, ndimSR, ndimSL, occsite[iflav], EvcGS_new, OprEvcGS);
      DenGS->addr[indxSL] += vect_Dotprod_dreal(EvcGS_new, OprEvcGS);
    }
    
    if(indxBR == 0) {
      for(iflav = 0; iflav < Nflavor; ++iflav) {
        OprEvcGS->Reset(OprEvcGS);
        dmrg_OprSRvect(ndimBL, ndimSL, ndimBR, ndimSR, ndimSR, occsite[iflav], EvcGS_new, OprEvcGS);
        DenGS->addr[indxSR] += vect_Dotprod_dreal(EvcGS_new, OprEvcGS);
	
	OprEvcGS->Reset(OprEvcGS);
	dmrg_OprBRvect(ndimBL, ndimSL, ndimBR, ndimSR, ndimBR, occsite[iflav], EvcGS_new, OprEvcGS);
        DenGS->addr[Length-1] += vect_Dotprod_dreal(EvcGS_new, OprEvcGS);
      }
    }
    
    dmrg_EvcGS_Save(EvcGS_new, EvcGS_old);
  }
  
  for(iflav = Nflavor-1; iflav > -1; --iflav) {
    mat2d_Del_dreal(occsite[iflav]);
    mat2d_Del_dreal(szsite[iflav]);
  }
  freeup(occsite);
  freeup(szsite);
  
  vect_Del_dreal(OprEvcGS);
  vect_Del_dreal(EvcGS_new);
  vect_Del_dreal(EvcGS_old);
}
    
// Save DMRG calculation results
static void dmrg_Output(void)
{
  FILE *fp = NULL;
  int  iblock, ndimBB, ndimW, mdimW;
  
  ndimBB = EvcGS->ndim;
  
  if(Rank == Root) {
    fp = fopen(evcgsbin, "w");
    fwrite(&ndimBB, sizeof(int), 1, fp);
    fwrite(EvcGS->addr, sizeof(dreal), ndimBB, fp);
    fclose(fp);
    
    fp = fopen(wblockbin, "w");
    for(iblock = 0; iblock < Nblock; ++iblock) {
      ndimW = Wblock[iblock]->nrow;
      mdimW = Wblock[iblock]->ncol;
      fwrite(&ndimW, sizeof(int), 1, fp);
      fwrite(&mdimW, sizeof(int), 1, fp);
      fwrite(Wblock[iblock]->addr, sizeof(dreal), ndimW*mdimW, fp);
    }
    fclose(fp);
  }
}
// Infinite one-dimensional DMRG calculations
static void dmrg_Infinite(void)
{
  int iblock, indxBL, indxBR, indxSL, indxSR;
  
  if(Rank == Root) {
    printf("%s\n", dmrglabel);
    printf("%s ======== DMRG Infinite ========\n", dmrglabel);
    printf("%s NblockInf = %5d\n", dmrglabel, NblockInf);
  }
  
  for(iblock = 0; iblock < NblockInf; ++iblock) {
    indxBL = iblock;   indxBR = iblock;
    indxSL = iblock+1; indxSR = iblock+2;
    
    if(indxBL == 0) {
      dmrg_Matsite(Echain[indxBL], Uchain[indxBL], SZblock[indxBL], SPblock[indxBL], Hblock[indxBL]);
    } else {
      dmrg_Matblock(SZblock[indxBL], SPblock[indxBL], Hblock[indxBL]);
    }
    dmrg_Superblock(indxBL, indxSL, indxBR, indxSR);
    dmrg_EvcGS_Init_Inf();
    davidson_Diag(NdimBB, &EigGS, EvcGS, dmrg_HBBvect, dmrg_PBBvect);
    if(Rank == Root) {
      printf("%s indxBL, indxSL = %5d %5d\n", dmrglabel, indxBL, indxSL);
      printf("%s indxBR, indxSR = %5d %5d\n", dmrglabel, indxBR, indxSR);
      printf("%s\n", dmrglabel);
      printf("%s EigGS = %15.5e\n", dmrglabel, EigGS);
      fflush(stdout);
    }
    dmrg_Dmat();
    dmrg_Wblock(Wblock[indxBL]);
  }
  
  if(Rank == Root) {
    printf("%s ===============================\n", dmrglabel);
    printf("%s\n", dmrglabel);
  }
}
// Remnant calculations after infinite, before finite DMRG calculations
static void dmrg_Remnant(void)
{
  int        ndimBL, ndimSL, ndimBR, ndimSR, ndimBB,
             iblock, indxBL, indxBR, indxSL, indxSR;
  vect_dreal *EvcGS_old = NULL;
  
  if(Rank == Root) {
    printf("%s\n", dmrglabel);
    printf("%s ========= DMRG Remnant ========\n", dmrglabel);
  }
  
  EvcGS_old = vect_New_dreal();
  
  ndimBB = EvcGS->ndim;
  EvcGS_old->Alloc(ndimBB, EvcGS_old);
  
  vect_Copy_dreal(EvcGS, EvcGS_old);
  
  for(iblock = NblockInf; iblock < Nblock; ++iblock) {
    indxBL = iblock;   indxBR = Length-iblock-4;
    indxSL = iblock+1; indxSR = iblock+2;
    
    if(indxBL == 0) {
      dmrg_Matsite(Echain[indxBL], Uchain[indxBL], SZblock[indxBL], SPblock[indxBL], Hblock[indxBL]);
    } else {
      dmrg_Matblock(SZblock[indxBL], SPblock[indxBL], Hblock[indxBL]);
    }
    dmrg_Superblock(indxBL, indxSL, indxBR, indxSR);
    dmrg_EvcGS_Init_Rem();
    davidson_Diag(NdimBB, &EigGS, EvcGS, dmrg_HBBvect, dmrg_PBBvect);
    
    dmrg_EvcGS_Save(EvcGS, EvcGS_old);
    if(Rank == Root) {
      printf("%s indxBL, indxSL = %5d %5d\n", dmrglabel, indxBL, indxSL);
      printf("%s indxBR, indxSR = %5d %5d\n", dmrglabel, indxBR, indxSR);
      printf("%s\n", dmrglabel);
      printf("%s EigGS = %15.5e\n", dmrglabel, EigGS);
      fflush(stdout);
    }
    dmrg_Dmat();
    dmrg_Wblock(Wblock[indxBL]);
  }
  
  if(Rank == Root) {
    printf("%s ===============================\n", dmrglabel);
    printf("%s\n", dmrglabel);
  }
  
  ndimBL = Mretain[0];
  ndimBR = Mretain[Nblock-1];
  ndimSL = Nsite;
  ndimSR = Nsite;
  ndimBB = ndimBL*ndimSL*ndimBR*ndimSR;
  
  EvcGS->Deall(EvcGS);
  EvcGS->Alloc(ndimBB, EvcGS);
  
  dmrg_Evc_Mirror(ndimBL, ndimSL, ndimBR, ndimSR, EvcGS_old, EvcGS);
  
  vect_Del_dreal(EvcGS_old);
}
// Finite one-dimensional DMRG sweep
static void dmrg_Sweep(void)
{
  int        ndimBL, ndimSL, ndimBR, ndimSR, ndimBB,
             iblock, indxBL, indxBR, indxSL, indxSR, len;
  vect_dreal *EvcGS_old = NULL;
  
  if(Rank == Root) {
    printf("%s\n", dmrglabel);
    printf("%s ========== DMRG Sweep =========\n", dmrglabel);
  }
  
  EvcGS_old = vect_New_dreal();
  
  for(iblock = 0; iblock < Nblock; ++iblock) {
    indxBL = iblock;   indxBR = Length-iblock-4;
    indxSL = iblock+1; indxSR = iblock+2;
    
    if(indxBL == 0) {
      dmrg_Matsite(Echain[indxBL], Uchain[indxBL], SZblock[indxBL], SPblock[indxBL], Hblock[indxBL]);
    } else {
      dmrg_Matblock(SZblock[indxBL], SPblock[indxBL], Hblock[indxBL]);
    }
    dmrg_Superblock(indxBL, indxSL, indxBR, indxSR);
    if(indxBL > 0) dmrg_EvcGS_Predict(indxBL, indxBR, EvcGS_old, EvcGS);
    davidson_Diag(NdimBB, &EigGS, EvcGS, dmrg_HBBvect, dmrg_PBBvect);
    
    dmrg_EvcGS_Save(EvcGS, EvcGS_old);
    if(Rank == Root) {
      printf("%s indxBL, indxSL = %5d %5d\n", dmrglabel, indxBL, indxSL);
      printf("%s indxBR, indxSR = %5d %5d\n", dmrglabel, indxBR, indxSR);
      printf("%s\n", dmrglabel);
      printf("%s EigGS = %15.5e\n", dmrglabel, EigGS);
      fflush(stdout);
    }
    dmrg_Dmat();
    dmrg_Wblock(Wblock[indxBL]);
  }
  
  ndimBL = Mretain[0];
  ndimBR = Mretain[Nblock-1];
  ndimSL = Nsite;
  ndimSR = Nsite;
  ndimBB = ndimBL*ndimSL*ndimBR*ndimSR;
  
  EvcGS->Deall(EvcGS);
  EvcGS->Alloc(ndimBB, EvcGS);
  
  dmrg_Evc_Mirror(ndimBL, ndimSL, ndimBR, ndimSR, EvcGS_old, EvcGS);
  
  dmrg_DenGS(EvcGS);
  
  if(Rank == Root) {
    printf("%s\n", dmrglabel);
    printf("%s DenGS =\n", dmrglabel);
    for(len = 0; len < Length; ++len)
      printf("%s %15.5e\n", dmrglabel, DenGS->addr[len]);
  }
  
  if(Rank == Root) {
    printf("%s ===============================\n", dmrglabel);
    printf("%s\n", dmrglabel);
  }
  
  vect_Del_dreal(EvcGS_old);
}
// Check convergence of ground state wavefunction with fixed basis
static void dmrg_Check(void)
{
  int        iblock, indxBL, indxBR, indxSL, indxSR;
  vect_dreal *EvcGS_old = NULL;
  
  if(Rank == Root) {
    printf("%s\n", dmrglabel);
    printf("%s ========== DMRG Check =========\n", dmrglabel);
  }
  
  EvcGS_old = vect_New_dreal();
  
  for(iblock = 0; iblock < Nblock; ++iblock) {
    indxBL = iblock;   indxBR = Length-iblock-4;
    indxSL = iblock+1; indxSR = iblock+2;
    
    if(indxBL == 0) {
      dmrg_Matsite(Echain[indxBL], Uchain[indxBL], SZblock[indxBL], SPblock[indxBL], Hblock[indxBL]);
    } else {
      dmrg_Matblock(SZblock[indxBL], SPblock[indxBL], Hblock[indxBL]);
    }
    dmrg_Superblock(indxBL, indxSL, indxBR, indxSR);
    if(indxBL > 0) dmrg_EvcGS_Predict(indxBL, indxBR, EvcGS_old, EvcGS);
    
    davidson_Diag(NdimBB, &EigGS, EvcGS, dmrg_HBBvect, dmrg_PBBvect);
    
    dmrg_EvcGS_Save(EvcGS, EvcGS_old);
    if(Rank == Root) {
      printf("%s indxBL, indxSL = %5d %5d\n", dmrglabel, indxBL, indxSL);
      printf("%s indxBR, indxSR = %5d %5d\n", dmrglabel, indxBR, indxSR);
      printf("%s\n", dmrglabel);
      printf("%s EigGS = %15.5e\n", dmrglabel, EigGS);
      fflush(stdout);
    }
    dmrg_Dmat_fixed(Wblock[indxBL]);
  }
  
  if(Rank == Root) {
    printf("%s ===============================\n", dmrglabel);
    printf("%s\n", dmrglabel);
  }
  
  vect_Del_dreal(EvcGS_old);
}

void dmrg_Finite(void)
{
  int isweep;
  
  dmrg_Infinite();
  dmrg_Remnant();
  for(isweep = 0; isweep < Nsweep; ++isweep) {
    if(Rank == Root)
      printf("%s ========== ISWEEP %3d ==========\n", dmrglabel, isweep);
    dmrg_Sweep();
  }
  dmrg_Output();
  dmrg_Check();
}
