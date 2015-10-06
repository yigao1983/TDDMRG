#ifndef dmrg_h
#define dmrg_h

#include "matrix.h"

extern const int   NSP;
extern int         isGS, isTD;
extern int         Nblock, Mret, Nsweep, Nflavor, Nsite, NblockInf;
extern int         NdimBL, NdimSL, NdimBR, NdimSR, NdimBB, NdimSS, NdimL, NdimR;
extern int         *Mretain;
extern vect_dreal  *DenGS, *EvcGS;
extern mat2d_dreal **Wblock;

extern const char *evcgsbin;
extern const char *wblockbin;
//extern const char *hblockbin;
//extern const char *szblockbin;
//extern const char *spblockbin;

typedef void (*Paulimat)(mat2d_dreal *);

void dmrg_Init(void);
void dmrg_Quit(void);
void dmrg_Finite(void);

inline int dmrg_indxBB(int, int, int, int, int, int, int, int);
void dmrg_SigmaZ(mat2d_dreal *);
void dmrg_SigmaP(mat2d_dreal *);
void dmrg_Ssite(mat2d_dreal **, Paulimat);
void dmrg_Occmat(mat2d_dreal **, mat2d_dreal **);
void dmrg_Hmat(dreal, dreal, mat2d_dreal **, mat2d_dreal *);
void dmrg_Vmat(dreal, int, int, mat2d_dreal **, mat2d_dreal **, mat2d_dreal *);

#endif
