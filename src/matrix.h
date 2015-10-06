#ifndef matrix_h
#define matrix_h

#include "precision.h"

typedef struct mat1d_dreal
{
  int   ndim;
  dreal *addr;
  void (*Alloc)(int, struct mat1d_dreal *);
  void (*Deall)(struct mat1d_dreal *);
  void (*Print)(const char *, const struct mat1d_dreal *);
  void (*Reset)(struct mat1d_dreal *);
} mat1d_dreal;

typedef struct mat1d_dcmplx
{
  int    ndim;
  dcmplx *addr;
  void (*Alloc)(int, struct mat1d_dcmplx *);
  void (*Deall)(struct mat1d_dcmplx *);
  void (*Print)(const char *, const struct mat1d_dcmplx *);
  void (*Reset)(struct mat1d_dcmplx *);
} mat1d_dcmplx;

typedef struct mat2d_dreal
{
  int   nrow, ncol;
  dreal *addr;
  dreal **ptr;
  void (*Alloc)(int, int, struct mat2d_dreal *);
  void (*Deall)(struct mat2d_dreal *);
  void (*Print)(const char *, const struct mat2d_dreal *);
  void (*Reset)(struct mat2d_dreal *);
} mat2d_dreal;

typedef struct mat2d_dcmplx
{
  int    nrow, ncol;
  dcmplx *addr;
  dcmplx **ptr;
  void (*Alloc)(int, int, struct mat2d_dcmplx *);
  void (*Deall)(struct mat2d_dcmplx *);
  void (*Print)(const char *, const struct mat2d_dcmplx *);
  void (*Reset)(struct mat2d_dcmplx *);
} mat2d_dcmplx;

mat1d_dreal  *mat1d_New_dreal(void);
mat1d_dcmplx *mat1d_New_dcmplx(void);

void mat1d_Alloc_dreal(int, mat1d_dreal *);
void mat1d_Alloc_dcmplx(int, mat1d_dcmplx *);

void mat1d_Deall_dreal(mat1d_dreal *);
void mat1d_Deall_dcmplx(mat1d_dcmplx *);

void mat1d_Print_dreal(const char *, const mat1d_dreal *);
void mat1d_Print_dcmplx(const char *, const mat1d_dcmplx *);

void mat1d_Reset_dreal(mat1d_dreal *);
void mat1d_Reset_dcmplx(mat1d_dcmplx *);

void mat1d_Copy_dreal(const mat1d_dreal *, mat1d_dreal *);
void mat1d_Copy_dcmplx(const mat1d_dcmplx *, mat1d_dcmplx *);

mat2d_dreal  *mat2d_New_dreal(void);
mat2d_dcmplx *mat2d_New_dcmplx(void);

void mat2d_Alloc_dreal(int, int, mat2d_dreal *);
void mat2d_Alloc_dcmplx(int, int, mat2d_dcmplx *);

void mat2d_Deall_dreal(mat2d_dreal *);
void mat2d_Deall_dcmplx(mat2d_dcmplx *);

void mat2d_Print_dreal(const char *, const mat2d_dreal *);
void mat2d_Print_dcmplx(const char *, const mat2d_dcmplx *);

void mat2d_Reset_dreal(mat2d_dreal *);
void mat2d_Reset_dcmplx(mat2d_dcmplx *);

void mat2d_Copy_dreal(const mat2d_dreal *, mat2d_dreal *);
void mat2d_Copy_dcmplx(const mat2d_dcmplx *, mat2d_dcmplx *);

void mat2d_Transpose_dreal(const mat2d_dreal *, mat2d_dreal *);
void mat2d_Transpose_dcmplx(const mat2d_dcmplx *, mat2d_dcmplx *);

#define mat1d_Del_dreal(A)  { mat1d_Deall_dreal(A);  freeup(A); }
#define mat1d_Del_dcmplx(A) { mat1d_Deall_dcmplx(A); freeup(A); }
#define mat2d_Del_dreal(A)  { mat2d_Deall_dreal(A);  freeup(A); }
#define mat2d_Del_dcmplx(A) { mat2d_Deall_dcmplx(A); freeup(A); }

#endif
