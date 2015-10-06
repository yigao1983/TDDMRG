#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dbg.h"
#include "sys.h"
#include "matrix.h"

mat1d_dreal *mat1d_New_dreal(void)
{
  mat1d_dreal *mat1d = NULL;
  
  mat1d = (mat1d_dreal *) calloc(1, sizeof(mat1d_dreal));
  check_mem(mat1d, "mat1d");
  
  mat1d->Alloc = mat1d_Alloc_dreal;
  mat1d->Deall = mat1d_Deall_dreal;
  mat1d->Print = mat1d_Print_dreal;
  mat1d->Reset = mat1d_Reset_dreal;
  
  return mat1d;
  
 error:
  if(mat1d) freeup(mat1d);
  return NULL;
}

mat1d_dcmplx *mat1d_New_dcmplx(void)
{
  mat1d_dcmplx *mat1d = NULL;
  
  mat1d = (mat1d_dcmplx *) calloc(1, sizeof(mat1d_dcmplx));
  check_mem(mat1d, "mat1d");
  
  mat1d->Alloc = mat1d_Alloc_dcmplx;
  mat1d->Deall = mat1d_Deall_dcmplx;
  mat1d->Print = mat1d_Print_dcmplx;
  mat1d->Reset = mat1d_Reset_dcmplx;
  
  return mat1d;
  
 error:
  if(mat1d) freeup(mat1d);
  return NULL;
}

void mat1d_Alloc_dreal(int ndim, mat1d_dreal *mat1d)
{
  dreal *addr = NULL;
  
  mat1d->ndim = ndim;
  
  addr = (dreal *) calloc(ndim, sizeof(dreal));
  check_mem(addr, "addr");
  
  mat1d->addr = addr;
  
  return;
  
 error:
  if(addr) freeup(addr);
  abort();
}

void mat1d_Alloc_dcmplx(int ndim, mat1d_dcmplx *mat1d)
{
  dcmplx *addr = NULL;
  
  mat1d->ndim = ndim;
  
  addr = (dcmplx *) calloc(ndim, sizeof(dcmplx));
  check_mem(addr, "addr");
  
  mat1d->addr = addr;
  
  return;
  
 error:
  if(addr) freeup(addr);
  abort();
}

void mat1d_Deall_dreal(mat1d_dreal *mat1d)
{
  if(mat1d) {
    if(mat1d->addr) {
      freeup(mat1d->addr);
    }
  }
}

void mat1d_Deall_dcmplx(mat1d_dcmplx *mat1d)
{
  if(mat1d) {
    if(mat1d->addr) {
      freeup(mat1d->addr);
    }
  }
}

void mat1d_Print_dreal(const char *filename, const mat1d_dreal *mat1d)
{
  FILE *fp = NULL;
  int  ndim, idim;
  
  ndim = mat1d->ndim;
  
  if(Rank == Root) {
    fp = fopen(filename, "w");
    for(idim = 0; idim < ndim; ++idim) {
      fprintf(fp, "%15.5e\n", mat1d->addr[idim]);
    }
    fclose(fp);
  }
}

void mat1d_Print_dcmplx(const char *filename, const mat1d_dcmplx *mat1d)
{
  FILE *fp = NULL;
  int  ndim, idim;
  
  ndim = mat1d->ndim;

  if(Rank == Root) {
    fp = fopen(filename, "w");
    for(idim = 0; idim < ndim; ++idim) {
      fprintf(fp, "(%15.5e %15.5e)\n",
              creal(mat1d->addr[idim]),
              cimag(mat1d->addr[idim]));
    }
    fclose(fp);
  }
}

void mat1d_Reset_dreal(mat1d_dreal *mat1d)
{
  int    ndim;
  size_t nlen;
  
  check_mem(mat1d, "mat1d");
  
  ndim = mat1d->ndim;
  nlen = ndim * sizeof(dreal);
  
  memset(mat1d->addr, 0, nlen);
  
  return;
  
 error:
  abort();
}

void mat1d_Reset_dcmplx(mat1d_dcmplx *mat1d)
{
  int    ndim;
  size_t nlen;
  
  check_mem(mat1d, "mat1d");
  
  ndim = mat1d->ndim;
  nlen = ndim * sizeof(dcmplx);
  
  memset(mat1d->addr, 0, nlen);
  
  return;
  
 error:
  abort();
}

void mat1d_Copy_dreal(const mat1d_dreal *mat1d_src, mat1d_dreal *mat1d_des)
{
  int    ndim_src, ndim_des;
  size_t nlen;
  
  check_mem(mat1d_src, "mat1d_src");
  check_mem(mat1d_des, "mat1d_des");
  ndim_src = mat1d_src->ndim;
  ndim_des = mat1d_des->ndim;
  check(ndim_src == ndim_des, "Inequivalent ndim");
  
  nlen = ndim_src * sizeof(dreal);
  
  memcpy(mat1d_des->addr, mat1d_src->addr, nlen);
  
  return;
  
 error:
  abort();
}

void mat1d_Copy_dcmplx(const mat1d_dcmplx *mat1d_src, mat1d_dcmplx *mat1d_des)
{
  int    ndim_src, ndim_des;
  size_t nlen;
  
  check_mem(mat1d_src, "mat1d_src");
  check_mem(mat1d_des, "mat1d_des");
  ndim_src = mat1d_src->ndim;
  ndim_des = mat1d_des->ndim;
  check(ndim_src == ndim_des, "Inequivalent ndim");
  
  nlen = ndim_src * sizeof(dcmplx);
  
  memcpy(mat1d_des->addr, mat1d_src->addr, nlen);
  
  return;
  
 error:
  abort();
}

mat2d_dreal *mat2d_New_dreal(void)
{
  mat2d_dreal *mat2d = NULL;
  
  mat2d = (mat2d_dreal *) calloc(1, sizeof(mat2d_dreal));
  check_mem(mat2d, "mat2d");
  
  mat2d->Alloc = mat2d_Alloc_dreal;
  mat2d->Deall = mat2d_Deall_dreal;
  mat2d->Print = mat2d_Print_dreal;
  mat2d->Reset = mat2d_Reset_dreal;
  
  return mat2d;
  
 error:
  if(mat2d) freeup(mat2d);
  return NULL;
}

mat2d_dcmplx *mat2d_New_dcmplx(void)
{
  mat2d_dcmplx *mat2d = NULL;
  
  mat2d = (mat2d_dcmplx *) calloc(1, sizeof(mat2d_dcmplx));
  check_mem(mat2d, "mat2d");
  
  mat2d->Alloc = mat2d_Alloc_dcmplx;
  mat2d->Deall = mat2d_Deall_dcmplx;
  mat2d->Print = mat2d_Print_dcmplx;
  mat2d->Reset = mat2d_Reset_dcmplx;
  
  return mat2d;
  
 error:
  if(mat2d) freeup(mat2d);
  return NULL;
}

void mat2d_Alloc_dreal(int nrow, int ncol, mat2d_dreal *mat2d)
{
  int   icol;
  dreal *addr = NULL;
  dreal **ptr = NULL;
  
  mat2d->nrow = nrow;
  mat2d->ncol = ncol;
  
  addr = (dreal *) calloc(nrow * ncol, sizeof(dreal));
  check_mem(addr, "addr");
  ptr  = (dreal **) calloc(ncol, sizeof(dreal *));
  check_mem(ptr, "ptr");
  
  for(icol = 0; icol < ncol; ++icol) {
    ptr[icol] = addr + icol * nrow; 
  }
  
  mat2d->addr = addr;
  mat2d->ptr  = ptr;
  
  return;
  
 error:
  if(addr) freeup(addr);
  if(ptr)  freeup(ptr);
  abort();
}

void mat2d_Alloc_dcmplx(int nrow, int ncol, mat2d_dcmplx *mat2d)
{
  int    icol;
  dcmplx *addr = NULL;
  dcmplx **ptr = NULL;
  
  mat2d->nrow = nrow;
  mat2d->ncol = ncol;
  
  addr = (dcmplx *) calloc(nrow * ncol, sizeof(dcmplx));
  check_mem(addr, "addr");
  ptr  = (dcmplx **) calloc(ncol, sizeof(dcmplx *));
  check_mem(ptr, "ptr");
  
  for(icol = 0; icol < ncol; ++icol) {
    ptr[icol] = addr + icol * nrow; 
  }
  
  mat2d->addr = addr;
  mat2d->ptr  = ptr;
  
  return;
  
 error:
  if(addr) freeup(addr);
  if(ptr)  freeup(ptr);
  abort();
}

void mat2d_Deall_dreal(mat2d_dreal *mat2d)
{
  if(mat2d) {
    if(mat2d->addr) freeup(mat2d->addr);
    if(mat2d->ptr)  freeup(mat2d->ptr);
  }
}

void mat2d_Deall_dcmplx(mat2d_dcmplx *mat2d)
{
  if(mat2d) {
    if(mat2d->addr) freeup(mat2d->addr);
    if(mat2d->ptr)  freeup(mat2d->ptr);
  }
}

void mat2d_Print_dreal(const char *filename, const mat2d_dreal *mat2d)
{
  FILE *fp = NULL;
  int  nrow, irow, ncol, icol;
  
  nrow = mat2d->nrow;
  ncol = mat2d->ncol;
  
  if(Rank == Root) {
    fp = fopen(filename, "w");
    for(irow = 0; irow < nrow; ++irow) {
      for(icol = 0; icol < ncol; ++icol) {
        fprintf(fp, " %15.5e ", mat2d->ptr[icol][irow]);
      }
      fprintf(fp, "\n");
    }
    fclose(fp);
  }
}

void mat2d_Print_dcmplx(const char *filename, const mat2d_dcmplx *mat2d)
{
  FILE *fp = NULL;
  int  nrow, irow, ncol, icol;
  
  nrow = mat2d->nrow;
  ncol = mat2d->ncol;
  
  if(Rank == Root) {
    fp = fopen(filename, "w");
    for(irow = 0; irow < nrow; ++irow) {
      for(icol = 0; icol < ncol; ++icol) {
        fprintf(fp, " (%15.5e %15.5e) ",
                creal(mat2d->ptr[icol][irow]),
                cimag(mat2d->ptr[icol][irow]));
      }
      fprintf(fp, "\n");
    }
    fclose(fp);
  } 
}

void mat2d_Reset_dreal(mat2d_dreal *mat2d)
{
  int    nrow, ncol;
  size_t nlen;
  
  check_mem(mat2d, "mat2d");
  
  nrow = mat2d->nrow;
  ncol = mat2d->ncol;
  nlen = nrow * ncol * sizeof(dreal);
  
  memset(mat2d->addr, 0, nlen);
  
  return;
  
 error:
  abort();
}

void mat2d_Reset_dcmplx(mat2d_dcmplx *mat2d)
{
  int    nrow, ncol;
  size_t nlen;
  
  check_mem(mat2d, "mat2d");
  
  nrow = mat2d->nrow;
  ncol = mat2d->ncol;
  nlen = nrow * ncol * sizeof(dcmplx);
  
  memset(mat2d->addr, 0, nlen);
  
  return;
  
 error:
  abort();
}

void mat2d_Copy_dreal(const mat2d_dreal *mat2d_src, mat2d_dreal *mat2d_des)
{
  int    nrow_src, nrow_des, ncol_src, ncol_des;
  size_t nlen;
  
  check_mem(mat2d_src, "mat2d_src");
  check_mem(mat2d_des, "mat2d_des");
  nrow_src = mat2d_src->nrow; ncol_src = mat2d_src->ncol;
  nrow_des = mat2d_des->nrow; ncol_des = mat2d_des->ncol;
  check(nrow_src == nrow_des, "Inequivalent nrow");
  check(ncol_src == ncol_des, "Inequivalent ncol");
  
  nlen = nrow_src * ncol_src * sizeof(dreal);
  
  memcpy(mat2d_des->addr, mat2d_src->addr, nlen);
  
  return;
  
 error:
  abort();
}

void mat2d_Copy_dcmplx(const mat2d_dcmplx *mat2d_src, mat2d_dcmplx *mat2d_des)
{
  int    nrow_src, nrow_des, ncol_src, ncol_des;
  size_t nlen;
  
  check_mem(mat2d_src, "mat2d_src");
  check_mem(mat2d_des, "mat2d_des");
  nrow_src = mat2d_src->nrow; ncol_src = mat2d_src->ncol;
  nrow_des = mat2d_des->nrow; ncol_des = mat2d_des->ncol;
  check(nrow_src == nrow_des, "Inequivalent nrow");
  check(ncol_src == ncol_des, "Inequivalent ncol");
  
  nlen = nrow_src * ncol_src * sizeof(dcmplx);
  
  memcpy(mat2d_des->addr, mat2d_src->addr, nlen);
  
  return;
  
 error:
  abort();
}

void mat2d_Transpose_dreal(const mat2d_dreal *mat2d_src, mat2d_dreal *mat2d_des)
{
  int nrow_src, nrow_des, ncol_src, ncol_des, irow, icol;
  
  check_mem(mat2d_src, "mat2d_src");
  check_mem(mat2d_des, "mat2d_des");
  nrow_src = mat2d_src->nrow; ncol_src = mat2d_src->ncol;
  nrow_des = mat2d_des->nrow; ncol_des = mat2d_des->ncol;
  check(nrow_src == ncol_des, "Invalid ncol_des");
  check(ncol_src == nrow_des, "Invalid nrow_des");
  
  for(icol = 0; icol < ncol_src; ++icol) {
    for(irow = 0; irow < nrow_src; ++irow) {
      mat2d_des->ptr[irow][icol] = mat2d_src->ptr[icol][irow];
    }
  }
  
  return;
  
 error:
  abort();
}

void mat2d_Transpose_dcmplx(const mat2d_dcmplx *mat2d_src, mat2d_dcmplx *mat2d_des)
{
  int nrow_src, nrow_des, ncol_src, ncol_des, irow, icol;
  
  check_mem(mat2d_src, "mat2d_src");
  check_mem(mat2d_des, "mat2d_des");
  nrow_src = mat2d_src->nrow; ncol_src = mat2d_src->ncol;
  nrow_des = mat2d_des->nrow; ncol_des = mat2d_des->ncol;
  check(nrow_src == ncol_des, "Inequivalent nrow");
  check(ncol_src == nrow_des, "Inequivalent ncol");
  
  for(icol = 0; icol < ncol_src; ++icol) {
    for(irow = 0; irow < nrow_src; ++irow) {
      mat2d_des->ptr[irow][icol] = mat2d_src->ptr[icol][irow];
    }
  }
  
  return;
  
 error:
  abort();
}

// For debug
/*
#include "lapack.h"

void mat2d_check(void)
{
  const int ndim = 150;
  int       idim, jdim;
  mat2d_dcmplx *mat_a = NULL, *mat_b = NULL, *mat_c = NULL;
  
  mat_a = mat2d_New_dcmplx();
  check_mem(mat_a, "mat_a");
  mat_b = mat2d_New_dcmplx();
  check_mem(mat_b, "mat_b");
  mat_c = mat2d_New_dcmplx();
  check_mem(mat_c, "mat_c");
  
  mat_a->Alloc(ndim, ndim, mat_a);
  mat_b->Alloc(ndim, ndim, mat_b);
  mat_c->Alloc(ndim, ndim, mat_c);
  
  for(idim = 0 ; idim < ndim; ++idim) {
    for(jdim = 0; jdim < ndim; ++jdim) {
      mat_a->ptr[idim][jdim] = (dreal) rand() / RAND_MAX;
    }
  }
  
  mat2d_Copy_dcmplx(mat_a, mat_b);
  
  lapack_zinvs(ndim, mat_a->addr);
  
  lapack_zgemm(ndim, ndim, ndim, 'N', 'N', 1.0, mat_a->addr, mat_b->addr, 0.0, mat_c->addr);
  
  mat_a->Print("mat_a.txt", mat_a);
  mat_b->Print("mat_b.txt", mat_b);
  mat_c->Print("mat_c.txt", mat_c);
  
  mat2d_del_dcmplx(mat_a);
  mat2d_del_dcmplx(mat_b);
  mat2d_del_dcmplx(mat_c);
  
  return;
  
 error:
  if(mat_a) freeup(mat_a);
  if(mat_b) freeup(mat_b);
  if(mat_c) freeup(mat_c);
  abort();
}
*/
