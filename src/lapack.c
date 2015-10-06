#include <stdio.h>
#include <stdlib.h>
#include "dbg.h"
#include "sys.h"
#include "matrix.h"
#include "lapack.h"

// Interface to lapack routine dgemm
// mm: nrow of A
// nn: ncol of B
// kk: ncol and nrow of C
void lapack_dgemm(int mm, int nn, int kk,
                  char transa, char transb,
                  dreal alpha, dreal *AA, dreal *BB,
                  dreal beta, dreal *CC)
{
  int lda, ldb, ldc;
  
  if(transa == 'N' || transa == 'n') {
    lda = (1 > mm) ? 1 : mm;
  } else {
    lda = (1 > kk) ? 1 : kk;
  }
  
  if(transb == 'N' || transb == 'n') {
    ldb = (1 > kk) ? 1 : kk;
  } else {
    ldb = (1 > nn) ? 1 : nn;
  }
  
  ldc = (1 > mm) ? 1 : mm;
  
  dgemm_(&transa, &transb, &mm, &nn, &kk, &alpha, AA, &lda, BB, &ldb, &beta, CC, &ldc); 
}

// Interface to lapack routine dgetrf
// mm: nrow of A
// nn: ncol of A
void lapack_dgetrf(int mm, int nn, dreal *AA, int *ipiv)
{
  int lda, info;
  
  lda = (1 > mm) ? 1 : mm;
  
  dgetrf_(&mm, &nn, AA, &lda, ipiv, &info);
  check(info == 0, "Failed dgetrf, info = %d", info);
  
  return;
  
 error:
  abort();
}

// Interface to lapack routine dgetri
// nn: nrow and ncol of A
void lapack_dgetri(int nn, dreal *AA, int *ipiv)
{
  int   lda, lwork, info;
  dreal *work = NULL;
  
  lda = (1 > nn) ? 1 : nn;
  lwork = lda;
  
  work = (dreal *) calloc(lwork, sizeof(dreal));
  check_mem(work, "work");
  
  dgetri_(&nn, AA, &lda, ipiv, work, &lwork, &info);
  check(info == 0, "Failed dgetri, info = %d", info);
  
  freeup(work);
  
  return;
  
 error:
  if(work) freeup(work);
  abort();
}

// Inverse of matrix AA
void lapack_dinvs(int nn, dreal *AA)
{
  int *ipiv = NULL;
  
  ipiv = (int *) calloc(nn, sizeof(int));
  check_mem(ipiv, "ipiv");
  
  lapack_dgetrf(nn, nn, AA, ipiv);
  lapack_dgetri(nn, AA, ipiv);
  
  freeup(ipiv);
  
  return;
  
 error:
  if(ipiv) freeup(ipiv);
  abort();
}

// Interface to lapack routine zgemm
// mm: nrow of A
// nn: ncol of B
// kk: ncol and nrow of C
void lapack_zgemm(int mm, int nn, int kk,
                  char transa, char transb,
                  dcmplx alpha, dcmplx *AA, dcmplx *BB,
                  dcmplx beta, dcmplx *CC)
{
  int lda, ldb, ldc;
  
  if(transa == 'N' || transa == 'n') {
    lda = (1 > mm) ? 1 : mm;
  } else {
    lda = (1 > kk) ? 1 : kk;
  }
  
  if(transb == 'N' || transb == 'n') {
    ldb = (1 > kk) ? 1 : kk;
  } else {
    ldb = (1 > nn) ? 1 : nn;
  }
  
  ldc = (1 > mm) ? 1 : mm;
  
  zgemm_(&transa, &transb, &mm, &nn, &kk, &alpha, AA, &lda, BB, &ldb, &beta, CC, &ldc); 
}

// Interface to lapack routine zgetrf
// mm: nrow of A
// nn: ncol of A
void lapack_zgetrf(int mm, int nn, dcmplx *AA, int *ipiv)
{
  int lda, info;
  
  lda = (1 > mm) ? 1 : mm;
  
  zgetrf_(&mm, &nn, AA, &lda, ipiv, &info);
  check(info == 0, "Failed zgetrf, info = %d", info);
  
  return;
  
 error:
  abort();
}

// Interface to lapack routine zgetri
// nn: nrow and ncol of A
void lapack_zgetri(int nn, dcmplx *AA, int *ipiv)
{
  int     lda, lwork, info;
  dcmplx  *work = NULL;
  
  lda = (1 > nn) ? 1 : nn;
  lwork = lda;
  
  work = (dcmplx *) calloc(lwork, sizeof(dcmplx));
  check_mem(work, "work");
  
  zgetri_(&nn, AA, &lda, ipiv, work, &lwork, &info);
  check(info == 0, "Failed zgetri, info = %d", info);
  
  freeup(work);
  
  return;
  
 error:
  if(work) freeup(work);
  abort();
}

// Inverse of matrix AA
void lapack_zinvs(int nn, dcmplx *AA)
{
  int *ipiv = NULL;
  
  ipiv = (int *) calloc(nn, sizeof(int));
  check_mem(ipiv, "ipiv");
  
  lapack_zgetrf(nn, nn, AA, ipiv);
  lapack_zgetri(nn, AA, ipiv);
  
  freeup(ipiv);
  
  return;
  
 error:
  if(ipiv) freeup(ipiv);
  abort();
}

void lapack_dsteqr(int nn, int ldz, dreal *alph, dreal *beta, dreal *zz)
{
  int   nwork, info;
  char  compz = 'I';
  dreal *work = NULL;
  
  nwork = (1 >= 2*nn-2) ? 1 : 2*nn-2;
  
  work = (dreal *) calloc(nwork, sizeof(dreal));
  check_mem(work, "work");
  
  dsteqr_(&compz, &nn, alph, beta, zz, &ldz, work, &info);
  
  freeup(work);
  
  return;
  
 error:
  if(work) freeup(work);
  abort();
}

void lapack_dsyev(int nn, dreal *AA, dreal *ww)
{
  int   lda, lwork, info;
  char  jobz = 'V', uplo = 'U';
  dreal *work = NULL;
  
  lda = (1 > nn) ? 1 : nn;
  lwork = (1 > 3*nn-1) ? 1 : 3*nn-1;
  work = (dreal *) calloc(lwork, sizeof(dreal));
  check_mem(work, "work");
  
  dsyev_(&jobz, &uplo, &nn, AA, &lda, ww, work, &lwork, &info);
  
  freeup(work);
  
  return;
  
 error:
  if(work) freeup(work);
  abort();
}

void lapack_dsyevr(int nn, char range, dreal vl, dreal vu, int il, int iu, dreal abstol,
                   dreal *AA, dreal *ww)
{
  int   lda, ldz, mm, lwork, liwork, info, ii;
  char  jobz = 'V', uplo = 'U';
  int   *isuppz = NULL, *iwork = NULL;
  dreal *ZZ = NULL, *work = NULL;
  
  lda = (1 > nn) ? 1 : nn;
  ldz = lda;
  lwork = (1 > 26*nn) ? 1 : 26*nn;
  liwork = (1 > 10*nn) ? 1 : 10*nn;
  ZZ = (dreal *) calloc(ldz*nn, sizeof(dreal));
  check_mem(ZZ, "ZZ");
  isuppz = (int *) calloc(2*lda, sizeof(int));
  check_mem(isuppz, "isuppz");
  work = (dreal *) calloc(lwork, sizeof(dreal));
  check_mem(work, "work");
  iwork = (int *) calloc(liwork, sizeof(int));
  
  dsyevr_(&jobz, &range, &uplo, &nn, AA, &lda, &vl, &vu, &il, &iu,
          &abstol, &mm, ww, ZZ, &ldz, isuppz, work, &lwork,
          iwork, &liwork, &info);
  
  for(ii = 0; ii < mm; ++mm)
    memcpy(AA+ii*nn, ZZ+ii*nn, nn*sizeof(dreal));
  
  freeup(work);
  freeup(isuppz);
  freeup(ZZ);
  
  return;
  
 error:
  if(work) freeup(work);
  if(isuppz) freeup(isuppz);
  if(ZZ) freeup(ZZ);
  abort();
}

void lapack_zheev(int nn, dcmplx *AA, dreal *ww)
{
  int    lda, lwork, rsize, info;
  char   jobz = 'V', uplo = 'U';
  dreal  *rwork = NULL;
  dcmplx *work = NULL;
  
  lda = (1 > nn) ? 1 : nn;
  lwork = (1 > 2*nn-1) ? 1 : 2*nn-1;
  rsize = (1 > 3*nn-2) ? 1 : 3*nn-2;
  work = (dcmplx *) calloc(lwork, sizeof(dcmplx));
  check_mem(work, "work");
  rwork = (dreal *) calloc(rsize, sizeof(dreal));
  check_mem(rwork, "rwork");
  
  zheev_(&jobz, &uplo, &nn, AA, &lda, ww, work, &lwork, rwork, &info);
  
  freeup(rwork);
  freeup(work);
  
  return;
  
 error:
  if(rwork) freeup(rwork);
  if(work) freeup(work);
  abort();
}
