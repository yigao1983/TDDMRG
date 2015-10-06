#ifndef precision_h
#define precision_h

#include <complex.h>

typedef float          sreal;
typedef double         dreal;
typedef float  complex scmplx;
typedef double complex dcmplx;

#ifdef MPI

#define MPI_SREAL  MPI_FLOAT
#define MPI_DREAL  MPI_DOUBLE
#define MPI_SCMPLX MPI_C_FLOAT_COMPLEX
#define MPI_DCMPLX MPI_C_DOUBLE_COMPLEX

#endif

#endif
