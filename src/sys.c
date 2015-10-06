#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <time.h>
#include "dbg.h"
#include "sys.h"

#ifdef MPI
#include <mpi.h>
#endif

int Size, Rank, Root;

static const char *syslabel = "sys:            ";

void sys_Init(int *argc, char ***argv)
{
  time_t current_time;
  char   time_start[SCHAR_MAX];
  
#ifdef MPI
  MPI_Init(argc, argv);
  MPI_Comm_size(MPI_COMM_WORLD, &Size);
  MPI_Comm_rank(MPI_COMM_WORLD, &Rank);
#else
  Size = 1;
  Rank = 0;
#endif
  Root = 0;
  
  time(&current_time);
  ctime_r(&current_time, time_start);
  
  if(Rank == Root) {
    printf("%s\n", syslabel);
    printf("%s start time: %s", syslabel, time_start);
    printf("%s\n", syslabel);
#ifdef MPI
    printf("%s running in parallel version\n", syslabel);
    printf("%s\n", syslabel);
#else
    printf("%s running in serial version\n", syslabel);
    printf("%s\n", syslabel);
#endif
    printf("%s Size: %d\n", syslabel, Size);
    printf("%s\n", syslabel);
  }
}

void sys_Quit(void)
{
  time_t current_time;
  char   time_end[SCHAR_MAX];
  
  sys_Sync();
  
  time(&current_time);
  ctime_r(&current_time, time_end);
  
  if(Rank == Root) {
    printf("%s\n", syslabel);
    printf("%s end time: %s", syslabel, time_end);
    printf("%s\n", syslabel);
  }
  
#ifdef MPI
  MPI_Finalize();
#endif
}

void sys_Sync(void)
{
#ifdef MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
}

void sys_Abort(void)
{
#ifdef MPI
  MPI_Abort(MPI_COMM_WORLD, 1);
#endif
  
  exit(EXIT_FAILURE);
}

void sys_Bcast_char(char *string, int nc, int root)
{
#ifdef MPI
  MPI_Bcast(string, nc, MPI_CHAR, root, MPI_COMM_WORLD);
#endif
}

void sys_Bcast_int(int *buffer, int nc, int root)
{
#ifdef MPI
  MPI_Bcast(buffer, nc, MPI_INT, root, MPI_COMM_WORLD);
#endif
}

void sys_Bcast_sreal(sreal *buffer, int nc, int root)
{
#ifdef MPI
  MPI_Bcast(buffer, nc, MPI_SREAL, root, MPI_COMM_WORLD);
#endif
}

void sys_Bcast_dreal(dreal *buffer, int nc, int root)
{
#ifdef MPI
  MPI_Bcast(buffer, nc, MPI_DREAL, root, MPI_COMM_WORLD);
#endif
}

void sys_Bcast_scmplx(scmplx *buffer, int nc, int root)
{
#ifdef MPI
  MPI_Bcast(buffer, nc, MPI_SCMPLX, root, MPI_COMM_WORLD);
#endif
}

void sys_Bcast_dcmplx(dcmplx *buffer, int nc, int root)
{
#ifdef MPI
  MPI_Bcast(buffer, nc, MPI_DCMPLX, root, MPI_COMM_WORLD);
#endif
}

void sys_Allreduce_sum_int(int *buffer, int nc)
{
#ifdef MPI
  int ic;
  int *buffer_sum = NULL;
  
  buffer_sum = (int *) calloc(nc, sizeof(int));
  check_mem(buffer_sum, "buffer_sum (int)");
  
  MPI_Allreduce(buffer, buffer_sum, nc, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  for(ic = 0; ic < nc; ++ic) {
    buffer[ic] = buffer_sum[ic];
  }
  
  freeup(buffer_sum);
  
  return;
  
 error:
  if(buffer_sum) freeup(buffer_sum);
  abort();
#endif
}

void sys_Allreduce_sum_sreal(sreal *buffer, int nc)
{
#ifdef MPI
  int   ic;
  sreal *buffer_sum = NULL;
  
  buffer_sum = (sreal *) calloc(nc, sizeof(sreal));
  check_mem(buffer_sum, "buffer_sum (sreal)");
  
  MPI_Allreduce(buffer, buffer_sum, nc, MPI_SREAL, MPI_SUM, MPI_COMM_WORLD);
  for(ic = 0; ic < nc; ++ic) {
    buffer[ic] = buffer_sum[ic];
  }
  
  freeup(buffer_sum);
  
  return;
  
 error:
  if(buffer_sum) freeup(buffer_sum);
  abort();
#endif
}

void sys_Allreduce_sum_dreal(dreal *buffer, int nc)
{
#ifdef MPI
  int   ic;
  dreal *buffer_sum = NULL;
  
  buffer_sum = (dreal *) calloc(nc, sizeof(dreal));
  check_mem(buffer_sum, "buffer_sum (dreal)");
  
  MPI_Allreduce(buffer, buffer_sum, nc, MPI_DREAL, MPI_SUM, MPI_COMM_WORLD);
  for(ic = 0; ic < nc; ++ic) {
    buffer[ic] = buffer_sum[ic];
  }
  
  freeup(buffer_sum);
  
  return;
  
 error:
  if(buffer_sum) freeup(buffer_sum);
  abort();
#endif
}

void sys_Allreduce_sum_scmplx(scmplx *buffer, int nc)
{
#ifdef MPI
  int    ic;
  scmplx *buffer_sum = NULL;
  
  buffer_sum = (scmplx *) calloc(nc, sizeof(scmplx));
  check_mem(buffer_sum, "buffer_sum (scmplx)");
  
  MPI_Allreduce(buffer, buffer_sum, nc, MPI_SCMPLX, MPI_SUM, MPI_COMM_WORLD);
  for(ic = 0; ic < nc; ++ic) {
    buffer[ic] = buffer_sum[ic];
  }
  
  freeup(buffer_sum);
  
  return;
  
 error:
  if(buffer_sum) freeup(buffer_sum);
  abort();
#endif
}

void sys_Allreduce_sum_dcmplx(dcmplx *buffer, int nc)
{
#ifdef MPI
  int    ic;
  dcmplx *buffer_sum = NULL;
  
  buffer_sum = (dcmplx *) calloc(nc, sizeof(dcmplx));
  check_mem(buffer_sum, "buffer_sum (dcmplx)");
  
  MPI_Allreduce(buffer, buffer_sum, nc, MPI_DCMPLX, MPI_SUM, MPI_COMM_WORLD);
  for(ic = 0; ic < nc; ++ic) {
    buffer[ic] = buffer_sum[ic];
  }
  
  freeup(buffer_sum);
  
  return;
  
 error:
  if(buffer_sum) freeup(buffer_sum);
  abort();
#endif
}
