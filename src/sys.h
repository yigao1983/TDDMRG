#ifndef sys_h
#define sys_h

#include <stdio.h>
#include "precision.h"

extern int Size, Rank, Root;

void sys_Init(int *, char ***);
void sys_Quit(void);
void sys_Sync(void);
void sys_Abort(void);
void sys_Bcast_char(char *, int, int);
void sys_Bcast_int(int *, int, int);
void sys_Bcast_dreal(dreal *, int, int);
void sys_Bcast_dcmplx(dcmplx *, int, int);
void sys_Allreduce_sum_int(int *, int);
void sys_Allreduce_sum_dreal(dreal *, int);
void sys_Allreduce_sum_dcmplx(dcmplx *, int);

#define abort() { printf("[ABORT] (%s:%s) Rank = %d\n", __FILE__, __FUNCTION__, Rank); sys_Abort(); }

#endif
