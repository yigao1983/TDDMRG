#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include "dbg.h"
#include "sys.h"
#include "numeric.h"
#include "constants.h"
#include "parameters.h"

int   Nspin, Length, Bridge;
dreal *Echain = NULL, *Uchain = NULL, *Tchain = NULL;

static const char *parlabel = "parameters:     ";
static const char *chainin  = "chain.in";

static void parameters_Print(void)
{
  int len, brd;
  
  if(Rank == Root) {
    
    printf("%s Fermi chain parameters\n", parlabel);
    printf("%s\n", parlabel);
    printf("%s Nspin     = %5d\n", parlabel, Nspin);
    printf("%s Length    = %5d\n", parlabel, Length);
    printf("%s Bridge    = %5d\n", parlabel, Bridge);
    printf("%s E, U (Ha) = \n",    parlabel);
    for(len = 0; len < Length; ++len) {
      switch(Nspin) {
      case 1:
        printf("%s %5d %15.5e\n", parlabel, len, Echain[len]);
        fflush(stdout);
        break;
      case 2:
        printf("%s %5d %15.5e %15.5e\n", parlabel, len, Echain[len], Uchain[len]);
        break;
      default:
        abort();
      }
      fflush(stdout);
    }
    printf("%s T (Ha)    = \n",    parlabel);
    for(brd = 0; brd < Bridge; ++brd) {
      printf("%s %5d %15.5e\n", parlabel, brd, Tchain[brd]);
      fflush(stdout);
    }
  }
}

void parameters_Init(void)
{
  FILE *fp = NULL;
  int  lenhalf, brdhalf, len, brd;
  char line[SHRT_MAX];
  
  if(Rank == Root) {
    fp = fopen(chainin, "r");
    
    if(fgets(line, sizeof(line), fp)) sscanf(line, "%d %d", &Nspin, &Length);
    
    Echain = (dreal *) calloc(Length, sizeof(dreal));
    Uchain = (dreal *) calloc(Length, sizeof(dreal));
    
    lenhalf = (Length + Length%2) / 2;
    for(len = 0; len < lenhalf; ++len) {
      if(fgets(line, sizeof(line), fp)) {
        switch(Nspin) {
        case 1:
          sscanf(line, "%lf", Echain+len);
          break;
        case 2:
          sscanf(line, "%lf %lf", Echain+len, Uchain+len);
          break;
        default:
          abort();
        }
        Echain[Length-len-1] = Echain[len];
        Uchain[Length-len-1] = Uchain[len];
      }
    }
    
    Bridge = Length - 1;
    
    Tchain = (dreal *) calloc(Bridge, sizeof(dreal));
    
    brdhalf = (Bridge + Bridge%2) / 2;
    for(brd = 0; brd < brdhalf; ++brd) {
      if(fgets(line, sizeof(line), fp)) {
        sscanf(line, "%lf", Tchain+brd);
        Tchain[Bridge-brd-1] = Tchain[brd];
      }
    }
    
    fclose(fp);
  }
  
  sys_Bcast_int(&Nspin,  1, Root);
  sys_Bcast_int(&Length, 1, Root);
  sys_Bcast_int(&Bridge, 1, Root);
  
  check(Length > 3, "Invalid Length");
  
  if(Rank != Root) {
    Echain = (dreal *) calloc(Length, sizeof(dreal));
    Uchain = (dreal *) calloc(Length, sizeof(dreal));
    Tchain = (dreal *) calloc(Bridge, sizeof(dreal));
  }
  
  sys_Bcast_dreal(Echain, Length, Root);
  sys_Bcast_dreal(Uchain, Length, Root);
  sys_Bcast_dreal(Tchain, Bridge, Root);
  
  for(len = 0; len < Length; ++len) {
    Echain[len] /= EV_HA;
    Uchain[len] /= EV_HA;
  }
  
  for(brd = 0; brd < Bridge; ++brd)
    Tchain[brd] /= EV_HA;
  
  parameters_Print();
  
  return;
  
 error:
  if(Tchain) freeup(Tchain);
  if(Uchain) freeup(Uchain);
  if(Echain) freeup(Echain);
  abort();
}

void parameters_Quit(void)
{
  if(Tchain) freeup(Tchain);
  if(Uchain) freeup(Uchain);
  if(Echain) freeup(Echain);
}

