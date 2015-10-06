#include <stdlib.h>
#include "sys.h"
#include "parameters.h"
#include "dmrg.h"
#include "tddmrg.h"

int main(int argc, char *argv[])
{
  sys_Init(&argc, &argv);
  
  parameters_Init();
  // DMRG computation
  dmrg_Init();
  if(isGS) dmrg_Finite();
  // TDDMRG computation
  if(isTD) {
    tddmrg_Init();
    tddmrg_Propagate();
    tddmrg_Quit();
  }
  // end of TDDMRG
  dmrg_Quit();
  // end of DMRG
  parameters_Quit();
  
  sys_Quit();
  
  return EXIT_SUCCESS;
}
