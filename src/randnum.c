#include <stdlib.h>
#include "randnum.h"

static const unsigned int seed = 1103515245 + 12345;

void inline randnum_Set(void)
{
  srand(seed);
}

dreal inline randnum_Get(void)
{
  dreal rr;
  
  rr = (dreal) rand() / RAND_MAX;
  
  return rr;
}
