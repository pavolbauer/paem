/* Propensity file for test "spin1"
*  P. Bauer 2012-03-26
*/

#include <stdlib.h>
#include "propensities.h"

#define NR 1

/* Reaction propensities. */
double rFun1(const int *x, double t, double vol, const double *data, int sd)
{
  return x[0]*x[1]/vol;
}

PropensityFun *ALLOC_propensities(void)
{
  PropensityFun *ptr = (PropensityFun *)malloc(sizeof(PropensityFun)*NR);
  
  ptr[0]=rFun1;
    
  return ptr;
}

void FREE_propensities(PropensityFun* ptr)
{
  free(ptr);
}