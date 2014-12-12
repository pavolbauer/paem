/* Propensity file for test "spin7"
*  P. Bauer 2012-03-28
*/

#include <stdlib.h>
#include "propensities.h"

#define NR 2

/* Reaction propensities. */
double rFun1(const int *x, double t, double vol, const double *data, int sd)
{
  return x[0]*0.15/vol;
}

double rFun2(const int *x, double t, double vol, const double *data, int sd)
{
  return x[1]*0.15/vol;
}

PropensityFun *ALLOC_propensities(void)
{
  PropensityFun *ptr = (PropensityFun *)malloc(sizeof(PropensityFun)*NR);
  
  ptr[0]=rFun1;
  ptr[1]=rFun2;
    
  return ptr;
}

void FREE_propensities(PropensityFun* ptr)
{
  free(ptr);
}