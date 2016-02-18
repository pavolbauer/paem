/* Propensity file for test "spinpaem14"
*  P. Bauer 2015-01-20
*/

#include <stdlib.h>
#include "propensities.h"

#define NR 2

/* Reaction propensities. */
double rFun1(const int *x, double t, double vol, const double *data, int sd)
{
  return x[0]*1/vol;
}

double rFun2(const int *x, double t, double vol, const double *data, int sd)
{
  return x[1]*1/vol;
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