/* Propensity file for test "spin7"
*  P. Bauer 2012-03-28
*/

#include <stdlib.h>
#include "propensities.h"

#define NR 2

/* Reaction propensities. */
double rFun1(const int *x, double t, double vol, const double *data, int sd)
{
  return x[0]*625; // R:D ~ 1:1
}

double rFun2(const int *x, double t, double vol, const double *data, int sd)
{
  return x[1]*625; // R:D ~ 1:1
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