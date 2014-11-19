/* Propensity file for the Jeschke-model
*  P. Bauer 2014-11-13
*/

#include <stdlib.h>
#include "propensities.h"

#define NR 1

/* Reaction propensities. */
double rFun1(const int *x, double t, double vol, const double *data, int sd)
{
  return 1000; // 2*1000*0.5
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