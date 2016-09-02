/* Propensity file for the Wang-model
*  P. Bauer 2012-12-05
*/

#include <stdlib.h>
#include "propensities.h"

#define NR 3

#define C1 0.01
#define C2 0.01
#define C3 10

/* state of A is conserved */
#define xA 1000

/* other states vary */
#define B 0
#define C 1


/* Reaction propensities. */
double rFun1(const int *x, double t, double vol, const double *data, int sd)
{
  return C1*xA*x[B]/vol; //A+B->2B+A
}

double rFun2(const int *x, double t, double vol, const double *data, int sd)
{
  return C2*x[B]*x[C]/vol; //C+B->2C
}

double rFun3(const int *x, double t, double vol, const double *data, int sd)
{
  return C3*x[C]; //C->0
}

PropensityFun *ALLOC_propensities(void)
{
  PropensityFun *ptr = (PropensityFun *)malloc(sizeof(PropensityFun)*NR);
  
  ptr[0]=rFun1;
  ptr[1]=rFun2;
  ptr[2]=rFun3;
    
  return ptr;
}

void FREE_propensities(PropensityFun* ptr)
{
  free(ptr);
}