#ifndef DE_P
#define DE_P

//DE takes the problem as the input.
//problem include all the information that DE needs for optimization
#include "problem.h"
#include "algorithm"
#include "utilities.h"
#include "EvalUT.h"
#include "printi.h"
#include <string>
//void constrain(matrix* pop,float min,float max,unsigned int popsize,unsigned int vecsize);
//unsigned int uniquerand(trng::yarn2 *R, unsigned int i1,unsigned int i2);
void DE(problem& Problem,matrix& out,int mynode,int* seed,double& val,double zt,double MuA);
#endif