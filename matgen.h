#ifndef MATGEN_H
#define MATGEN_H




#include "matrix.h"
#include <iostream>
#include <math.h>
#include <vector>

typedef std::vector<matrix> mVector;


int dim(const int n);
void Ladder(matrix& a,int n,int lattice);
void NumberOp(matrix& a,int n,int lattice);

#endif