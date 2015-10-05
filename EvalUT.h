#include<iostream>
#include "matrix.h"
#include "problem.h"
#include "printi.h"

void Hamiltonian(problem& Problem,matrix& _new,matrix& _old,double zt,double MuA,int mynode);
double norm2(matrix vec);
void old_new_norm(matrix& _new,matrix& _old,matrix& _out);