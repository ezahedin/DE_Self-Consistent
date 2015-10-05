#include "utilities.h"
//=======================================================================
void rand_matrix(int a,int* seed,matrix& t1)
{
    
    int rows=t1.rows;
    int cols=t1.cols;
    int size=rows*cols;
    dlarnv (&a, seed, &size,t1.array1d);
    
}

//=======================================================================
void random_shuffle1(int* mat,int size,int mynode)
{
    srand (unsigned(time(0)*mynode));
    std::random_shuffle( &mat[0], &mat[size]);
   /* for (int i=0; i<size; i++)
        std::cout<<mat[i]<<'\n';*/
}

//=======================================================================
void constrain(matrix& pop,double min,double max)
{
    int rows=pop.rows;
    int cols=pop.cols;
    int size=cols*rows;
    for (int i=0; i<size; i++) {
        if(pop.array1d[i]>max){
            pop.array1d[i]=max;
        } else if (pop.array1d[i]<min) {
            pop.array1d[i]=min;
        }
    }
}
//=======================================================================
matrix projection(matrix& t1)
{
    int G[20]={0 ,1, 2, 3, 4, 5, 6, 8, 9, 12, 16, 17, 18, 20, 21, 24, 32, 33, 36, 48};
    matrix tmp(20,20);
    
    for (int k=0; k<20; k++)
        for (int j=0; j<20; j++)
            tmp(k,j)=t1(G[k],G[j]);
    
    return tmp;
}
//=======================================================================
void seed_fun(int* seed,int mynode) {
    if (fmod((double)mynode,2)==0) {
        seed[0]=mynode+61;seed[1]=mynode+234;seed[2]=mynode+3021;seed[3]=mynode+2401;
    } else {
        seed[0]=mynode+100;seed[1]=mynode+2074;seed[2]=mynode+1101;seed[3]=mynode+904;
    }
}
