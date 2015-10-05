#include "matgen.h"

int dim(const int n) {
    int num=0;
    for (int k=0;k<=n;k++) {
        for(int i=0;i<=n;i++) {
            for(int j=0;j<=n;j++) {
                if((i+j)==k) {
                    num++;
                }
            }
        }
    }
    
    return num;
}


void Ladder(matrix& a,int n,int lattice) {
    
    int num=0;
    mVector Vector;
    matrix tmp(1,2);
    matrix  p1(1,2);
    for (int k=0;k<=n;k++) {
        for(int i=0;i<=n;i++) {
            for(int j=0;j<=n;j++) {
                if((i+j)==k) {
                    tmp(0,0)=i;
                    tmp(0,1)=j;
                    Vector.push_back(tmp);
                }
            }
        }
    }
    
    int Vsize=Vector.size();
    
    //std::cout<<Vsize<<'\n';
    
    zeros(a);
    
    switch (lattice) {
        case 1:
            for (int i=0; i<Vsize; i++) {
                p1(0,0)=Vector[i](0,0)-1;
                p1(0,1)=Vector[i](0,1);
                for (int k=0; k<Vsize; k++) {
                    matrix p2=Vector[k];
                    if (p1==p2) {
                        a(k,i)=sqrt(Vector[i](0,0));
                    }
                }
            }

            break;
        case 2:

            for (int i=0; i<Vsize; i++) {
                p1(0,0)=Vector[i](0,0);
                p1(0,1)=Vector[i](0,1)-1;
                for (int k=0; k<Vsize; k++) {
                    matrix p2=Vector[k];
                    if (p1==p2) {
                        a(k,i)=sqrt(Vector[i](0,1));
                    }
                }
            }
            
            break;

        default:
            std::cout<<"There is no matrix to be generated \n";
            break;
    }
    
}


void NumberOp(matrix& a,int n,int lattice) {
    
    int num=1;
    mVector Vector;
    matrix tmp(1,2);
    for (int k=0;k<=n;k++) {
        for(int i=0;i<=n;i++) {
            for(int j=0;j<=n;j++) {
                if((i+j)==k) {
                    tmp(0,0)=i;
                    tmp(0,1)=j;
                    Vector.push_back(tmp);
                }
            }
        }
    }
    
    int Vsize=Vector.size();
    zeros(a);
    
    switch (lattice) {
        case 1:
            for (int i=0; i<Vsize; i++)
                for (int k=0; k<Vsize; k++)
                    if (i==k)
                        a(k,i)=Vector[i](0,0);

            
            break;
        case 2:
            
            for (int i=0; i<Vsize; i++)
                for (int k=0; k<Vsize; k++)
                    if (i==k)
                        a(k,i)=Vector[i](0,1);
            
            break;
            
        default:
            std::cout<<"There is no matrix to be generated \n";
            break;
    }
    
    
    
}





