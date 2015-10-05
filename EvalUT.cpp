#include "EvalUT.h"
#include "iostream"

void Hamiltonian(problem& Problem,matrix& _new,matrix& _old,double zt,double MuA,int mynode) {
    
    int hx=0;
    int theta=0;
    double U12=0.6;
    double zV=0.6;
    int U0=1;
    double MuB;
    MuB=MuA;
    matrix A=Problem.H[0];
    matrix B=Problem.H[1];
    matrix Na=Problem.H[2];
    matrix Nb=Problem.H[3];
    matrix Ap=transpose(A);
    matrix Bp=transpose(B);
    matrix Va(Problem.Com_D,Problem.Com_D),Ea(1,Problem.Com_D);
    matrix Vb(Problem.Com_D,Problem.Com_D),Eb(1,Problem.Com_D);

    

    
    
    matrix I(Problem.Com_D,Problem.Com_D);
    I.eye();
    matrix Ha(Problem.Com_D,Problem.Com_D),Hb(Problem.Com_D,Problem.Com_D);

    
    for (int i=0; i<_new.rows; i++) {
        zeros(Ha);
        zeros(Hb);
        
        Ha=-zt*((_old(i,2)*Ap+_old(i,2)*A)+(_old(i,3)*Bp+_old(i,3)*B))+hx*(Ap*B+Bp*A)-(MuA*Na+MuB*Nb)+(0.5*U0)*(Na*(Na-I))+(0.5*U0)*(Nb*(Nb-I))+U12*Na*Nb+zV*(_old(i,6)+_old(i,7))*(Na+Nb)+I*((zt/2)*(_old(i,2)*_old(i,0)+_old(i,3)*_old(i,1)))-I*(0.5*zV*(_old(i,4)+_old(i,5))*(_old(i,6)+_old(i,7)));
        
        

        
        
        Hb=-zt*((_old(i,0)*Ap+_old(i,0)*A)+(_old(i,1)*Bp+_old(i,1)*B))+hx*(Ap*B+Bp*A)-(MuB*Na+MuB*Nb)+(0.5*U0)*(Na*(Na-I))+(0.5*U0)*(Nb*(Nb-I))+U12*Na*Nb+zV*(_old(i,4)+_old(i,5))*(Na+Nb)+I*((zt/2)*(_old(i,2)*_old(i,0)+_old(i,3)*_old(i,1)))-I*(0.5*zV*(_old(i,4)+_old(i,5))*(_old(i,6)+_old(i,7)));
        
        
        eig(Ea,Va,Ha);
        eig(Eb,Vb,Hb);
        matrix Va0c=Va.coli(0);
        matrix Va0c_conj=vec_transpose(Va.coli(0));
        matrix Vb0c=Vb.coli(0);
        matrix Vb0c_conj=vec_transpose(Vb.coli(0));

        
        matrix tmp(1,1);
        tmp=Va0c_conj*A*Va0c;
        _new(i,0)=tmp(0,0);
        tmp=Va0c_conj*B*Va0c;
        _new(i,1)=tmp(0,0);
        tmp=Va0c_conj*Na*Va0c;
        _new(i,4)=tmp(0,0);
        tmp=Va0c_conj*Nb*Va0c;
        _new(i,5)=tmp(0,0);
        
        tmp=Vb0c_conj*A*Vb0c;
        _new(i,2)=tmp(0,0);
        tmp=Vb0c_conj*B*Vb0c;
        _new(i,3)=tmp(0,0);
        tmp=Vb0c_conj*Na*Vb0c;
        _new(i,6)=tmp(0,0);
        tmp=Vb0c_conj*Nb*Vb0c;
        _new(i,7)=tmp(0,0);
    }
    
    

}



void old_new_norm(matrix& _new,matrix& _old,matrix& _out) {
    
    int r=_new.rows;
    int c=_new.cols;
//    matrix tmp1(1,c),tmp2(1,c);
    for (int i=0; i<r; i++) {
        //tmp1=_new.rowi(i);
        //tmp2=_old.rowi(i);
        _out.array1d[i]=norm2(_new.rowi(i)-_old.rowi(i));
    }
}


double norm2(matrix vec) {
    
    int one=1;
    int r=vec.cols;
    double f;

    return (f=dnrm2(&r,vec.array1d, &one));
    
}


