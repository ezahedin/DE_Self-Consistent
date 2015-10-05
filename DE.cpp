#include "DE.h"
void DE(problem& Problem,matrix& out,int mynode,int* seed,double& val,double zt,double MuA) {
    //===============================Initializing DE parameters===========================================
    //Population size NPDIM
    int dist1=1;
    int dist2=2;
    unsigned int NPDIM=(Problem.NPDim)*(Problem.K);
    unsigned int each_popsize=Problem.K*Problem.M;
    //mutatted population
    matrix _oldpop(NPDIM,each_popsize);
    matrix _newpop(NPDIM,each_popsize);
    matrix _mpop(NPDIM,each_popsize);
    matrix _newcpop(NPDIM,each_popsize);
    matrix _postnewcpop(NPDIM,each_popsize);
    //Crossed over elements
    matrix index(1,each_popsize);
    matrix fed(1,NPDIM);
    matrix tmpfed(1,NPDIM);
    
    

    
    //Initialization of populations
    rand_matrix(dist1,seed,_oldpop);
    Hamiltonian(Problem,_newpop,_oldpop,zt,MuA,mynode);

    old_new_norm(_newpop,_oldpop,fed);

    double F=Problem.F;
    double CR=Problem.CR;
    
    int min_id=fed.min();
    int G=1;
    int indexid[NPDIM];

    for (int i=0; i<NPDIM; i++)
        indexid[i]=i;
    
    //out=_newpop.rowi(min_id);

    for (int i=0; i<each_popsize; i++) {
        out(0,i)=_newpop(min_id,i);
    }
    
    
    
    while (G<Problem.G) {
        
        
        
        if (fmod(double(G),5000)==0) {
            rand_matrix(dist1,seed,_oldpop);
            Hamiltonian(Problem,_newpop,_oldpop,zt,MuA,mynode);
            
            old_new_norm(_newpop,_oldpop,fed);
            
        }
        

        min_id=fed.min();


        
        if(fed.array1d[min_id]<1e-5)
            break;
        
        //Mutation part
        for (int i=0; i<NPDIM; i++) {
            //rand_matrix(dist1,seed,index);
            random_shuffle1(indexid,NPDIM,seed[2]);
            
            
            if ((indexid[0]==i)|(indexid[1]==i)|(indexid[2]==i)) {
                int tmp1=indexid[0];int tmp2=indexid[1];int tmp3=indexid[2];
                indexid[0]=indexid[3];indexid[1]=indexid[4];indexid[2]=indexid[5];
                indexid[3]=tmp1;indexid[4]=tmp2;indexid[5]=tmp3;
            }
            
            for (int j=0; j<each_popsize; j++) {
                _mpop(i,j)=_newpop(min_id,j)+F*(_newpop(min_id,j)-_newpop(indexid[1],j))+F*(_newpop(min_id,j)-_newpop(indexid[0],j))+F*(_newpop(indexid[2],j)-_newpop(indexid[0],j));
            }
        }
        
        
        constrain(_mpop,0,20);
        
        _newcpop=_newpop;

        
        for (int i=0; i<NPDIM; i++) {
            rand_matrix(dist1,seed,index);
            for (int j=0; j<each_popsize; ++j)
                if (index.array1d[j]<CR)
                    _newcpop(i,j)=_mpop(i,j);
        }
        
        
        Hamiltonian(Problem,_postnewcpop,_newcpop,zt,MuA,mynode);
        old_new_norm(_postnewcpop,_newcpop,tmpfed);

        for (int i=0; i<NPDIM; i++) {
            if (tmpfed.array1d[i]<fed.array1d[i]) {
                fed.array1d[i]=tmpfed.array1d[i];
                    for (int j=0; j<each_popsize; j++)
                        _newpop(i,j)=_newcpop(i,j);}
        
            }
        
        G++;
    }
    //out=_newpop.rowi(min_id);
    
    for (int i=0; i<each_popsize; i++) {
        out(0,i)=_newpop(min_id,i);
    }
    

    
    val=fed.array1d[min_id];
    
}


