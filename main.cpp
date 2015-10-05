#include "cmatrix.h"
#include "kron.h"
#include <math.h>
#include <iostream>
#include "utilities.h"
#include <time.h>
#include "DE.h"
#include "problem.h"
#include "EvalUT.h"
#include "printi.h"
#include <mpi.h>
#include "matgen.h"

//0  0  0  0  0.5  0.5  0.5  0.5  0  0  0.865169

int main()
{
    //==============================================================Initializing MPI Gloabl Variable=================================================================

    int totalnode,mynode;
    MPI::Init();
    MPI_Comm_size(MPI_COMM_WORLD,&totalnode);
    MPI_Comm_rank(MPI_COMM_WORLD,&mynode);
    MPI_Status status;

    //==============================================================Initializing Variables for DE====================================================================
    unsigned int Generation=15000;
    float Crossover_CR=0.5;
    float Mutation_F=0.9;
    double NPDim_p=1;
    double ctrpar=8;
    //defining system paremeters
    unsigned int D=20;   //Dimension of each qubits
    unsigned int boson_pop=10;
    unsigned int Com_D=dim(10);
    unsigned int PM=1;  //numner of control Hamiltonain
    int seed[4];
    seed_fun(seed,mynode);

    

    
    //========================================================Declaring individual Hamiltonains (Real)===================================================


    
    matrix A(Com_D,Com_D),B(Com_D,Com_D),Na(Com_D,Com_D),Nb(Com_D,Com_D);
    Ladder(A,boson_pop,1);
    Ladder(B,boson_pop,2);
    NumberOp(Na,boson_pop,1);
    NumberOp(Nb,boson_pop,2);
    
    
  

    mVector Htot(4,matrix(Com_D,Com_D));Htot[0]=A;Htot[1]=B;Htot[2]=Na;Htot[3]=Nb;
    //=====================================Define the problem to be sent to DE function, Pre-initialization of DE optimization===========================
    problem Problem(PM,Crossover_CR,Mutation_F,D,NPDim_p,ctrpar,Generation,Htot,Com_D);
    
    double lengthzt=0.4;
    double lengthmu=30;
    double grid1=0.002;//.002
    double grid2=0.075;//.05
    int mu_dim=lengthmu/grid2+1; //11
    matrix out(1,8);
    int num=0;
    double val;
    
    int Num_of_points=lengthzt/grid1+1; //41
    int Points_per_node=Num_of_points/totalnode; //5
    double Interval_per_node=Points_per_node*grid1; //.05
    
    
    
    
    int lastnode_zt_dim=round((lengthzt-(double(totalnode-1))*(Interval_per_node))/grid1+1); //6

    
    int checkpoint=lastnode_zt_dim-Points_per_node; //1
    
    
    //make correction for the last node
    int total_points_per_node[totalnode];
    matrix interval_mat(1,totalnode);
    for (int i=0; i<totalnode; i++)
        if (i<checkpoint) {
            total_points_per_node[i]=1+Points_per_node;
            interval_mat(0,i)=total_points_per_node[i]*grid1;
        } else {
            total_points_per_node[i]=Points_per_node;
            interval_mat(0,i)=total_points_per_node[i]*grid1;
        }
    
//total_points_per_node matrix=[6 5 5 5 5 5 5 5]
//interval_mat=[.06 .05 .05 .05 .05 .05 .05 .05]
    
    
    
    Points_per_node=total_points_per_node[mynode];
    matrix result(Points_per_node*mu_dim,11);
    
    double total_dis_sofar=0.0;
    for (int i=0; i<mynode; i++) {
        total_dis_sofar+=(interval_mat(0,i));
    }
    
    
    
    int total_dis_sofar_int=0;
    for (int i=0; i<mynode; i++) {
        total_dis_sofar_int+=total_points_per_node[i];
    }
    
    for (int i=total_dis_sofar_int; i<(total_dis_sofar_int+Points_per_node);i++) {
        
        double zt=i*grid1;
        
        
        for (int j=0; j<mu_dim; j++) {
            double mu=j*grid2;


            DE(Problem,out,mynode,seed,val,zt,mu);
            
            
            
            for (int k=0; k<8; k++)
                result(num,k)=out(0,k);
            
            result(num,8)=zt;
            result(num,9)=mu;
            result(num,10)=val;
            num++;
        }
    }

    
    
    
    
    
    //matrix to be sent to node zero
    matrix send1(total_points_per_node[0]*mu_dim,11); //for processors with higer nodes
    matrix send2(total_points_per_node[totalnode-1]*mu_dim,11); //for processores with less nodes

    //print out the results
    if (mynode!=0) {
        
            MPI_Send(result.array1d,Points_per_node*mu_dim*11,MPI_DOUBLE, 0, 0,MPI_COMM_WORLD);
        
    } else {
        printi(result,"Results.dat");

        for (int j=1; j<totalnode; j++) {
            if (j<checkpoint) {
                MPI_Recv(send1.array1d,total_points_per_node[0]*mu_dim*11, MPI_DOUBLE, j, 0, MPI_COMM_WORLD,&status);
                printi(send1,"Results.dat");
            } else {
                MPI_Recv(send2.array1d,total_points_per_node[totalnode-1]*mu_dim*11, MPI_DOUBLE, j, 0, MPI_COMM_WORLD,&status);
                printi(send2,"Results.dat");
            }
    }
    }


    MPI::Finalize();


    return 0;
}


