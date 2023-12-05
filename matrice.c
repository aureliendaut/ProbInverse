#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "function.h"
#include "parametres.h"

void BuildSecondMember( double* u, double* b, int cas){


    int n = Nx*Ny;
    double dx = 1.0/(Nx+1);
    double dy = 1.0/(Ny+1); 
    
    double beta = -1.0/(dx*dx) ; 
    double gamma = -1.0/(dy*dy);

    for (int i = 0 ; i < Nx ; i++)
    {
        for (int j = 0 ; j < Ny ; j++)
        {
            b[i+ j*Nx] = Source((i+1)*dx, (j+1)*dy, cas) + u[i+j*Nx]; 
        }
    }
    for(int I=1; I<n-1; I++)
  {

    if(I<Nx-1)//Frontière du bas sans les deux coins inférieurs
    {
      b[I] = b[I] -gamma*UpBottomCondition((I+1)*dx,0.0,cas);
    }

    else if(I>n-Nx)//Frontière du haut sans les deux coins supérieurs
    {
      int i = I%Nx;
      b[I] = b[I] -gamma*UpBottomCondition((i+1)*dx,Ly,cas);

    }


    else if(((I+1)%Nx==0)&&(I!=Nx-1))//Coté droit sans les deux coins de droite
    {
        int j=I/Nx;
        b[I] = b[I] - beta*LeftRightCondition(Lx,(j+1)*dy,cas);
    }
    
    else if(((I+1)%Nx==1)&&(I!=n-Nx))//Coté gauche sans les deux coins de gauche
    {
        int j=I/Nx;
         b[I] = b[I] - beta*LeftRightCondition(0.0,(j+1)*dy,cas);
    }
}
    b[0] = b[0] - beta*LeftRightCondition(0.0,dy,cas) - gamma*UpBottomCondition(dx,0.0,cas);
    b[Nx-1] =b[Nx-1] - beta*LeftRightCondition(Lx,dy,cas) - gamma*UpBottomCondition(Nx*dx,0.0,cas);
    b[n-Nx] = b[n-Nx]  - beta*LeftRightCondition(0.0,Ny*dy,cas) - gamma*UpBottomCondition(dx,Ly,cas);
    b[n-1] =  b[n-1] - beta*LeftRightCondition(Lx,Ny*dy,cas) - gamma*UpBottomCondition(Nx*dx,Ly,cas);


 }   

void matvec(double *x , double*y) {
    double dx = Lx/(Nx+1);
    double dy = Ly/(Ny+1); 
    double alpha = (2.0/(dx*dx) + 2.0/(dy*dy)) ;  
    double beta = -1.0/(dx*dx) ; 
    double gamma = -1.0/(dy*dy);

    //Traitement des Nx premiers termes 
    y[0] =alpha*x[0] + beta*x[1] + gamma*x[Nx];
    for (int i = 1; i < Nx-1 ; i++){
        y[i] = alpha*x[i] + beta*(x[i-1] + x[i+1]) + gamma*x[i+Nx];
    }
    
    //On traite à part les termes en Nx-1 à cause du 0
    y[Nx-1] = alpha*x[Nx-1] + beta*x[Nx-2] + gamma*x[2*Nx-1];

    //Multiplication pour les lignes à 5 coefficients = pas dans les Nx premieres ou Nx dernieres lignes
    for (int i = Nx ; i < (Nx*(Ny-1)) ; i++){
        if (i%Nx == 0){
            //traite les lignes à 4 élements non nuls
            y[i] = alpha*x[i] + beta*x[i+1] + gamma*(x[i-Nx] + x[i+Nx]); 
        } 
        else if (i%Nx == Nx-1){
        //traite les lignes à 4 élements non nuls

            y[i] = alpha*x[i] + beta*x[i-1] + gamma*(x[i-Nx] + x[i+Nx]); 
        }
        else{
            y[i] = alpha*x[i] + beta*(x[i-1] + x[i+1]) + gamma*(x[i-Nx] + x[i+Nx]);
        }
    }
    //Traitement des Nx derniers termes 
     y[Nx*(Ny-1)] = alpha*x[Nx*(Ny-1)] +beta*x[Nx*(Ny-1) + 1] + gamma*x[Nx*(Ny-2)];
    for (int i = Nx*(Ny-1) + 1 ; i < Nx*Ny-1 ; i++){
        y[i] = alpha*x[i] + beta*(x[i-1] + x[i+1]) +gamma*x[i-Nx];
    }
    //On traite à part les termes en Nx-1 à cause du 0
    y[Nx*Ny -1] = alpha*x[Nx*Ny -1] + beta*x[Nx*Ny -2]  +gamma*x[Nx*(Ny-1) -1];
}