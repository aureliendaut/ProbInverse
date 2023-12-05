#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrice.h"
#include "parametres.h"


// produit scalaire
double dot_product(double* a, double* b) {
    int n = Nx*Ny;
    double result = 0.0;
    for (int i = 0; i < n; i++) {
        result  = result+ a[i] * b[i];
    }
    return result;
}

// Norme 
double norm(double* vector, int n) {
    double result = 0.0;
    for (int i = 0; i < n; i++) {
        result  =  result + vector[i] * vector[i];
    }
    return sqrt(result);
}


// Gradient Conjugué
void GC(double* b, double* X, double eps, int kmax) {
  int n = Nx*Ny ;
  double* r_k = (double*)malloc(n * sizeof(double)); 
  double* r_kplus = (double*)malloc(n * sizeof(double));
  double* z = (double*)malloc(n * sizeof(double));
  double* Y = (double*)malloc(n * sizeof(double));
  double* p = (double*)malloc(n * sizeof(double));
  double beta ; 
  double alpha ; 
  double gamma ; 
  
  for (int i = 0; i < n; i++) {
    X[i] = 0.0;
  }
  matvec(X,Y);
  for (int i = 0; i < n; i++) {
    r_k[i] = b[i] - Y[i]; 
    p[i] = r_k[i] ; 
  }
  beta =  norm ( r_k , n ); 
  int k =0;
  
  while (beta > eps && k <= kmax) 
  {
    matvec(p,z); 
    alpha = dot_product(r_k, r_k) / dot_product(z, p);
      
    for (int i=0; i < n ; i++ ){
      X[i]=X[i] + alpha*p[i];
      r_kplus[i]=r_k[i]-alpha*z[i];
    }

    gamma = dot_product(r_kplus, r_kplus) / dot_product(r_k, r_k);

    for (int i=0 ; i < n ; i++){
      p[i]=r_kplus[i]+ gamma*p[i];

    }
    for (int i=0 ; i < n ; i++){
      r_k[i]=r_kplus[i];

    }
    beta = norm(r_k, n);
    k+=1;
  }


  free(r_k);
  free (r_kplus);
  free(z) ; 
  free(p) ; 
  free(Y);
}