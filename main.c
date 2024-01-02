#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "matrice.h"
#include "function.h"
#include "GC.h"
#include "parametres.h"

typedef double (*vecteur2D)[2];
typedef double vecteur[2];

//calcul de f avec les polynômes de Lagrange définis dans function.c
void calcul_f(vecteur2D a, vecteur2D P_a, vecteur X) 
{

    vecteur2D f;

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++) 
        {
            f = lagrange(a,P_a,X);
        }
    }
}




int main(int argc, char *argv[]) {
   if (argc != 2) {
      printf("Entrer le nom du fichier de paramètres \n");
      return 1;
   }

   FILE *file = fopen(argv[1], "r");
   if (file == NULL) {
      printf("Le fichier ne peut pas être ouvert'%s'\n", argv[1]);
      return 1;
   }
    //Lecture des paramètres dans le fichier
    char parameter[20];
    while (fscanf(file, "%s", parameter) != EOF) {
        if (strcmp(parameter, "Nx") == 0) {
            fscanf(file, "%d", &Nx);
        } else if (strcmp(parameter, "Ny") == 0) {
            fscanf(file, "%d", &Ny);
        } else if (strcmp(parameter, "Lx") == 0) {
            fscanf(file, "%lf", &Lx);
        } else if (strcmp(parameter, "Ly") == 0) {
            fscanf(file, "%lf", &Ly);
        }
        else if (strcmp(parameter, "cas") == 0) {
            fscanf(file, "%d", &cas);
        }
    }

    fclose(file);

    int N = Nx*Ny ; 
   //Vecteur solution à tn-1 et tn
   double* U0 = (double*) malloc(N* sizeof(double)); 
   double* U = (double*) malloc(N* sizeof(double));
   //vecteur 2nd membre
   double* b = (double*) malloc(N* sizeof(double));

   //vecteur erreur 
   double* erreur = (double*) malloc(N* sizeof(double));

    //vecteur sol exacte
    double* U_exacte = (double*) malloc(N* sizeof(double));

   //Initialisation
    double dx = 1.0/(Nx+1) ; 
    double dy = 1.0/(Ny+1) ; 
   for (int k =0 ; k<N ; k++){
      U0[k] = 0.0 ; 
      U[k] = 0.0 ; 
   }

   //Boucle principale
    int kmax = 10000;
    double eps = 1.0e-4 ; 
    double t = 0.0 ; 

    BuildSecondMember( U0,b, cas);  
    GC(b , U , eps, kmax) ; 



    FILE *fp_num, *fp_exacte;
    // Ouverture du fichier 
    char filename_num[100];
    char filename_exacte[100];
    sprintf(filename_num, "sol.dat"); 
    sprintf(filename_exacte, "sol_exacte.dat"); 
    fp_num = fopen(filename_num, "w");
    fp_exacte = fopen(filename_exacte, "w");

    // Écriture de la solution dans les fichiers et de la sol exacte

    for(int i=0;i<Nx;i++)
        {
        for(int j=0;j<Ny;j++)
            {
            fprintf(fp_num, "%f %f %f \n", (i+1)*dx,(j+1)*dy , U[j*Nx+i]);
            if (cas ==1){
                fprintf(fp_exacte, "%f %f %f \n", (i+1)*dx, (j+1)*dy ,  (((i+1)*dx)*(i+1)*dx -(i+1)*dx)*(((j+1)*dx)*(j+1)*dx -(j+1)*dx));
            }
            else if (cas ==2){
                fprintf(fp_exacte, "%f %f %f \n", (i+1)*dx, (j+1)*dy ,  sin((i+1)*dx) + cos((j+1)*dy)  );
                }
            }
        }

    // Fermeture du fichier
    fclose(fp_num);
    fclose(fp_exacte);

    //Mise à jour de Un
    for (int i = 0; i< N ; i++){
        U0[i] = U[i] ; 
    } 


    //Calcul de l'erreur
    // for (int i =0; i<N ; i++){
    //     erreur[i] = 0.0 ; 
    // }
    // for (int k =0; k<N ; k++){
    //     int i = k%Nx +1 ;
    //     int j = floor(k/Nx) + 1 ; 
    //     double x= i*dx ; 
    //     double y= j*dy ; 

    //     U_exacte[k] = (x*x -x)*(y*y -y);
    //     //U_exacte[k] = sin(x) + cos(y);
        
    //     erreur[k] = U[k] - U_exacte[k] ; 
    // }
    // printf("log erreur %f \n",log((norm(erreur,N)))) ; 
    // printf("log dx %f \n", log(dx)) ; 
    free (U0);
    free(U);
    free(b) ;  
    free(U_exacte) ; 
    free(erreur) ; 
   return 0;

}