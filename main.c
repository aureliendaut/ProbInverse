#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "matrice.h"
#include "function.h"
#include "GC.h"
#include <time.h>
#include "parametres.h"

#define NB_RESTARTS 1
#define NB_ITER_MAX 1000
#define BATCH_SIZE 5
#define GRADIENT_EPSILON 1e-5
#define GIVE_UP_THRESHOLD 1e6
#define ADAM_EPSILON 1e-8
#define LEARNING_RATE 0.01

double lagrange2d(int n, double *alpha, double x, double y) {
    double result = 0.0;
    for (int k = 0; k < n; k++) {
        double basis_polynomial = 1.0;
        for (int j = 0; j < n; j++) {
            if (!(alpha[3 * j] == alpha[3 * k])) {
                basis_polynomial *= (x - alpha[3 * j]) / (alpha[3 * k] - alpha[3 * j]);
            }
            if (!(alpha[3 * j + 1] == alpha[3 * k + 1])) {
                basis_polynomial *= (y - alpha[3 * j + 1]) / (alpha[3 * k + 1] - alpha[3 * j + 1]);
            }
        }
        result += basis_polynomial * alpha[3 * k + 2];
    }
    return result;
}

// Fonction pour calculer f en utilisant le polynôme de Lagrange
void compute_f_lagrange(int n, double *alpha, double dx, double dy, int Nx, int Ny, double *f) {
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            f[i + Nx * j] = lagrange2d(n, alpha, (i + 1) * dx, (j + 1) * dy);
        }
    }
}

// Fonction pour obtenir des échantillons aléatoires
void get_random_samples(double *U, int m, double dx, double dy, int Nx, int Ny, double NOISE_MU, double NOISE_SIGMA, double *T, double *X, double *Y, double *Gamma) {
    int N = Nx * Ny;
    double *Xtilde = (double *)malloc(m * sizeof(double));
    double *Ytilde = (double *)malloc(m * sizeof(double));
    int *I = (int *)malloc(m * sizeof(int));
    int *J = (int *)malloc(m * sizeof(int));

    // Remplir X et Y avec des valeurs aléatoires
    for (int i = 0; i < m; i++) {
        X[i] = (double)rand() / RAND_MAX;
        Y[i] = (double)rand() / RAND_MAX;
    }

    // Calcul des matrices de interpolation linéaire
    for (int i = 0; i < m; i++) {
        Xtilde[i] = fmod(X[i], dx) / dx;
        Ytilde[i] = fmod(Y[i], dy) / dy;
        I[i] = (int)(X[i] / dx) - 1;
        J[i] = (int)(Y[i] / dy) - 1;
    }

    
    for (int nn = 0; nn < m; nn++) {
        int i = I[nn];
        int j = J[nn];
        double xt = Xtilde[nn];
        double yt = Ytilde[nn];

        if (i != -1) {
            if (j != -1) {
                Gamma[nn * N + i + j * Nx] = xt * (yt - 1.0) + 1 - yt;
            }
            if (j != Ny - 1) {
                Gamma[nn * N + i + (j + 1) * Nx] = yt * (1.0 - xt);
            }
        }

        if (i != Nx - 1) {
            if (j != -1) {
                Gamma[nn * N + (i + 1) + j * Nx] = xt * (1.0 - yt);
            }
            if (j != Ny - 1) {
                Gamma[nn * N + (i + 1) + (j + 1) * Nx] = xt * yt;
            }
        }
    }
    for (int nn = 0; nn < m; nn++) {
        T[nn] = 0.0 ; 
        for (int k = 0; k < N; k++) {
            T[nn] += Gamma[nn * N + k] * U[k];
        }
        T[nn] += NOISE_SIGMA * ((double)rand() / RAND_MAX) + NOISE_MU;
    }

    // Libération de la mémoire
    free(Xtilde);
    free(Ytilde);
    free(I);
    free(J);
}
//Fonction pour générer un échantillon aléatoire sans remplacement
void generate_batch(int n, int batch_size, int *batch_indices) {
    // Initialisation du générateur de nombres aléatoires avec une graine
    srand(time(NULL));

    // Création de l'ensemble d'indices [0, 1, ..., 3*n-1]
    int all_indices[3 * n];

    for (int i = 0; i < 3 * n; i++) {
        all_indices[i] = i;
    }

    // Sélection d'un échantillon aléatoire sans remplacement
    for (int i = 0; i < batch_size; i++) {
        // Sélection d'un indice aléatoire
        int random_index = rand() % (3 * n - i);

        // Ajout de l'indice sélectionné au lot
        batch_indices[i] = all_indices[random_index];

        // Remplacement de l'indice sélectionné par le dernier indice non sélectionné
        all_indices[random_index] = all_indices[3 * n - i - 1];
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
    int m = 50 ;  
 
   //Vecteur solution à tn-1 et tn
    double* U0 = (double*) malloc(N* sizeof(double)); 
    double* Ue = (double*) malloc(N* sizeof(double));
    double* U = (double*) malloc(N* sizeof(double));
    double* gradI_T = (double *)malloc(m * sizeof(double));

   //vecteur 2nd membre
    double* b = (double*) malloc(N* sizeof(double));
    double* f = (double*) calloc(N, sizeof(double));

   //vecteur erreur 
    double* erreur = (double*) malloc(N* sizeof(double));

    //vecteur sol exacte
    double* U_exacte = (double*) malloc(N* sizeof(double));

   //Initialisation
    double dx = 1.0/(Nx+1) ; 
    double dy = 1.0/(Ny+1) ;
    double NOISE_SIGMA = 0.0 , NOISE_MU = 0.0 ; //noise

    for (int k =0 ; k<N ; k++){
      U0[k] = 0.0 ; 
      U[k] = 0.0 ; 
    }

   //Boucle principale
    int kmax = 10000;
    double eps = 1.0e-4 ; 
    double t = 0.0 ; 
    int n = 2;
    int ADAM = 1;  // 1 pour vrai, 0 pour faux 

    //Définition des pas pour approcher la dérivée
    double da, db, dw;
    double dalpha[3 * n];
    double best_alpha[3 * n];

    da = dx;
    db = dy;
    dw = dx;
    for (int i = 0; i < 3 * n; i++) {
        dalpha[i] = 0.1 * ((i % 3 == 0) ? da : ((i % 3 == 1) ? db : dw));   
    }
    double *x = (double *)malloc(Nx * sizeof(double));
    double *y = (double *)malloc(Ny * sizeof(double));
    double *ralpha = (double *)malloc(3 * n * sizeof(double));
    double *fe = (double *)malloc(N * sizeof(double));
    double *T = (double *)malloc(m * sizeof(double));
    double *samplesX = (double *)malloc(m * sizeof(double));
    double *samplesY = (double *)malloc(m * sizeof(double));
    double *Gamma = (double *)malloc(m * N * sizeof(double));
    for (int i = 0; i < Nx; i++) {
        x[i] = (i + 1) * dx;
    }

    for (int j = 0; j < Ny; j++) {
        y[j] = (j + 1) * dy;
    }

   // Génération du champ aléatoire
    for (int i = 0; i < 3 * n; i++) {
        ralpha[i] = (double)rand() / RAND_MAX;
    }

    // Calcul du champ exact
    compute_f_lagrange(n, ralpha, dx, dy, Nx, Ny, fe);
    BuildSecondMember( U0,b, cas,fe);  
    GC(b , Ue , eps, kmax) ; 

    // Échantillonnage du champ avec bruit
    get_random_samples(Ue, m, dx, dy, Nx, Ny, NOISE_MU, NOISE_SIGMA, T, samplesX, samplesY, Gamma);
    
    double histories[NB_RESTARTS][NB_ITER_MAX];
    double lowerDiff = 1e16;

    // Boucle pour les redémarrages
    for (int i = 0; i < NB_RESTARTS; i++) {
        double history[NB_ITER_MAX];
        double alpha[3 * n];
        for (int i = 0; i < 3*n; i++) {
            alpha[i] = ((double)rand() / RAND_MAX); // Nombre aléatoire entre 0 et 1
        }

        double normG = 1e16;
        double G2sum = 0.0;
        double adam_mw[3 * n];
        double adam_vw[3 * n] ;
        memset(adam_mw, 0, sizeof(adam_mw));
        memset(adam_vw, 0, sizeof(adam_vw));
        double adam_beta1 = 0.9;
        double adam_beta2 = 0.999;

        // Boucle pour la descente de gradient
        for (int nb_iter = 1; nb_iter <= NB_ITER_MAX; nb_iter++) {
            // Conditions d'arrêt
            if (normG < GRADIENT_EPSILON) {
                break;
            }
            if (nb_iter >= 100 && normG > GIVE_UP_THRESHOLD) {
                break;
            }
            printf("Avancement nb iter : %d \n",nb_iter) ; 
            // Calcul de f et U
            compute_f_lagrange(n, alpha, dx, dy, Nx, Ny, f);
            BuildSecondMember( U0,b, cas,f);  
            GC(b , U , eps, kmax) ; 
            for (int nn = 0; nn < m; nn++) {
                gradI_T[nn] = 0.0;
            }
            // Calculer gradI_T
            for (int nn = 0; nn < m; nn++) {
                for (int i = 0; i < N; i++) {
                    gradI_T[nn] += Gamma[nn * N + i] * U[nn];  
                }
            }
            // Calcul du gradient
            int batch[BATCH_SIZE];
            generate_batch(n, BATCH_SIZE, batch);
            printf("TEST \n") ; 

            double G[3 * n];

            memset(G, 0, sizeof(G));
            // Calcul de G
            for (int k_index = 0; k_index < BATCH_SIZE; k_index++) {
                int k = batch[k_index];
                double alphaL[3 * n], alphaR[3 * n];
                // Copie de alpha dans alphaL et alphaR
                memcpy(alphaL, alpha, 3 * n * sizeof(double));
                memcpy(alphaR, alpha, 3 * n * sizeof(double));
    
                // Modification de alphaL et alphaR
                alphaL[k] -= dalpha[k];
                alphaR[k] += dalpha[k];
    
                // Calcul de df = compute_f_lagrange(n, alphaR, dx, dy, Nx, Ny) - compute_f_lagrange(n, alphaL, dx, dy, Nx, Ny)
                double dfL[N], dfR[N],df[N];
                compute_f_lagrange(n, alphaR, dx, dy, Nx, Ny, dfR);
                compute_f_lagrange(n, alphaL, dx, dy, Nx, Ny, dfL);

                for (int j = 0; j < N; j++) {
                    df[j] = (dfR[j] - dfL[j]) / (2. * dalpha[k]);
                }

            // Résolution du système linéaire H * dU = df
                double dU[N];
                BuildSecondMember( U0,b, cas,df);  
                GC(b , dU , eps, kmax) ;

            // Calcul de G[k] = gradI_T @ (Gamma @ dU)
                G[k] = 0.0;
                for (int j = 0; j < m; j++) {
                    for (int i = 0; i < N; i++) {
                        G[k] += gradI_T[j] * Gamma[j * N + i] * dU[i];
                        }
                    }
            }

            // Mise à jour de l'historique
            normG = 0.0;
            for (int k = 0; k < 3 * n; k++) {
                normG += G[k] * G[k];
            }
            normG = sqrt(normG);
            printf("normG : %f \n ",normG) ; 
            history[nb_iter - 1] = normG;

            // Descente de gradient (Adam ou AdaGrad)
            if (ADAM==1) {
                for (int i = 0; i < 3 * n; i++) {
                    adam_mw[i] = adam_beta1 * adam_mw[i] + (1.0 - adam_beta1) * G[i];
                    double m_hat = adam_mw[i] / (1.0 - pow(adam_beta1, nb_iter));

                    adam_vw[i] = adam_beta2 * adam_vw[i] + (1.0 - adam_beta2) * pow(G[i], 2);
                    double v_hat = adam_vw[i] / (1.0 - pow(adam_beta2, nb_iter));

                    alpha[i] -= LEARNING_RATE * m_hat / (sqrt(v_hat) + ADAM_EPSILON);
                }
            }
             else {
                for (int i = 0; i < 3 * n; i++) {
                    G2sum += pow(G[i], 2);
                    alpha[i] -= LEARNING_RATE * G[i] / sqrt(G2sum);
                }            
            }
        }

        // Calcul du dernier état
        compute_f_lagrange(n, alpha, dx, dy, Nx, Ny, f);
        BuildSecondMember( U0,b, cas,f);  
        GC(b , U , eps, kmax) ; 

        // Comparaison avec la meilleure solution
        double finalDiff = 0.0;
        for (int i = 0; i < m; i++) {
            double diff = 0.0;
            for (int j = 0; j < N; j++) {
                diff += (Gamma[i * N + j] * U[j] - T[i]) * (Gamma[i * N + j] * U[j] - T[i]);
            }
            finalDiff += diff;
        }

        finalDiff = sqrt(finalDiff);

        // Sauvegarde si c'est la meilleure solution
        if (lowerDiff > finalDiff) {
            lowerDiff = finalDiff;
            memcpy(best_alpha, alpha, 3 * n * sizeof(double));
        }

        // Sauvegarde de l'historique
        memcpy(histories[i], history, NB_ITER_MAX * sizeof(double));

        // Conditions d'arrêt
        if (normG < GRADIENT_EPSILON) {
            break;
        }
    }





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
                fprintf(fp_exacte, "%f %f %f \n", (i+1)*dx, (j+1)*dy , Ue[j*Nx +i]);
            }
            else if (cas ==2){
                fprintf(fp_exacte, "%f %f %f \n", (i+1)*dx, (j+1)*dy ,  Ue[j*Nx+i]  );
                }
            }
        }

    // Fermeture du fichier
    fclose(fp_num);
    fclose(fp_exacte);
    free(U0);
    free(U);
    free(b) ;  
    free(U_exacte) ; 
    free(erreur) ; 
    free(x);
    free(y) ;
    free(ralpha);
    free(fe);
    free(T);
    free(samplesX);
    free(samplesY);
    free(Gamma);     
    free(f) ;
    free(gradI_T); 
   return 0;

}