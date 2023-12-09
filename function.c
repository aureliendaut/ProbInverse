#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "parametres.h"

//contient le terme source
double Source(double x, double y, int cas){
        if (cas == 1)
    {
        return 2.0*(y-(y*y) +x - (x*x));  ; 
    }
    
    else if(cas==2)
  {
        return sin(x) + cos(y);
  }

    else if(cas==3)//Instationnaire
  {
      return exp(-pow((x-Lx/2),2.0)) * exp(-pow((y-Ly/2.0),2)) * cos(acos(-1.)/2.0);
  }

  else
  {
        return 0.0;
  }
}
//Condition bord haut/bas
double UpBottomCondition(double x, double y, int cas){ 
        if (cas == 1)
    {
            return 0.0 ; 
    }
    
    else if(cas==2)
  {
        return sin(x) + cos(y);

  }

    else if(cas==3)//Instationnaire
  {
        return 0.0;

  }

    else
  {
        return 0.0;
  }
}
//condition de bord gauche/droit
double LeftRightCondition(double x, double y, int cas){
    if (cas == 1)
    {
        return 0.0 ; 
    }
    
    else if(cas==2)
  {
        return sin(x) + cos(y);
  }

    else if(cas==3)//Instationnaire
  {
        return 1.0;
  }

  else
  {
        return 0.0;
  }
}


typedef double* vecteur;
typedef vecteur* vecteur2D;

//Prend en entrée les points 2D d'interpolation pt_x, les valeurs 2D à matcher pt_y et un point x=(x1,x2) où évaluer ce polynôme et met la valeur du polynôme en ce point dans y=(y1,y2), n le nombre de points d'interpolation
void LagrangeCoeff(vecteur2D pt_x, vecteur2D pt_y, vecteur x, vecteur y, int n)
{
	//Allocations
	/*x = (double*)malloc(2*sizeof(double));
	y = (double*)malloc(2*sizeof(double));
	pt_x = (vecteur*)malloc(n*sizeof(vecteur));
	pt_y = (vecteur*)malloc(n*sizeof(vecteur));
	for (int i = 0; i < n; i++)
	{
        	pt_x[i] = (double*)malloc(2*sizeof(double));
        	pt_y[i] = (double*)malloc(2*sizeof(double));
    	}

	free(x);
	free(y);*/
	
	y[0]=0;
	y[1]=0;
	for(int i = 0 ; i < n ; i++)
	{
		for(int j = 0 ; j < n ; j++)
		{
			double Lx=1;
		        double Ly=1;
			for(int k = 0 ; k < n ; k++)
			{
				if( k == i && k == j)
				{
				}
				else if(k == i)
				{
					Ly *= (x[1] - pt_x[k][1])/(pt_x[j][1] - pt_x[k][1]);
				}
				else if(k == j)
				{
					Lx *= (x[0] - pt_x[k][0])/(pt_x[i][0] - pt_x[k][0]);
				}
				else
				{
					Lx *= (x[0] - pt_x[k][0])/(pt_x[i][0] - pt_x[k][0]);
					Ly *= (x[1] - pt_x[k][1])/(pt_x[j][1] - pt_x[k][1]);
				}
			}
			y[0] += pt_y[i][0] * Lx * Ly;
			y[1] += pt_y[j][1] * Lx * Ly;
		}
	}

}
