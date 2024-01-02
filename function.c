#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "function.h"
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

//Partie Lagrange

const int n = 10;

double *lagrange(vecteur2D a, vecteur2D P_a, vecteur X)
{
      //Les points a_i d'interpolation
      //Les valeurs en les points d'interpolation P(a_i)
      //Le vecteur X des points d'étude
      //Le vecteur résultat res
      //n le nombre de pts d'interpolation

      vecteur res;

	res[0]=0;
	res[1]=0;
      

	for(int i = 0 ; i < n ; i++)
	{
		for(int j = 0 ; j < n ; j++)
		{
			double L_i = 1;
		      double L_j = 1;

			for(int k = 0 ; k < n ; k++)
			{
				if( k == i && k == j)
				{
                              L_i = 1;
		                  L_j = 1;
				}
				else if(k == i)
				{
                              L_i = 1;
				      L_j = L_j*(X[1] - a[k][1])/(a[j][1] - a[k][1]);
				}
				else if(k == j)
				{
					L_i = L_i*(X[0] - a[k][0])/(a[i][0] - a[k][0]);
                              L_j = 1;
				}
				else
				{
					L_i = L_i*(X[0] - a[k][0])/(a[i][0] - a[k][0]);
					L_j = L_j*(X[1] - a[k][1])/(a[j][1] - a[k][1]);
				}
			}
			res[0] += P_a[i][0] * L_i * L_j;
			res[1] += P_a[j][1] * L_i * L_j;
		}
	}

      return res;

}