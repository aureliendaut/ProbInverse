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