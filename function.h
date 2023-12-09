typedef double* vecteur;
typedef vecteur* vecteur2D;

double Source(double x, double y, int cas);
double UpBottomCondition(double x, double y, int cas);
double LeftRightCondition(double x, double y, int cas);

void LagrangeCoeff(vecteur2D pt_x, vecteur2D pt_y, vecteur x, vecteur y); //Prend en entrée les points 2D d'interpolation pt_x, les valeurs 2D à matcher pt_y et un point x=(x1,x2) où évaluer ce polynôme et met la valeur du polynôme en ce point dans y=(y1,y2)
