double Source(double x, double y, int cas);
double UpBottomCondition(double x, double y, int cas);
double LeftRightCondition(double x, double y, int cas);

extern const int n;
typedef double (*vecteur2D)[2];
typedef double vecteur[2];

double *lagrange(vecteur2D a, vecteur2D P_a, vecteur X);