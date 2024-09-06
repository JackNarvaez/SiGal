// --------------------------------------------------------------- //
// The Plummer's model. Plummer (1911)                             //
// --------------------------------------------------------------- //

#include <math.h>
#include "utils.h"

const double TPI        = 2*M_PI;
const double sqrt2      = sqrt(2.);

//---------------------------------------------------------------------- //
// Probability distribution of velocity.                                 //
//---------------------------------------------------------------------- //
double g(double x) {
    double x2 = x*x;
    return x2*pow(1-x2, 3.5);
}

void plummer_dist(double *Pos, double *Vel, double *Mass, const int Nl, const int N, const double R, const double M) {
    double X1, X2, X3;
    double r, Ve;
    double dm = M/N;
    int ii;
    for (ii=0; ii<Nl; ii++) {
        Mass[ii] = dm;
        // Position
        do {
            X1 = RAND();
            r = 1./sqrt(pow(X1, -2./3.)-1.);
        } while (r > R);

        spher2cartes(Pos + 3*ii, R);

        // Velocity
        do {
            X2 = RAND();
            X3 = RAND();
        } while (0.1*X3 < g(X2));
        r *= 1./R;
        Ve = sqrt2*sqrt(M/R)*pow(1.+r*r, -0.25)*X2;
        spher2cartes(Vel + 3*ii, Ve);
    }
    frm2com(Pos, Vel, Mass, Nl);
}