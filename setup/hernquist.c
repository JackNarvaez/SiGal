// --------------------------------------------------------------- //
// The Hernquist’s model. Hernquist (1990)                         //
// --------------------------------------------------------------- //

#include <math.h>
#include "utils.h"

const double TPI        = 2*M_PI;
const double sqrt2      = sqrt(2.);


//---------------------------------------------------------------------- //
// Probability distribution of velocity.                                 //
//---------------------------------------------------------------------- //
double h_g(double x) {
    double q2   = 1.-x*x;
    double q    = sqrt(q2);
    return x*x / pow(1.0 - q2, 2.5) * (3.0 * asin(q) + q * sqrt(1.0 - q2) * (1.0 - 2.0 * q2) * (8.0 * q2 * q2 - 8.0 * q2 - 3.0));
}

double hernquist_density(double r, double M, double a) {
    return (M * a) / (r * pow(r + a, 3));
}

void hernquist_dist(double *Pos, double *Vel, double *Mass, const int Nl, const int N, const double R, const double M) {
    double X1, X2, X3;
    double r, Ve;
    double dm = M / N;
    int ii;
    for (ii = 0; ii < Nl; ii++) {
        Mass[ii] = dm;
    
        do {
            X1 = RAND();
            r = 1./ (1./sqrt(X1) - 1.);
        } while (r > R);

        spher2cartes(Pos + 3 * ii, r);

        // Velocity [Warning]!!
        do {
            X2 = RAND();
            X3 = RAND();
        } while (X3 > h_g(X2));
        r *= 1./R;
        Ve =  sqrt2*sqrt(M)*sqrt(1.0 / (1.0 +  r))*X2;
        spher2cartes(Vel + 3 * ii, Ve);
    }

    frm2com(Pos, Vel, Mass, Nl);
}
