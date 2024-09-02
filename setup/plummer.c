#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
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

//---------------------------------------------------------------------- //
// Generate a Plummer sphere of mass M and parametrized radius R.        //
// rmax=5R, which is 94% of the total mass for the Plummer distribution. //
//---------------------------------------------------------------------- //
void plummer_dist(double *Pos, double *Vel, double *Mass, const int N, const double R, const double M, const double seed) {
    srand(seed);
    double X1, X2, X3;
    double r, Ve;
    double m = 1.0/N;
    double dm = m*M;
    double cum_mass_min = 0.0;
    double cum_mass_max = m;
    int ii;
    for (ii=0; ii<N; ii++) {
        Mass[ii] = dm;
        // Position
        do {
            do {
                X1 = rndm(cum_mass_min, cum_mass_max);
                cum_mass_min = cum_mass_max;
                cum_mass_max += m;
                X1 = RAND();
            } while (X1 < 1.0e-10);
            r = 1./sqrt(pow(X1, -2./3.)-1.);
        } while (r > 5);

        r *= R;
        spher2cartes(Pos + 3*ii, r);

        // Velocity
        X2   = RAND();
        X3   = RAND();
        while (0.1*X3 < g(X2)) {
            X2 = RAND();
            X3 = RAND();
        }
        Ve = sqrt(M/R)*sqrt2 * pow(1.+r*r, -0.25)*X2;
        spher2cartes(Vel + 3*ii, Ve);
    }
    frm2com(Pos, Vel, Mass, N);
}