#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "utils.h"

const double TPI        = 2*M_PI;
const double sqrt2      = sqrt(2.);


// Probability distribution of velocity.
double h_g(double q) {
    double q2 = q * q;
    double E = q2;
    double sqrE = q;
    return (1.0 / pow(1.0 - E, 5.0 / 2.0) * (3.0 * asin(sqrE) + sqrE * sqrt(1.0 - E) * (1.0 - 2.0 * E) * (8.0 * E * E - 8.0 * E - 3.0)));
}


double hernquist_density(double r, double M, double a) {
    return (M * a) / (r * pow(r + a, 3));
}


// Generate a Hernquist sphere of mass M and parametrized radius R.
// rmax = 5R, which is 94% of the total mass for the Hernquist distribution.
void hernquist_dist(double *Pos, double *Vel, double *Mass, const int N, const double R, const double M, const double seed) {
    srand(seed);
    double X1, X2, X3;
    double r, Ve;
    double m = 1.0 / N;
    double dm = m * M;
    int ii;
    for (ii = 0; ii < N; ii++) {
        Mass[ii] = dm;
    
        // Position: Acceptance-Rejection Method
        do {
            X1 = RAND();
            X2 = RAND();
            r = R * pow(X1, 1.0 / 3.0) / (1.0 - pow(X1, 1.0 / 3.0));  // Transformation inverse
        } while (X2 > hernquist_density(r, M, R) );  // Acceptance criterion
        spher2cartes(Pos + 3 * ii, r);

        // Velocity
        do {
            X2 = RAND();
            X3 = RAND();
        } while (X3 > h_g(X2));

        //q = sqrt(-r * G * M / (r + R));
        Ve =  sqrt2*sqrt(M/R)* sqrt(1.0 / (1.0 +  r));
        spher2cartes(Vel + 3 * ii, Ve);
    }

    frm2com(Pos, Vel, Mass, N);
}
