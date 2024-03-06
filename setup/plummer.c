#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define RAND() ((double)rand()/(double)(RAND_MAX))

const double TPI = 2*M_PI;
const double sqrt2 = sqrt(2.);

void spher2cartes(double * Vec, double r, double X2, double X3) {
    Vec[2]  = (1.-2.*X2)*r;
    double aux = pow(r*r - Vec[2]*Vec[2], 0.5);
    Vec[0]  = aux * cos(TPI*X3);
    Vec[1]  = aux * sin(TPI*X3);
}

double g(double x) {
    double x2 = x*x;
    return x2*pow(1-x2, 3.5);
}

void plummer_dist(double *Pos, double *Vel, double *Mass, const int N) {
    double X1, X2, X3, X4, X5, X6, X7;
    double r, Ve;
    double m = 1./N;
    int ii;
    for (ii=0; ii<N; ii++) {
        Mass[ii] = m;
        X1   = RAND();
        X2   = RAND();
        X3   = RAND();
        // Position
        r    = pow(pow(X1, -2./3.)-1., -0.5);
        spher2cartes(Pos + 3*ii, r, X2, X3);

        // Velocity
        X4   = RAND();
        X5   = RAND();
        while (0.1*X5 < g(X4)) {
            X4 = RAND();
            X5 = RAND();
        }
        Ve = sqrt2 * pow(1.+r*r, -0.25)*X4;
        X6   = RAND();
        X7   = RAND();
        spher2cartes(Vel + 3*ii, Ve, X6, X7);
    }
}