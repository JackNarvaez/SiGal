#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define RAND() ((double)rand()/(double)(RAND_MAX))

const double TPI        = 2*M_PI;
const double sqrt2      = sqrt(2.);

//---------------------------------------------------------------------- //
// Generate a random number in [min, max] with uniform distribution.     //
//---------------------------------------------------------------------- //
double rndm(double min, double max) {
    return min + (max-min)*RAND();
}

//---------------------------------------------------------------------- //
// Assign a vector with magnitude 'r' a random direction.                //
//---------------------------------------------------------------------- //
void spher2cartes(double * Vec, double r) {
    double X1   = RAND();
    double X2   = RAND();
    double aux1 = (1.-2.*X1)*r;
    double aux2 = sqrt(r*r - aux1*aux1);
    Vec[0]  = aux2 * cos(TPI*X2);
    Vec[1]  = aux2 * sin(TPI*X2);
    Vec[2]  = aux1;
}

//---------------------------------------------------------------------- //
// Probability distribution of velocity.                                 //
//---------------------------------------------------------------------- //
double g(double x) {
    double x2 = x*x;
    return x2*pow(1-x2, 3.5);
}

//---------------------------------------------------------------------- //
// Convert to the CoM reference frame                                    //
//---------------------------------------------------------------------- //
void frm2com(double *Pos, double *Vel, double *Mass, const int N) {
    double x_com    = 0.0;
    double y_com    = 0.0;
    double z_com    = 0.0;
    double vx_com   = 0.0;
    double vy_com   = 0.0;
    double vz_com   = 0.0;
    int ii;
    for (ii=0; ii<N; ii++) {
        x_com   += Mass[ii]*Pos[3*ii];
        y_com   += Mass[ii]*Pos[3*ii+1];
        z_com   += Mass[ii]*Pos[3*ii+2];
        vx_com  += Mass[ii]*Vel[3*ii];
        vy_com  += Mass[ii]*Vel[3*ii+1];
        vz_com  += Mass[ii]*Vel[3*ii+2];
    }
    for (ii=0; ii<N; ii++) {
        Pos[3*ii]   -= x_com;
        Pos[3*ii+1] -= y_com;
        Pos[3*ii+2] -= z_com;
        Vel[3*ii]   -= vx_com;
        Vel[3*ii+1] -= vy_com;
        Vel[3*ii+2] -= vz_com;
    }
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