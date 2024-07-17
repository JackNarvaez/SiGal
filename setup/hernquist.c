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
// Probability distribution of velocity.
double g(double q) {
    double q2 = q * q;
    double E = q2;
    double sqrE = q;
    return (1.0 / pow(1.0 - E, 5.0 / 2.0) * (3.0 * asin(sqrE) + sqrE * sqrt(1.0 - E) * (1.0 - 2.0 * E) * (8.0 * E * E - 8.0 * E - 3.0)));
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
// Generate a Hernquist sphere of mass M and parametrized radius R.        //
// rmax=5R, which is 94% of the total mass for the Hernquist distribution. //
//---------------------------------------------------------------------- //

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
        } while (X3 > g(X2));

        //q = sqrt(-r * G * M / (r + R));
        Ve =  sqrt2*sqrt(M/R)* sqrt(1.0 / (1.0 +  r));
        spher2cartes(Vel + 3 * ii, Ve);
    }

    frm2com(Pos, Vel, Mass, N);
}
