#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#define G 1.0     
#define M_central 1.0
#define RAND() ((double)rand()/(double)(RAND_MAX))

//---------------------------------------------------------------------- //
// Generate a random number in [min, max] with uniform distribution.     //
//---------------------------------------------------------------------- //
double rndm(double min, double max) {
    return min + (max-min)*RAND();
}
 
double keplerian_velocity(double r) {
    return sqrt(G * M_central / r);
}

//---------------------------------------------------------------------- //
// Probability distribution of velocity.                                 //
//---------------------------------------------------------------------- //
double uniform(double r_max) {
    return r_max * sqrt(rndm(0.0, 1.0));
}

//---------------------------------------------------------------------- //
// Assign vel ond pos with magnitude 'r' a random direction.                //
//---------------------------------------------------------------------- //
void cilin2cartes(double r, double theta, double v, double *Pos, double *Vel) {
    Pos[0] = r * cos(theta);      // x
    Pos[1] = r * sin(theta);      // y
    Pos[2] = 1e-3*rndm(0.1,0.2);                 // z;

    Vel[0] = -v * sin(theta) + 0.01 * v * (RAND() - 0.5);   // vx with dispersion
    Vel[1] = v * cos(theta) + 0.01 * v * (RAND() - 0.5);    // vy with dispersion
    Vel[2] = 1e-3*rndm(0.1,0.2);                                           // vz
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
// Generate a uniform keplerian disk sphere of mass M and parametrized radius R.        //
//---------------------------------------------------------------------- //
void keplerian_disk(double *Pos, double *Vel, double *Mass, const int N, const double R, const double M, const double seed) {
    srand(seed);
    double m = M / N;
    double r_max = R;
    int ii;
    for (ii = 0; ii < N; ii++) {
        Mass[ii] = m;

        double r = uniform(r_max);
        double theta = rndm(0, 2.0 * M_PI);
        double v = keplerian_velocity(r);
        cilin2cartes(r, theta, v, Pos + 3*ii, Vel + 3*ii);
    }
    frm2com(Pos, Vel, Mass, N);

}