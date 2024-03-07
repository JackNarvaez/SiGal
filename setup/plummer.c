#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define RAND() ((double)rand()/(double)(RAND_MAX))

const double TPI        = 2*M_PI;
const double sqrt2      = sqrt(2.);
const double scaleftr   = 16.0/(3.0*M_PI);
const double invsclftr  = 1./scaleftr;
const double sqrtscldrt = sqrt(scaleftr);

void spher2cartes(double * Vec, double r, double sclftr) {
    double X1   = RAND();
    double X2   = RAND();
    double aux1 = (1.-2.*X1)*r;
    double aux2 = sqrt(r*r - aux1*aux1)*sclftr;
    Vec[0]  = aux2 * cos(TPI*X2);
    Vec[1]  = aux2 * sin(TPI*X2);
    Vec[2]  = aux1*sclftr;
}

double g(double x) {
    double x2 = x*x;
    return x2*pow(1-x2, 3.5);
}

double rndm(double min, double max) {
    return min + (max-min)*RAND();
}

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

void adjust_units(double *Pos, double*Vel, double *Mass, const int N) {
    double epot = 0.0;
    double ekin = 0.0;
    double xrel = 0.0;
    double yrel = 0.0;
    double zrel = 0.0;
    int ii, jj;
    for (ii=0; ii<N; ii++) {
        ekin += Mass[ii]*(Vel[3*ii]*Vel[3*ii] + Vel[3*ii+1]*Vel[3*ii+1] + Vel[3*ii+2]*Vel[3*ii+2]);
        for (jj = ii+1; jj<N; jj++) {
            xrel =  Pos[3*ii]-Pos[3*jj];
            yrel =  Pos[3*ii+1]-Pos[3*jj+1];
            zrel =  Pos[3*ii+2]-Pos[3*jj+2];
            epot += -Mass[jj]*Mass[ii] / sqrt(xrel*xrel+yrel*yrel+zrel*zrel);
        }
    }
    double alpha = -2*epot;
    double beta  = 0.5/sqrt(ekin);
    for (ii=0; ii<N; ii++) {
        Pos[3*ii]   *= alpha;
        Pos[3*ii+1] *= alpha;
        Pos[3*ii+2] *= alpha;
        Vel[3*ii]   *= beta;
        Vel[3*ii+1] *= beta;
        Vel[3*ii+2] *= beta;
    }
}

void plummer_dist(double *Pos, double *Vel, double *Mass, const int N) {
    double X1, X2;
    double cum_mass, r, Ve;
    double m = 1.0/N;
    double cum_mass_min = 0.0;
    double cum_mass_max = m;
    int ii;
    for (ii=0; ii<N; ii++) {
        Mass[ii] = m;
        // Position
        cum_mass = rndm(cum_mass_min, cum_mass_max);
        cum_mass_min = cum_mass_max;
        cum_mass_max += m;

        r    = 1./sqrt(pow(cum_mass, -2./3.)-1.);
        spher2cartes(Pos + 3*ii, r, invsclftr);

        // Velocity
        X1   = RAND();
        X2   = RAND();
        while (0.1*X2 < g(X1)) {
            X1 = RAND();
            X2 = RAND();
        }
        Ve = sqrt2 * pow(1.+r*r, -0.25)*X1;
        spher2cartes(Vel + 3*ii, Ve, sqrtscldrt);
    }
    frm2com(Pos, Vel, Mass, N);
    adjust_units(Pos, Vel, Mass, N); 
}