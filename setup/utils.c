#include <math.h>
#include "utils.h"

#define Q   1.5     // Toomre Parameter for Kuzmin galaxy
#define TPI (2 * M_PI)
#define SQRT2 sqrt(2)

using F_ROOT = double(double, double);

double rndm(double min, double max) {
    return min + (max - min) * RAND();
}

double rand_normal(double mean, double stddev) {
    double u1 = RAND();
    double u2 = RAND();
    double z0 = sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
    return z0 * stddev + mean;
}

void spher2cartes(double *Vec, double r) {
    double X1 = RAND();
    double X2 = RAND();
    double aux1 = (1.0 - 2.0 * X1) * r;
    double aux2 = sqrt(r * r - aux1 * aux1);
    Vec[0] = aux2 * cos(TPI * X2);
    Vec[1] = aux2 * sin(TPI * X2);
    Vec[2] = aux1;
}

void frm2com(double *Pos, double *Vel, double *Mass, const int N) {
    double x_com = 0.0, y_com = 0.0, z_com = 0.0;
    double vx_com = 0.0, vy_com = 0.0, vz_com = 0.0;
    double total_mass = 0.0;
    for (int ii = 0; ii < N; ii++) {
        total_mass += Mass[ii];
        x_com += Mass[ii] * Pos[3 * ii];
        y_com += Mass[ii] * Pos[3 * ii + 1];
        z_com += Mass[ii] * Pos[3 * ii + 2];
        vx_com += Mass[ii] * Vel[3 * ii];
        vy_com += Mass[ii] * Vel[3 * ii + 1];
        vz_com += Mass[ii] * Vel[3 * ii + 2];
    }

    for (int ii = 0; ii < N; ii++) {
        Pos[3 * ii] -= x_com / total_mass;
        Pos[3 * ii + 1] -= y_com / total_mass;
        Pos[3 * ii + 2] -= z_com / total_mass;
        Vel[3 * ii] -= vx_com / total_mass;
        Vel[3 * ii + 1] -= vy_com / total_mass;
        Vel[3 * ii + 2] -= vz_com / total_mass;
    }
}

double ZerosBisection(double x, int *success, F_ROOT f_root) {
    double it_max = 100;
    double err = 1.e-12;
    double a = 0, b = 100;
    double m, fa, fm;
    fa = f_root(a, x);
    int it = 0;
    while(b-a>= err && it < it_max) {
        m = 0.5*(a+b);
        fm = f_root(m, x);
	    if (fm*fa<0.) {
            b = m;
        } else {
            a = m;
            fa = fm;
        }
        it ++;
    }
    if (it ==100) {
        *success = 0;
    }
    return 0.5*(a+b);
}

double epicyclic_frequency(double r, double vel) {
    return SQRT2 * vel / r;
}

double radial_velocity_dispersion(double r,double vel,double sigma) {
    double kappa = epicyclic_frequency(r, vel);
    return (Q * 3.36  * sigma) / kappa;
}

double tangential_velocity_dispersion(double r, double vel,double sigma) {
    return (3.36 * sigma) / (SQRT2 * epicyclic_frequency(r, vel));
}

double vertical_velocity_dispersion(double r, double vel,double sigma) {
    return 0.5 * tangential_velocity_dispersion(r, vel, sigma);
}