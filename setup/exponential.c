#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "utils.h"  

#define G 1.0
#define M_c 1.0         // Masa del cuerpo central
#define Q 1.5           // Parámetro de Toomre

double exponential_velocity(double r) {
    return sqrt(G * M_c * r * exp(-r / 1.0));
}

double surface_density_exponential(double r) {
    return (M_c / (2 * M_PI )) * exp(-r / 1.0);
}

double sample_exponential_radius(double R) {
    double r; 
    do {
        r = - 1.0* log(rand_uniform(0, 1)); // Método de muestreo por transformada inversa
    } while (r > 5*R);
    return r;
}

void exponential_disk(double *Pos, double *Vel, double *Mass, const int N, const double M, const double R, const double seed) {
    srand(seed);
    double r, theta, sigma_rad, sigma_tan, sigma_z;

    for (int ii = 0; ii < N; ii++) {
        theta = rand_uniform(0, 2 * M_PI);
        r = sample_exponential_radius(R);  
        double vel = exponential_velocity(r);
        double sigma = surface_density_exponential(r);

        Pos[3 * ii] = r * cos(theta);
        Pos[3 * ii + 1] = r * sin(theta);
        Pos[3 * ii + 2] = 0.0;

        sigma_rad = radial_velocity_dispersion(r, vel, sigma);
        sigma_tan = tangential_velocity_dispersion(r, vel, sigma);
        sigma_z   = vertical_velocity_dispersion(r, vel, sigma);

        double vr = rand_normal(0, sigma_rad);
        double vtheta = exponential_velocity(r) + rand_normal(0, sigma_tan);
        double vz = rand_normal(0, sigma_z);

        Vel[3 * ii] = vr * cos(theta) - vtheta * sin(theta);
        Vel[3 * ii + 1] = vr * sin(theta) + vtheta * cos(theta);
        Vel[3 * ii + 2] = vz;

        Mass[ii] = M / N;
    }
    frm2com(Pos, Vel, Mass, N);
}