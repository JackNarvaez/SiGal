#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "utils.h"
#include <mpi.h>

#define G 1.0
#define M_c 1.0         // Mass central body
#define Q 1.5           // Toomre Parameter


double kuzmin_velocity(double r) {
    return sqrt( (G * M_c * r * r) * pow(r * r + 1.0, - 3.0 / 2.0) );
}

double surface_density(double r) {
    return (M_c / (2 * M_PI)) * pow(1 + r * r, -1.5);
}

double sample_kuzmin_radius(double R) {
    double r, u, density, density_max;
    do {
        r = rand_uniform(0, 15 * R);
        u = rand_uniform(0, 1);
        density = surface_density(r);
        density_max = surface_density(0); // Densidad máxima en r = 0
    } while (u > density / density_max);
    return r;
}

void kuzmin_disk(double *Pos, double *Vel, double *Mass, const int N, const double R, const double M, const double seed) {
    srand(seed);
    double r, theta;

    for (int ii = 0; ii < N; ii++) {
        theta = rand_uniform(0, 2.0 * M_PI);
        r = sample_kuzmin_radius(R); // Use inverse CDF sampling for radius
        
        Pos[3 * ii] = r * cos(theta);
        Pos[3 * ii + 1] = r * sin(theta);
        Pos[3 * ii + 2] = 0.0;


        double vr = 0.0;
        double vtheta = rand_uniform(0, kuzmin_velocity(r));
        double vz = rand_uniform(-0.1, 0.1);

        Vel[3 * ii] = vr * cos(theta) - vtheta * sin(theta);
        Vel[3 * ii + 1] = vr * sin(theta) + vtheta * cos(theta);
        Vel[3 * ii + 2] = vz;

        Mass[ii] = M / N;
    }
    frm2com(Pos, Vel, Mass, N);
}