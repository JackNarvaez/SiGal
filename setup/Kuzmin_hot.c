#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#define G 1.0
#define M_c 1.0         // Mass central body
#define Q 1.5           // Toomre Parameter
#include "Kuzmin.h"
#include "utils.h"



void kuzmin_hot_disk(double *Pos, double *Vel, double *Mass, const int N, const double R, const double M, const double seed) {
    srand(seed);
    double r, theta, sigma_rad, sigma_z, sigma_tan, vel, sigma;

    for (int ii = 0; ii < N; ii++) {
        theta = rand_uniform(0, 2 * M_PI);
        r = sample_kuzmin_radius(R); 
        vel = kuzmin_velocity(r);
        sigma = surface_density(r);

        Pos[3 * ii] = r * cos(theta);
        Pos[3 * ii + 1] = r * sin(theta);
        Pos[3 * ii + 2] = 0.0;

        sigma_rad = radial_velocity_dispersion(r, vel, sigma);
        sigma_tan = tangential_velocity_dispersion(r, vel,  sigma);
        sigma_z   = vertical_velocity_dispersion(r, vel, sigma);

        double vr = rand_normal(0, sigma_rad);
        double vtheta = rand_normal(kuzmin_velocity(r), sigma_tan);
        //double vtheta = rand_uniform(0, kuzmin_velocity(r)) ;
        double vz = rand_normal(0, sigma_z);

        Vel[3 * ii] = vr * cos(theta) - vtheta * sin(theta);
        Vel[3 * ii + 1] = vr * sin(theta) + vtheta * cos(theta);
        Vel[3 * ii + 2] = vz;

        Mass[ii] = M / N;
    }
    frm2com(Pos, Vel, Mass, N);
}