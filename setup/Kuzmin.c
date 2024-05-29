#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#define G 1.0
#define M_c 1.0         // Masa del cuerpo central
#define Q 1.5           // Parámetro de Toomre

double rand_uniform(double a, double b) {
    return a + (b - a) * rand() / (double)RAND_MAX;
}

double rand_normal(double mean, double stddev) {
    double u1 = rand_uniform(0, 1);
    double u2 = rand_uniform(0, 1);
    double z0 = sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
    return z0 * stddev + mean;
}

double kuzmin_velocity(double r) {
    return sqrt(G * M_c / r);
}

double epicyclic_frequency(double r) {
    return sqrt(2) * kuzmin_velocity(r) / r;
}

double surface_density(double r) {
    return (M_c / (2 * M_PI)) * pow(1 + r*r, -1.5);
}

double radial_velocity_dispersion(double r) {
    double kappa = epicyclic_frequency(r);
    double sigma = surface_density(r);
    return (Q * 3.36 * G * sigma) / kappa;
}

double tangential_velocity_dispersion(double r) {
    double sigma = surface_density(r);
    return (3.36 * G * M_c* sigma ) / (2 * r * epicyclic_frequency(r));
}

double vertical_velocity_dispersion(double r) {
    return 0.5 * tangential_velocity_dispersion(r);
}

//---------------------------------------------------------------------- //
// Convert to the CoM reference frame                                    //
//---------------------------------------------------------------------- //
void frm2com(double *Pos, double *Vel, double *Mass, const int N) {
    double x_com = 0.0;
    double y_com = 0.0;
    double z_com = 0.0;
    double vx_com = 0.0;
    double vy_com = 0.0;
    double vz_com = 0.0;
    int ii;
    double total_mass = 0.0;
    for (ii = 0; ii < N; ii++) {
        total_mass += Mass[ii];
        x_com += Mass[ii] * Pos[3 * ii];
        y_com += Mass[ii] * Pos[3 * ii + 1];
        z_com += Mass[ii] * Pos[3 * ii + 2];
        vx_com += Mass[ii] * Vel[3 * ii];
        vy_com += Mass[ii] * Vel[3 * ii + 1];
        vz_com += Mass[ii] * Vel[3 * ii + 2];
    }

    for (ii = 0; ii < N; ii++) {
        Pos[3 * ii] -= x_com;
        Pos[3 * ii + 1] -= y_com;
        Pos[3 * ii + 2] -= z_com;
        Vel[3 * ii] -= vx_com;
        Vel[3 * ii + 1] -= vy_com;
        Vel[3 * ii + 2] -= vz_com;
    }
}

// Función para inicializar un disco de Kuzmin con dispersión de velocidades
void kuzmin_disk(double *Pos, double *Vel, double *Mass, const int N, const double R, const double M, const double seed) {
    srand(seed);
    double r, theta, sigma_rad, sigma_tan, sigma_z;

    for (int ii = 0; ii < N; ii++) {
        theta = rand_uniform(0, 2 * M_PI);
        r = 1.0 * sqrt(rand_uniform(0, 1) / (1 - rand_uniform(0, 1)));  // Inverse transform sampling method

        Pos[3 * ii] = r * cos(theta);
        Pos[3 * ii + 1] = r * sin(theta);
        Pos[3 * ii + 2] = 0.0;

        sigma_rad = radial_velocity_dispersion(r);
        sigma_tan = tangential_velocity_dispersion(r);
        sigma_z = vertical_velocity_dispersion(r);

        double vr = rand_normal(0, sigma_rad);
        double vtheta = kuzmin_velocity(r) + rand_normal(0, sigma_tan);
        double vz = rand_normal(0, sigma_z);

        Vel[3 * ii] = vr * cos(theta) - vtheta * sin(theta);
        Vel[3 * ii + 1] = vr * sin(theta) + vtheta * cos(theta);
        Vel[3 * ii + 2] = vz;

        Mass[ii] = 1.0 / N;
    }
    frm2com(Pos, Vel, Mass, N);

}
