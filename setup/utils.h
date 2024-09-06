// utils.h
#ifndef UTILS_H_
#define UTILS_H_

#define RAND() ((double)rand()/(double)(RAND_MAX))
using F_ROOT = double(double, double);

// Prototipos de funciones comunes
double rand_normal(double, double);
double rndm(double, double);
void spher2cartes(double *, double);
void frm2com(double *, double *, double *, const int);
double ZerosBisection(double, int *, F_ROOT);
double epicyclic_frequency(double, double);
double radial_velocity_dispersion(double, double, double);
double tangential_velocity_dispersion(double, double, double);
double vertical_velocity_dispersion(double, double, double );

#endif // UTILS_H_
