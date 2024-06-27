#ifndef EXPONENTIAL_H_
#define EXPONENTIAL_H_

double rand_uniform(double, double);
double exponential_velocity(double);
void frm2com(double *, double *, double *, const int);
void exponential_disk(double *, double *, double *, const int, const double, const double, const double);
double rand_normal(double, double);
double epicyclic_frequency(double);
double radial_velocity_dispersion(double);
double tangential_velocity_dispersion(double);
double vertical_velocity_dispersion(double);
double sample_exponential_radius(double );
double surface_density_exponential(double);

#endif // EXPONENTIAL_H_