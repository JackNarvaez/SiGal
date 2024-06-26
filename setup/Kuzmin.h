#ifndef KUZMIN_H_
#define KUZMIN_H_

double rand_uniform(double, double);
double kuzmin_velocity(double);
void frm2com(double *, double *, double *, const int);
void kuzmin_disk(double *, double *, double *, const int, const double, const double, const double);
double rand_normal(double, double);
double epicyclic_frequency(double);
double sample_kuzmin_radius(double );
double radial_velocity_dispersion( double);
double tangential_velocity_dispersion( double);
double vertical_velocity_dispersion( double);
double surface_density(double);

#endif // KUZMIN_H_