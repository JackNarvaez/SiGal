#ifndef KUZMIN_H_
#define KUZMIN_H_

double rand_uniform(double, double);
double rand_normal(double, double);
double kuzmin_velocity(double);
double surface_density(double);
void frm2com(double *, double *, double *, const int);
double sample_kuzmin_radius(double );
void kuzmin_disk(double *, double *, double *, const int, const double, const double, const double);

#endif // KUZMIN_H_