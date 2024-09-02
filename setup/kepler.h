#ifndef KEPLER_H_
#define KEPLER_H_

double uniform( double );
double keplerian_velocity(double );
void cilin2cartes(double , double , double , double *, double *);
void keplerian_disk(double *, double *, double *, const int , const double, const double , const double);

#endif // KEPLER_H_
