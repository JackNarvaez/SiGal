#ifndef KEPLER_H_
#define KEPLER_H_

double rndm(double, double);
double uniform( double );
double keplerian_velocity(double );
void cilin2cartes(double , double , double , double *, double *);
void frm2com(double *, double *, double *, const int);
void keplerian_disk(double *, double *, double *, const int , const double, const double , const double);

#endif // KEPLER_H_
