#ifndef HERNQUIST_H_
#define HERNQUIST_H_

double rndm(double, double);
void spher2cartes(double *, double, double);
double hernquist_density(double , double , double );
double g(double);
void frm2com(double *, double *, double *, const int);
void hernquist_dist(double *, double *, double *, const int, const double, const double, const double);

#endif // HERNQUIST_H_