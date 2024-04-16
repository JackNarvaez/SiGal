#ifndef HERNQUIST_H_
#define HERNQUIST_H_

void spher2cartes(double *, double, double);
double g(double);
double rndm(double, double);
void frm2com(double *, double *, double *, const int);
void adjust_units(double *, double*, double *, const int);
void plummer_dist(double *, double *, double *, const int);
void hernquist_dist(double *, double *, double *, const int);

#endif // HERNQUIST_H_