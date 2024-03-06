#ifndef NBODIES_H_
#define NBODIES_H_

typedef void (* function)(const double *, const double *, double *, const int *, const int, const int, const int, const int, const int, MPI_Status);
typedef void (* Integrator)(double *, double *, const double *, double *, double *, double *, const double, const int, function, const int *, const int, const int, const int, const int, MPI_Status);

void Gravitational_Acc(double *, const double *, const double *, const double *, const int *, const int);

void Acceleration(const double *, const double *, double *, const int *, const int, const int, const int, const int, const int, MPI_Status);

void Save_vec(FILE *, const double *, const double *, const int);

void Save_data(const double *, const double *, const int *, const int, const int, const int, const int, const int, MPI_Status);

void Euler(double *, double *, const double *, double *, double *, double *, const double, const int, function, const int *, const int, const int, const int, const int, MPI_Status);

void PEFRL(double *, double *, const double *, double *, double *, double *, const double, const int, function, const int *, const int, const int, const int, const int, MPI_Status);

void Evolution(double *, double *, const double *, double *, double *, double *, const int *, const int, const int, const int, const int, const int, MPI_Status, const int, const double, const int, function, Integrator);

#endif // NBODIES_H_