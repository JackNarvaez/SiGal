#ifndef NBODIES_H_
#define NBODIES_H_
#include "structures.h"
#include "tree.h"

bool FarDistance(const double*, const Node*);
void Forcei(double*, const double*, Node*);
void GravitationalAccTree(Node*, const double*, double*, const int);
void Acceleration(Node*, const double*, double*, const int*, const int, const int, const int, const int, const int, MPI_Status);
void Save_data(const double*, const double*, const double*, const int, const int);
void EULER(double*, const double*, const double, const double, const int);
void PEFRL(double*, double*, const double*, double*, const double*, const double*, Node*, const double, const int, const int*, const int, const int, const int, const int, MPI_Status);
void Evolution(body*, body*, double*, double*, Node*, const int, int*, int*, int*, int*, double*, double*, const double, const double, const int, const int, const int, const int, MPI_Status, const int, const double, const int);

#endif // NBODIES_H_