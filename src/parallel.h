#ifndef PARALLEL_H_
#define PARALLEL_H_
#include "structures.h"

void DistributeParticles(body* GlobBds, const double* Buffer_r, const double* Buffer_v, const double* Buffer_m, const double* Buffer_i, int* len, int* counts, int* displacements1, int* displacements3, double* GlobMin, double* GlobMax, const double R, const int N, const int nP, const int maxdeep);
void GlobalDistri(body* bd, body* GlobBds, double* rootmin, double* rootmax, int* len, int* counts, int* displacements1, int* displacements3, double* GlobMin, double* GlobMax, const double R, const int N, const int nP, const int pId, const int root, const int maxdeep);

#endif // PARALLEL_H_