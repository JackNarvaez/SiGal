#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "structures.h"
#include "tree.h"

void DistributeParticles(body* GlobBds, const double* Buffer_r, const double* Buffer_v, const double* Buffer_m, const int* Buffer_i, int* len, int* counts, int* displacements1, int* displacements3, double* GlobMin, double* GlobMax, const double R, const int N, const int nP, const int maxdeep)
{
    double inimin[3], inimax[3];
    int ii;
    for (ii = 0; ii<3; ii++) {
        inimax[ii] = 5.1*R;
        inimin[ii] = -inimax[ii];
    }

    int* ORBTreeid = BuiltGlobalTree(Buffer_r, Buffer_m, N, inimin, inimax, GlobMin, GlobMax, nP, maxdeep);

    int jj, pp;
    for (ii=0; ii<N; ii++) {
        pp = ORBTreeid[nP+ii];
        GlobBds->m[ii] = Buffer_m[pp];
        GlobBds->i[ii] = Buffer_i[pp];
        for (jj=0; jj<3; jj++) {
            GlobBds->r[3*ii+jj] = Buffer_r[3*pp+jj];
            GlobBds->v[3*ii+jj] = Buffer_v[3*pp+jj];
        }
    }

    // Save lens
    for(ii=0; ii<nP; ii++){
        len[ii] = ORBTreeid[ii];
        counts[ii] = 3*len[ii];
    }
    displacements1[0] = 0;
    displacements3[0] = 0;
    for (ii=1; ii<nP; ii++) {
        displacements1[ii] = displacements1[ii-1] + len[ii-1];
        displacements3[ii] = displacements3[ii-1] + counts[ii-1];
    }
    free(ORBTreeid);
}

void GlobalDistri(body* bd, body* GlobBds, double* rootmin, double* rootmax, int* len, int* counts, int* displacements1, int* displacements3, double* GlobMin, double* GlobMax, const double R, const int N, const int nP, const int pId, const int root, const int maxdeep)
{
    free(bd->r);
    free(bd->v);
    free(bd->a);
    free(bd->m);
    free(bd->i);

    if (pId == root) {
        // Buffer
        double* Buffer_r = (double *) malloc(3*N*sizeof(double));   // [x, y, z]
        double* Buffer_v = (double *) malloc(3*N*sizeof(double));   // [vx, vy, vz]
        double* Buffer_m = (double *) malloc(N*sizeof(double));     // [m]
        int* Buffer_i = (int *) malloc(N*sizeof(int));     // [m]
        memcpy(Buffer_r, GlobBds->r, 3*N*sizeof(double));
        memcpy(Buffer_v, GlobBds->v, 3*N*sizeof(double));
        memcpy(Buffer_m, GlobBds->m, N*sizeof(double));
        memcpy(Buffer_i, GlobBds->i, N*sizeof(int));

        DistributeParticles(GlobBds, Buffer_r, Buffer_v, Buffer_m, Buffer_i, len, counts, displacements1, displacements3, GlobMin, GlobMax, R, N, nP,maxdeep);

        free(Buffer_r);
        free(Buffer_v);
        free(Buffer_m);
        free(Buffer_i);
    }

    MPI_Bcast(len, nP, MPI_INT, root, MPI_COMM_WORLD);
    MPI_Bcast(counts, nP, MPI_INT, root, MPI_COMM_WORLD);
    MPI_Bcast(displacements1, nP, MPI_INT, root, MPI_COMM_WORLD);
    MPI_Bcast(displacements3, nP, MPI_INT, root, MPI_COMM_WORLD);

    int max = len[0];
    for (int i = 1; i < nP; i++) { 
        if (max < len[i]) max = len[i]; 
    }

    bd->r = (double *) malloc(3*len[pId]*sizeof(double));        // [x, y, z]
    bd->v = (double *) malloc(3*len[pId]*sizeof(double));        // [vx, vy, vz]
    bd->a = (double *) malloc(3*max*sizeof(double));        // [ax, ay, az]
    bd->m = (double *) malloc(len[pId]*sizeof(double));          // [m]
    bd->i = (int *) malloc(len[pId]*sizeof(int));          // [m]

    int ii;
    for (ii = 0; ii < 3*len[pId]; ii++) {
        bd->a[ii] = 0.0;
    }

    MPI_Scatterv(GlobBds->r, counts, displacements3, MPI_DOUBLE, bd->r, 3*len[pId], MPI_DOUBLE, root, MPI_COMM_WORLD);
    MPI_Scatterv(GlobBds->v, counts, displacements3, MPI_DOUBLE, bd->v, 3*len[pId], MPI_DOUBLE, root, MPI_COMM_WORLD);
    MPI_Scatterv(GlobBds->m, len, displacements1, MPI_DOUBLE, bd->m, len[pId], MPI_DOUBLE, root, MPI_COMM_WORLD);
    MPI_Scatterv(GlobBds->i, len, displacements1, MPI_INT, bd->i, len[pId], MPI_INT, root, MPI_COMM_WORLD);

    MPI_Scatter(GlobMin, 3, MPI_DOUBLE, rootmin, 3, MPI_DOUBLE, root, MPI_COMM_WORLD);
    MPI_Scatter(GlobMax, 3, MPI_DOUBLE, rootmax, 3, MPI_DOUBLE, root, MPI_COMM_WORLD);
}