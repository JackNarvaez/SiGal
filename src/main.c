/*-----------------------------------------------------------------------------
The N-Body problem using MPI
-----------------------------------------------------------------------------*/

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nbodies.h"
#include "structures.h"
#include "tree.h"
#include "data.h"
#include "parallel.h"

#define SETUP_NAME_SIZE 32


int main(int argc, char** argv) {
    int pId;                        // Process rank
    int nP;                         // Number of processes
    int tag{0};                     // Tag message
    int root{0};                    // Root process
    double prmts[6];
    char setup_name[SETUP_NAME_SIZE];  // Buffer para almacenar el nombre del setup


    read_parameters("./input", prmts, setup_name);
    
    int N     = (int) prmts[0];     // Total number of bodies
    double R  = prmts[2];           // Radius
    double dt = prmts[3];           // Time step
    int steps = (int)prmts[4];      // Evolution steps
    int jump  = (int) prmts[5];     // Data storage interval
    MPI_Status status;

    body bd;                        // Bodies

    /*Initializes MPI*/
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nP);
    MPI_Comm_rank(MPI_COMM_WORLD, &pId);

    int maxdeep = (int) log2(nP);
    
    // Length of each proccess
    int len[nP], counts[nP], displacements3[nP], displacements1[nP];

    body GlobBds;
    GlobBds.r = NULL;
    GlobBds.v = NULL;
    GlobBds.m = NULL;

    bd.r = NULL;
    bd.v = NULL;
    bd.a = NULL;
    bd.m = NULL;
    
    double* GlobMin = NULL;
    double* GlobMax = NULL;

    double* rootmin = (double *) malloc(3*sizeof(double));        // [x, y, z]
    double* rootmax = (double *) malloc(3*sizeof(double));        // [x, y, z]

    if (pId == root) {
        GlobMin = (double *) malloc(3*nP*sizeof(double));        // [xmin, ymin, zmin]
        GlobMax = (double *) malloc(3*nP*sizeof(double));        // [xmax, ymax, zmax]

        GlobBds.r = (double *) malloc(3*N*sizeof(double));          // [x, y, z]
        GlobBds.v = (double *) malloc(3*N*sizeof(double));          // [vx, vy, vz]
        GlobBds.m = (double *) malloc(N*sizeof(double));            // [m]

        char input[32] = "./Data/IniData.txt";
        read_data(input, GlobBds.r, GlobBds.v, GlobBds.m);
    }

    GlobalDistri(&bd, &GlobBds, rootmin, rootmax, len, counts, displacements1, displacements3, GlobMin, GlobMax, R, N, nP, pId, root, maxdeep);

    Node* Tree = NULL; //(bd.r, bd.m, len[pId], rootmin, rootmax);

    double start = MPI_Wtime();

    // Save positions
    Evolution(&bd, &GlobBds, rootmin, rootmax, Tree, N, len, counts, displacements1, displacements3, GlobMin, GlobMax, R, maxdeep, tag, pId, nP, root, status, steps, dt, jump);
    
    if (pId == root) {
        free(GlobBds.r);
        free(GlobBds.v);
        free(GlobBds.m);
        free(GlobMin);
        free(GlobMax);
        double end = MPI_Wtime();
        double time_ms = end-start;
        printf("\n%lf\t%d\n", time_ms, nP);
    }

    /*Finalizes MPI*/
    MPI_Finalize();

    free (bd.r);
    free (bd.v);
    free (bd.a);
    free (bd.m);
    free (rootmax);
    free (rootmin);

    return 0;
}