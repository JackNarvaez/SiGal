#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../src/structures.h"
#include "../src/data.h"
#include "utils.h"  
#include "plummer.h"       
#include "hernquist.h"     
#include "isothermal.h"     
#include "exponential.h"  
#include "kepler.h"       
#include "kuzmin.h"
#include "miyamoto.h"
#include "merger.h"
#include "mwg1.h"

#define SETUP_NAME_SIZE 32


int main(int argc, char** argv) {
    int pId;  // Process rank
    int nP;   // Number of processes
    double prmts[7];
    char setup_name[SETUP_NAME_SIZE];  //buffer for setup name

    // Read input parameters
    read_parameters("./input", prmts, setup_name);

    int N     = (int) prmts[0];     // Total number of bodies
    double M  = prmts[1];           // Total Mass of Galaxy
    double R  = prmts[2];           // Parametrized Radius
    int root{0};                    // Root process
    int seed  = 1234;

    // Inicializar MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nP);
    MPI_Comm_rank(MPI_COMM_WORLD, &pId);

    srand(seed+pId);
    
    // Length of each proccess
    int len[nP], counts[nP], displacements3[nP], displacements1[nP];
    int ii, end, begin;
    for (ii = 0; ii < nP; ii++) {
        end = (double)N / nP * (ii + 1);
        begin = (double)N / nP * ii;
        len[ii] = end - begin;
        counts[ii] = 3 * len[ii];
    }
    displacements1[0] = 0;
    displacements3[0] = 0;
    for (ii = 1; ii < nP; ii++) {
        displacements1[ii] = displacements1[ii - 1] + len[ii - 1];
        displacements3[ii] = displacements3[ii - 1] + counts[ii - 1];
    }

    body bd;

    // Allocate memory for each process
    bd.r = (double *)malloc(3 * len[pId] * sizeof(double));  // [x, y, z]
    bd.v = (double *)malloc(3 * len[pId] * sizeof(double));  // [vx, vy, vz]
    bd.m = (double *)malloc(len[pId] * sizeof(double));      // [m]
    bd.i = (int *)malloc(len[pId] * sizeof(int));      // [m]

    body GlobBds;
    
    if (pId == root) {
        GlobBds.r = (double *) malloc(3*N*sizeof(double));        // [x, y, z]
        GlobBds.v = (double *) malloc(3*N*sizeof(double));        // [vx, vy, vz]
        GlobBds.m = (double *) malloc(N*sizeof(double));          // [m]
        GlobBds.i = (int *) malloc(N*sizeof(int));          // [m]
    } else {
        GlobBds.r = NULL;
        GlobBds.v = NULL;
        GlobBds.m = NULL;
        GlobBds.i = NULL;
    }

    // Initialize the specific setup based on the setup name    
    if (strcmp(setup_name, "PLUMMER") == 0) {
        plummer_dist(bd.r, bd.v, bd.m, bd.i, len[pId], N, R, M, 0);
    } else if (strcmp(setup_name, "HERNQUIST") == 0) {
        hernquist_dist(bd.r, bd.v, bd.m, bd.i, len[pId], N, R, M, 0);
    } else if (strcmp(setup_name, "ISOTHERMAL") == 0) {
        isothermal_dist(bd.r, bd.v, bd.m, bd.i, len[pId], N, R, M, 0);
    } else if (strcmp(setup_name, "EXPONENTIAL") == 0) {
        exponential_disk(bd.r, bd.v, bd.m, bd.i, len[pId], N, R, M, 0);
    } else if (strcmp(setup_name, "KEPLER") == 0) {
        keplerian_disk(bd.r, bd.v, bd.m, bd.i, len[pId], N, R, M, 0);
        if (pId==root) {
            bd.m[0] = M;
            bd.r[0] = 0.0; bd.r[1] = 0.0; bd.r[2] = 0.0;
            bd.v[0] = 0.0; bd.v[1] = 0.0; bd.v[2] = 0.0;
        }
    } else if (strcmp(setup_name, "KUZMIN") == 0) {
        kuzmin_disk(bd.r, bd.v, bd.m, bd.i, len[pId], N, R, M, 0);
    } else if (strcmp(setup_name, "MIYAMOTO") == 0) {
        miyamoto_disk(bd.r, bd.v, bd.m, bd.i, len[pId], N, R, M, 0);
    }else if (strcmp(setup_name, "MERGER") == 0) {
        galaxy_collision(bd.r, bd.v, bd.m, bd.i, len[pId], N, R, M, 0);
    }else if (strcmp(setup_name, "MWG1") == 0) {
        mwg1(bd.r, bd.v, bd.m, bd.i, len[pId], N, R, M, 0);
    }else {
        fprintf(stderr, "Unknown setup: %s\n", setup_name);
        MPI_Finalize();
        return 1;
    }
    // Collect processed data at the root process
    MPI_Gatherv(bd.r, 3*len[pId], MPI_DOUBLE, GlobBds.r, counts, displacements3, MPI_DOUBLE, root, MPI_COMM_WORLD);
    MPI_Gatherv(bd.v, 3*len[pId], MPI_DOUBLE, GlobBds.v, counts, displacements3, MPI_DOUBLE, root, MPI_COMM_WORLD);
    MPI_Gatherv(bd.m, len[pId], MPI_DOUBLE, GlobBds.m, len, displacements1, MPI_DOUBLE, root, MPI_COMM_WORLD);
    MPI_Gatherv(bd.i, len[pId], MPI_INT, GlobBds.i, len, displacements1, MPI_INT, root, MPI_COMM_WORLD);

    // Save data at the root process
    if (pId == root) {
        char output[32] = "./Data/IniData.txt";
        write_data(output, GlobBds.r, GlobBds.v, GlobBds.m, GlobBds.i, N);
        free(GlobBds.r);
        free(GlobBds.v);
        free(GlobBds.m);
        free(GlobBds.i);
    }

    // Finalize MPI
    MPI_Finalize();

    // Free memory on each process
    free(bd.r);
    free(bd.v);
    free(bd.m);
    free(bd.i);

    return 0;
}
