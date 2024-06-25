/*-----------------------------------------------------------------------------
The N-Body problem using MPI
-----------------------------------------------------------------------------*/

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
//#include "plummer.h"
//#include "hernquist.h"
//#include "kepler.h"
#include "Kuzmin.h"
//#include "exponential.h"
#include "../src/structures.h"
#include "../src/data.h"

int main(int argc, char** argv) {
    int pId;                        // Process rank
    int nP;                         // Number of processes
    double prmts[7];

    read_parameters("./input", prmts);
    
    int N     = (int) prmts[0];     // Total number of bodies
    double M  = prmts[1];           // Total Mass of Galaxy
    double R  = prmts[2];           // Parametrized Radius
    int root{0};                    // Root process
    int seed  = 1234;

    body bd;                        // Bodies

    /*Initializes MPI*/
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nP);
    MPI_Comm_rank(MPI_COMM_WORLD, &pId);
    
    // Length of each proccess
    int len[nP], counts[nP], displacements3[nP], displacements1[nP];

    int ii, end, begin;
    for(ii=0; ii<nP; ii++){
        end = double(N)/nP*(ii+1);
        begin = double(N)/nP*ii;
        len[ii] = end-begin;
        counts[ii] = 3*len[ii];
    }
    displacements1[0] = 0;
    displacements3[0] = 0;
    for (ii=1; ii<nP; ii++) {
        displacements1[ii] = displacements1[ii-1] + len[ii-1];
        displacements3[ii] = displacements3[ii-1] + counts[ii-1];
    }

    bd.r = (double *) malloc(3*len[pId]*sizeof(double));        // [x, y, z]
    bd.v = (double *) malloc(3*len[pId]*sizeof(double));        // [vx, vy, vz]
    bd.m = (double *) malloc(len[pId]*sizeof(double));          // [m]
    
    //plummer_dist(bd.r, bd.v, bd.m, len[pId], R, M, seed + pId);
    //hernquist_dist(bd.r, bd.v, bd.m, len[pId], R, M, seed + pId);
    //keplerian_disk(bd.r, bd.v, bd.m, len[pId], R, M, seed + pId);
    //exponential_disk(bd.r, bd.v, bd.m, len[pId], R, M, seed + pId);
    kuzmin_disk(bd.r, bd.v, bd.m, len[pId], R, M, seed + pId);



    body GlobBds;
    GlobBds.r = NULL;
    GlobBds.v = NULL;
    GlobBds.m = NULL;

    if (pId == root) {
        GlobBds.r = (double *) malloc(3*N*sizeof(double));        // [x, y, z]
        GlobBds.v = (double *) malloc(3*N*sizeof(double));        // [vx, vy, vz]
        GlobBds.m = (double *) malloc(N*sizeof(double));          // [m]
    }
    
    MPI_Gatherv(bd.r, 3*len[pId], MPI_DOUBLE, GlobBds.r, counts, displacements3, MPI_DOUBLE, root, MPI_COMM_WORLD);
    MPI_Gatherv(bd.v, 3*len[pId], MPI_DOUBLE, GlobBds.v, counts, displacements3, MPI_DOUBLE, root, MPI_COMM_WORLD);
    MPI_Gatherv(bd.m, len[pId], MPI_DOUBLE, GlobBds.m, len, displacements1, MPI_DOUBLE, root, MPI_COMM_WORLD);

    // Read local particles information
    if (pId == root) {
        char output[32] = "./Data/IniData.txt";
        write_data(output, GlobBds.r, GlobBds.v, GlobBds.m, N);
        free(GlobBds.r);
        free(GlobBds.v);
        free(GlobBds.m);
    }

    /*Finalizes MPI*/
    MPI_Finalize();

    free (bd.r);
    free (bd.v);
    free (bd.m);

    return 0;
}