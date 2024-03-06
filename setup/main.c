/*-----------------------------------------------------------------------------
The N-Body problem using MPI
-----------------------------------------------------------------------------*/

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "plummer.h"
#include "../src/structures.h"
#include "../src/data.h"

int main(int argc, char** argv) {
    int pId;                        // Process rank
    int nP;                         // Number of processes
    double prmts[7];

    read_parameters("./input", prmts);
    
    int N     = (int) prmts[0];      // Total number of bodies

    body bd;                        // Bodies

    /*Initializes MPI*/
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nP);
    MPI_Comm_rank(MPI_COMM_WORLD, &pId);
    
    // Length of each proccess
    int * len = (int *) malloc(nP * sizeof(int));

    int ii, end, begin;
    for(ii=0; ii<nP; ii++){
        end = double(N)/nP*(ii+1);
        begin = double(N)/nP*ii;
        len[ii] = end-begin;
    }

    bd.r = (double *) malloc(3*len[pId]*sizeof(double));        // [x, y, z]
    bd.v = (double *) malloc(3*len[pId]*sizeof(double));        // [vx, vy, vz]
    bd.m = (double *) malloc(len[pId]*sizeof(double));          // [m]
    
    plummer_dist(bd.r, bd.v, bd.m, len[pId]);

    // Read local particles information
    char output[32] = "./Data/data";
    sprintf(output + strlen(output), "%d.txt", pId);
    write_data(output, bd.r, bd.v, bd.m, len[pId]);

    /*Finalizes MPI*/
    MPI_Finalize();

    free (len);
    free (bd.r);
    free (bd.v);
    free (bd.m);

    return 0;
}