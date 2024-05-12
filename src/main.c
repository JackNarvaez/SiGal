/*-----------------------------------------------------------------------------
The N-Body problem using MPI
-----------------------------------------------------------------------------*/

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "NBodies.h"
#include "structures.h"
#include "Tree.h"
#include "data.h"

int main(int argc, char** argv) {
    int pId;                        // Process rank
    int nP;                         // Number of processes
    int tag{0};                     // Tag message
    int root{0};                    // Root process
    double prmts[6];

    read_parameters("./input", prmts);
    
    int N     = (int) prmts[0];     // Total number of bodies
    double dt = prmts[3];           // Time step
    int steps = (int)prmts[4];      // Evolution steps
    int jump  = (int) prmts[5];     // Data storage interval
    MPI_Status status;
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
    bd.a = (double *) malloc(3*len[pId]*sizeof(double));        // [ax, ay, az]
    bd.m = (double *) malloc(len[pId]*sizeof(double));          // [m]
    
    for (ii = 0; ii < 3*len[pId]; ii++) {
        bd.a[ii] = 0.0;
    }

    // Read local particles information
    char input[32] = "./Data/data";
    sprintf(input + strlen(input), "%d.txt", pId);
    read_data(input, bd.r, bd.v, bd.m);

    double * rootmin = (double *) malloc(3*sizeof(double));        // [x, y, z]
    double * rootmax = (double *) malloc(3*sizeof(double));        // [x, y, z]
    
    for (ii = 0; ii<3; ii++) {
        rootmax[ii] = 10.0;
    }

    for (ii = 0; ii<3; ii++) {
        rootmin[ii] = -1.1*rootmax[ii];
        rootmax[ii] = -rootmin[ii];
    } 

    Node* Tree = BuiltTree(bd.r, bd.m, N, rootmin, rootmax);

    double start = MPI_Wtime();

    // Save positions
    Evolution(bd.r, bd.v, bd.m, bd.a, rootmin, rootmax, Tree, N, len, tag, pId, nP, root, status, steps, dt, jump);
    
    if (pId == root) {
        double end = MPI_Wtime();
        double time_ms = end-start;
        printf("\n%lf\t%d\n", time_ms, nP);
    }

    /*Finalizes MPI*/
    MPI_Finalize();

    free (len);
    free (bd.r);
    free (bd.v);
    free (bd.a);
    free (bd.m);
    free (rootmax);
    free (rootmin);

    return 0;
}