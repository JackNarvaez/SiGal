#include <mpi.h>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "structures.h"
#include "parallel.h"
#include "tree.h"

// Global Constants
const double G      = 1.0;
const double EPSG   = 0.0001;
const double EPSG2  = EPSG*EPSG;
const double THETA  = 0.75;
const double THETA2 = THETA*THETA;

// PEFRL Parameters
const double XI     = 0.1786178958448091;
const double LAMBDA = -0.2123418310626054;
const double CHI    = -0.06626458266981849;

const double UM2LAMBDAU2 = 0.5*(1.0-2*LAMBDA);
const double UM2CHIXI = 1.0-2*(CHI+XI);

// Distance criterion for computing of gravity force
bool FarDistance(const double* Pos, const Node* node)
{
    double xrel = node->CoM[0] - Pos[0];
    double yrel = node->CoM[1] - Pos[1];
    double zrel = node->CoM[2] - Pos[2];

    double d2   = xrel * xrel + yrel * yrel + zrel * zrel;

    // Find the longest side of the node
    double dx = node->max[0] - node->min[0];
    double dy = node->max[1] - node->min[1];
    double dz = node->max[2] - node->min[2];
    double side = std::max({dx, dy, dz});

    return side*side < THETA2*d2;
}

// Compute gravitation force for Body(Pos) due to Node
void Forcei(double* Acc, const double* Pos, Node* node)
{
    double dx, dy, dz, dq2, inv_rtd2, cb_d2, F;
    dx  = node->CoM[0] - Pos[0];
    dy  = node->CoM[1] - Pos[1];
    dz  = node->CoM[2] - Pos[2];
    dq2 = dx*dx + dy*dy + dz*dz + EPSG2;
    inv_rtd2 = 1./sqrt(dq2);
    cb_d2 = inv_rtd2*inv_rtd2*inv_rtd2;
    F = G* *node->Mass*cb_d2;
    Acc[0] += F*dx;
    Acc[1] += F*dy;
    Acc[2] += F*dz;
}

// Compute gravitational force for all local bodies
void GravitationalAccTree(Node* root, const double* Pos, double* Acc, const int Nloc)
{
    int ii;
    for (ii = 0; ii < Nloc; ++ii) {
        Node* node = root;
        Node* lastVisited = NULL;  // Last visited node

        // Store positions and acceleration pointers
        double * acc = &Acc[3 * ii];
        const double * pos = &Pos[3 * ii];

        while (node != NULL) {
            if (node == root && lastVisited != NULL) {
                break; 
            }
            //Leaf Node
            if (node->type) {
                if (*node->Mass > 0) {     // Not empty leaf
                    Forcei(acc, pos, node);
                }
            lastVisited = node;  
            node = node->next;         // Move to sibling or parent node
            } 
            // Internal node
            else {
                if (lastVisited == node->sibling) {  
                    lastVisited = node;    // Last visited is the sibling (down to up)
                    node = node->next;
                } else {
                    if (FarDistance(pos, node)) {
                        Forcei(acc, pos, node); 
                        lastVisited = node;  
                        node = node->next;
                    }else {
                        lastVisited = node;  
                        node = node->child;  // Move to the first child
                    }
                }
            }
        }
    }
}

// Ring Method for global calculation of the accelerations
void Acceleration(Node* Tree, const double* Pos, double* Acc, const int* len, const int N, const int tag, const int pId, const int nP, const int root, MPI_Status status)
{
    int max = len[0];
    for (int i = 1; i < nP; i++) { 
        if (max < len[i]) max = len[i]; 
    }


    //  Temporal arrays for saving data of shared bodies along the ring
    double *BufferPos = (double *) malloc(3*max*sizeof(double));

    int ii, jj;
    for (ii=0; ii < len[pId]; ii++) {
        for (jj=0; jj < 3; jj++) {
            BufferPos[3*ii + jj] = Pos[3*ii + jj];
        }
    }
    for (ii=len[pId]; ii < max; ii++) {
        for (jj=0; jj < 3; jj++) {
            BufferPos[3*ii + jj] = 0.0;
        }
    }
  
    /*The Ring Loop*/
    int dst = (pId+1)%nP;
    int src = (pId-1+nP)%nP;

    for (jj=0; jj<nP; jj++){
        MPI_Sendrecv_replace(BufferPos, 3*max, MPI_DOUBLE, dst, tag, src, tag, MPI_COMM_WORLD, &status);
        MPI_Sendrecv_replace(Acc, 3*max, MPI_DOUBLE, dst, tag, src, tag, MPI_COMM_WORLD, &status);
        GravitationalAccTree(Tree, BufferPos, Acc, len[(src-jj+nP)%nP]);
    }
    
    free(BufferPos);
}

void Save_data(const double *Pos, const double *Vec, const double *Mass, const int *i, const int N, const int it)
{
    /*---------------------------------------------------------------------------
    Saves positions, velocities, and masses of all bodies in Evolution.
    ---------------------------------------------------------------------------*/

    FILE *File;
    char filename[32] = "./Data/Ev_";
    sprintf(filename + strlen(filename), "%d", it);  // Add .bin extension for binary file
    File = fopen(filename, "wb");

    // Check if file opened successfully
    if (File == NULL) {
        fprintf(stderr, "Error opening file %s\n", filename);
        return;
    }

    fwrite(Pos, sizeof(double), 3 * N, File);
    fwrite(Vec, sizeof(double), 3 * N, File);
    fwrite(Mass, sizeof(double), N, File);
    fwrite(i, sizeof(int), N, File);

    // Close the file
    fclose(File);
}


void EULER(double* R, const double* DR, const double coeff, const double dt, const int Nloc)
{
    int ii, jj;
    for (ii = 0; ii < Nloc; ii++) {
        for (jj = 0; jj < 3; jj++) {
            R[3*ii+jj] +=  coeff*dt*DR[3*ii+jj];
        }
    }
}

void PEFRL(double* Pos, double* Vel, const double* Mass, double*Acc, const double* rootmin, const double* rootmax, Node* Tree, const double dt, const int N, const int * len, const int tag, const int pId, const int nP, const int root, MPI_Status status)
{
    /*---------------------------------------------------------------------------
    Position Extended Forest-Ruth Like method to calculate position and velocity
    at next time step.
    ---------------------------------------------------------------------------*/
    int Nloc = len[pId];

    EULER(Pos, Vel, XI, dt, Nloc);
    memset(Acc, 0, 3 * Nloc * sizeof(double));

    // Move 1
    Tree = BuiltTree(Pos, Mass, Nloc, rootmin, rootmax);
    Acceleration(Tree, Pos, Acc, len, N, tag, pId, nP, root, status);
    freeNode(Tree);
    EULER(Vel, Acc, UM2LAMBDAU2, dt, Nloc);
    EULER(Pos, Vel, CHI, dt, Nloc);
    memset(Acc, 0, 3 * Nloc * sizeof(double));
  
    // Move 2
    Tree = BuiltTree(Pos, Mass, Nloc, rootmin, rootmax);
    Acceleration(Tree, Pos, Acc, len, N, tag, pId, nP, root, status);
    freeNode(Tree); 
    EULER(Vel, Acc, LAMBDA, dt, Nloc);
    EULER(Pos, Vel, UM2CHIXI, dt, Nloc);
    memset(Acc, 0, 3 * Nloc * sizeof(double));

    // Move 3
    Tree = BuiltTree(Pos, Mass, Nloc, rootmin, rootmax);
    Acceleration(Tree, Pos, Acc, len, N, tag, pId, nP, root, status);
    freeNode(Tree); 
    EULER(Vel, Acc, LAMBDA, dt, Nloc);
    EULER(Pos, Vel, CHI, dt, Nloc);
    memset(Acc, 0, 3 * Nloc * sizeof(double));

    // Move 4
    Tree = BuiltTree(Pos, Mass, Nloc, rootmin, rootmax);
    Acceleration(Tree, Pos, Acc, len, N, tag, pId, nP, root, status);
    freeNode(Tree);
    EULER(Vel, Acc, UM2LAMBDAU2, dt, Nloc);
    EULER(Pos, Vel, XI, dt, Nloc);
}

void Evolution(body* bd, body* GlobBds, double* rootmin, double* rootmax, Node *Tree, const int N, int * len, int* counts, int* displacements1, int* displacements3, double* GlobMin, double* GlobMax, const double R, const double maxdeep, const int tag, const int pId, const int nP, const int root, MPI_Status status, const int steps, const double dt, const int jump)
{
    /*---------------------------------------------------------------------------
    Evolution of the system of bodies under gravitational interactions.
    ---------------------------------------------------------------------------*/
    if (pId==root) Save_data(GlobBds->r, GlobBds->v, GlobBds->m, GlobBds->i, N, 0);
    int ii;
    for (ii = 1; ii < steps+1; ii++){    
        PEFRL(bd->r, bd->v, bd->m, bd->a, rootmin, rootmax, Tree, dt, N, len, tag, pId, nP, root, status);
        MPI_Gatherv(bd->r, 3*len[pId], MPI_DOUBLE, GlobBds->r, counts, displacements3, MPI_DOUBLE, root, MPI_COMM_WORLD);
        MPI_Gatherv(bd->v, 3*len[pId], MPI_DOUBLE, GlobBds->v, counts, displacements3, MPI_DOUBLE, root, MPI_COMM_WORLD);
        MPI_Gatherv(bd->m, len[pId], MPI_DOUBLE, GlobBds->m, len, displacements1, MPI_DOUBLE, root, MPI_COMM_WORLD);
        MPI_Gatherv(bd->i, len[pId], MPI_INT, GlobBds->i, len, displacements1, MPI_INT, root, MPI_COMM_WORLD);
        GlobalDistri(bd, GlobBds, rootmin, rootmax, len, counts, displacements1, displacements3, GlobMin, GlobMax, R, N, nP, pId, root, maxdeep);
        if (pId==root) {
            if ( ii%jump == 0) {
                printf("It: %d\t t: %.4f\t dt: %.4f\n", ii, ii*dt, dt);
                Save_data(GlobBds->r, GlobBds->v, GlobBds->m, GlobBds->i, N, ii/jump);
            }
        }
    }
}