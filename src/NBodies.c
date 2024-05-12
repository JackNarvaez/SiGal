#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "structures.h"
#include "Tree.h"

// Global Constants
const double G      = 1.0;
const double EPSG   = 0.025;
const double EPSG2  = EPSG*EPSG;
const double THETA  = 0.5;
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

    // Find the maximum side length of the node
    double side = node->max[0] - node->min[0];
    if (node->max[1] - node->min[1] > side) side = node->max[1] - node->min[1];
    if (node->max[2] - node->min[2] > side) side = node->max[2] - node->min[2];

    bool cond = (side*side / d2) < THETA2;

    return cond;
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

        while (node != NULL) {
            if (node == root && lastVisited != NULL) {
                break; 
            }
            //Leaf Node
            if (node->type) {
                if (*node->Mass > 0) {     // Not empty leaf
                    Forcei(&Acc[3*ii], &Pos[3*ii], node);
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
                    if (FarDistance(&Pos[3*ii], node)) {
                        Forcei(&Acc[3*ii], &Pos[3*ii], node); 
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
    int max = N/nP+1;

    //  Temporal arrays for saving data of shared bodies along the ring
    double *BufferPos = (double *) malloc(3*(len[pId]+1)*sizeof(double));

    int ii, jj;
    for (ii=0; ii < len[pId]; ii++) {
        for (jj=0; jj < 3; jj++) {
            BufferPos[3*ii + jj] = Pos[3*ii + jj];
        }
    }
    for (jj=0; jj < 3; jj++) {
        BufferPos[3*len[pId] + jj] = 0.0;
    }
  
    /*The Ring Loop*/
    int dst = (pId+1)%nP;
    int scr = (pId-1+nP)%nP;
  
    for (jj=0; jj<nP; jj++){
        MPI_Sendrecv_replace(BufferPos, 3*max, MPI_DOUBLE, dst, tag, scr, tag, MPI_COMM_WORLD, &status);
        MPI_Sendrecv_replace(Acc, 3*max, MPI_DOUBLE, dst, tag, scr, tag, MPI_COMM_WORLD, &status);
        GravitationalAccTree(Tree, BufferPos, Acc, len[(scr-jj+nP)%nP]);
    }

    free(BufferPos);
}

void Save_vec(FILE *File, const double *Pos, const double *Mass, const int N)
{
    /*---------------------------------------------------------------------------
    Saves info from Pos and Mass in File.
    -----------------------------------------------------------------------------
    Arguments:
      File  :   File pointer where data is saved.
      Pos   :   Positions.
      Mass  :   Masses.
      N     :   Size of Mass.
    ---------------------------------------------------------------------------*/
    int ii;
    for (ii = 0; ii < N; ii++) {
        fprintf(File, "%.15lf\t%.15lf\t%.15lf\t%.15lf\n", Pos[3 * ii], Pos[3 * ii + 1], Pos[3 * ii + 2], Mass[ii]);
    }
}

void Save_data(const double * Pos, const double * Mass, const int *len, const int N, const int tag, const int pId, const int nP, const int root, MPI_Status status)
{
    /*---------------------------------------------------------------------------
    Saves positions and masses of all bodies in Evolution.txt.
    -----------------------------------------------------------------------------
    Arguments:
      File  :   File where data is saved.
      Pos   :   Position of bodies (1D vector).
      Mass  :   Mass of bodies (1D vector).
      len   :   Array with the number of bodies per node.
      N     :   Number of bodies.
      tag   :   Message tag.
      pId   :   Process identity.
      nP    :   Number of processes.
      root  :   Root process.
      status:   Status object.
    ---------------------------------------------------------------------------*/

    // Collect results in <root> process.
    if(pId==root){
        FILE *File;
        File = fopen("./Data/Evolution.txt", "a");
        Save_vec(File, Pos, Mass, len[root]);
        double * TempP = (double *) malloc(3*(N/nP+1)*sizeof(double));
        double * TempM = (double *) malloc((N/nP+1)*sizeof(double));
        int ii;
        for (ii = 0; ii < N/nP+1; ii++) {
            TempP[3*ii] = 0.0;
            TempP[3*ii+1] = 0.0;
            TempP[3*ii+2] = 0.0;
            TempM[ii] = 0.0;
        }
        for (ii = 0; ii < nP; ii++){
            if (ii != pId){
            MPI_Recv(TempP, 3*len[ii], MPI_DOUBLE, ii, tag, MPI_COMM_WORLD, &status);
            MPI_Recv(TempM, len[ii], MPI_DOUBLE, ii, tag, MPI_COMM_WORLD, &status);
            Save_vec(File, TempP, TempM,  len[ii]);
            }
        }
        fclose(File);
        free(TempP);
        free(TempM);
    } else {
        MPI_Send(Pos, 3*len[pId], MPI_DOUBLE, root, tag, MPI_COMM_WORLD);
        MPI_Send(Mass, len[pId], MPI_DOUBLE, root, tag, MPI_COMM_WORLD);
    }
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
    -----------------------------------------------------------------------------
    Arguments:
      Pos   :   Position of bodies (1D vector).
      Vel   :   Velocity of bodies (1D vector).
      Mass  :   Mass of bodies (1D vector).
      Acc   :   Accelerations (1D vector).
      dt    :   Time step.
      N     :   Total number of bodies.
      len   :   Array with the number of bodies per node.
      tag   :   Message tag.
      pId   :   Process identity.
      nP    :   Number of processes.
      root  :   Root process.
      status:   Status object.
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

void Evolution(double* Pos, double* Vel, const double* Mass, double* Acc, const double* rootmin, const double* rootmax, Node *Tree, const int N, const int * len, const int tag, const int pId, const int nP, const int root, MPI_Status status, const int steps, const double dt, const int jump)
{
    /*---------------------------------------------------------------------------
    Evolution of the system of bodies under gravitational interactions.
    -----------------------------------------------------------------------------
    Arguments:
      File  :   File where data is saved.
      Pos   :   Position of bodies (1D vector).
      Vel   :   Velocity of bodies (1D vector).
      Mass  :   Mass of bodies (1D vector).
      Acc   :   Accelerations (1D vector).
      len   :   Array with the number of bodies per node.
      N     :   Total number of bodies.
      tag   :   Message tag.
      pId   :   Process identity.
      nP    :   Number of processes.
      root  :   Root process.
      status:   Status object.
      steps :   Evolution steps.
      dt    :   Time step.
      jump  :   Data storage interval.
      Accel :   Function to calculate acceleration.
      evol  :   Integrator.
    ---------------------------------------------------------------------------*/
  
    if (pId==root){
        remove("Evolution.txt");
    }
    Save_data(Pos, Mass, len, N, tag, pId, nP, root, status);
    int ii, jj;
    for (ii = 0; ii < steps; ii++){
        PEFRL(Pos, Vel, Mass, Acc, rootmin, rootmax, Tree, dt, N, len, tag, pId, nP, root, status);
        if ( ii%jump == 0) {
            Save_data(Pos, Mass, len, N, tag, pId, nP, root, status);
            if (pId==root) printf("%d\n", ii/jump);
        }
        for (jj = 0; jj < 3*len[pId]; jj++) Acc[jj] = 0.0;
    }
}