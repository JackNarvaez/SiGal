/*-----------------------------------------------------------------------------
The N-Body problem
-----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "structures.h"
#include "Tree.h"

const double THETA  = 0.5;
const double G      = 1.0;


void treeforce(Node* root, double* bdr, double* bda,const double* bdm, const int N);

bool separated(const double * bdr, const Node* node);

void force(double * bda, const double * bdr, const double * bdm, Node* node);

void Euler(double *bdr, double *bdv, double *bda, const int N, const double dt);

void PEFRL(double *bdr, double *bdv, double *bda, double *bdm, const int N, double dt, double *rootmin, double *rootmax, double eps, Node *Tree);

void read_data(const char *File_address, double *Pos, double *Vel, double *Mass);

void writeTXT(double *bdr,double *bdm,  int N, const char *filename);

void printtree(Node* node);


int main(int argc, char** argv) {
    
    const int N = 100;     // Total number of bodies
    body bd;                        // Bodies
    double eps = 1.e-05;    // stop parameter for node's side

    bd.r = (double *) malloc(3*N*sizeof(double));        // [x, y, z]
    bd.rtemp = (double *) malloc(3*N*sizeof(double));    // [x, y, z]
    bd.v = (double *) malloc(3*N*sizeof(double));        // [vx, vy, vz]
    bd.vtemp = (double *) malloc(3*N*sizeof(double));    // [vx, vy, vz]
    bd.a = (double *) malloc(3*N*sizeof(double));        // [ax, ay, az]
    bd.m = (double *) malloc(N*sizeof(double));          // [m]

    int ii, jj;
    
    for (ii = 0; ii < 3*N; ii++) {
        bd.rtemp[ii] = 0.0;
        bd.vtemp[ii] = 0.0;
        bd.a[ii] = 0.0;
    }

    // Read local particles information
    char input[32] = "data0.txt";
    read_data(input, bd.r, bd.v, bd.m);

    double * rootmin = (double *) malloc(3*sizeof(double));        // [x, y, z]
    double * rootmax = (double *) malloc(3*sizeof(double));        // [x, y, z]
    
    for (jj = 0; jj<3; jj++) {
        rootmax[jj] = 0.0;
    }
    
    // Ajust limits of RootNode
    for (ii = 0; ii<N; ii++) {
        for (jj = 0; jj<3; jj++) {
            if (fabs(bd.r[3*ii+jj]) > rootmax[jj]) {
                rootmax[jj] = fabs(bd.r[3*ii+jj]);
            }
        }
    }
    for (jj = 0; jj<3; jj++) {
        rootmin[jj] = -1.1*rootmax[jj];
        rootmax[jj] = -rootmin[jj];
    } 

    Node* Tree = BuiltTree(bd.r, bd.m, N, rootmin, rootmax, eps);

    // Integration Paremeters
    double dt = 0.001; 
    int Nsteps= 10000;
    int jump = 100;

    // File Evolution.txt
    char filename[256] = "Evolution.txt";
    FILE* init_file = fopen(filename, "w");
    if (init_file != NULL) {
        fprintf(init_file, "X,Y,Z,m\n");
        fclose(init_file);
    }


    // Calculate force
    for (int ii = 0; ii < Nsteps; ii++) { 
        memset(bd.a, 0, 3 * N * sizeof(double));
        treeforce(Tree, bd.r, bd.a, bd.m, N);
        Euler(bd.r, bd.v, bd.a, N, dt);

        //PEFRL(bd.r, bd.v, bd.a,  bd.m, N, dt, rootmin,  rootmax, eps, Tree);

        if (ii % jump == 0) { 
            writeTXT(bd.r,bd.m,  N, "Evolution.txt");
        }
        printf("\n\n step = %d\n\n", ii);
        
        // Rebuild tree with new positions
        freeNode(Tree);
        Tree = BuiltTree(bd.r, bd.m, N, rootmin, rootmax, eps);
    }

    freeNode(Tree);
    

    // Save positions
    free (bd.r);
    free (bd.rtemp);
    free (bd.v);
    free (bd.vtemp);
    free (bd.a);
    free (bd.m);
    free (rootmax);
    free (rootmin);

    return 0;
}


void treeforce(Node* root, double* bdr, double* bda,const double* bdm, const int N) {

    for (int ii = 0; ii < N; ++ii) {
        Node* node = root;
        Node* lastVisited = NULL;  // Storage for the last visited node

        while (node != NULL) {
            if (node == root && lastVisited != NULL) {
                break; 
            }
            //Leaf Node
            if (node->type) {
                if (*node->Mass > 0) {     // Not empty leaf
                    force(&bda[3*ii], &bdr[3*ii], &bdm[ii], node);
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
                    if (separated(&bdr[3*ii], node)) {
                        force(&bda[3*ii], &bdr[3*ii], &bdm[ii], node); 
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


bool separated(const double * bdr, const Node* node) {

    // Calculate the relative position from the particle to the node's center of mass
    double xrel = node->CoM[0] - bdr[0];
    double yrel = node->CoM[1] - bdr[1];
    double zrel = node->CoM[2] - bdr[2];

    double d = sqrt(xrel * xrel + yrel * yrel + zrel * zrel);

    // Find the maximum side length of the node
    double side = fmax(node->max[0] - node->min[0],
                      fmax(node->max[1] - node->min[1],
                           node->max[2] - node->min[2]));

    bool Separated = (side / d < THETA);

    return Separated;
}


void force(double * bda, const double * bdr, const double * bdm, Node* node) {

    double xrel = node->CoM[0] - bdr[0];
    double yrel = node->CoM[1] - bdr[1];
    double zrel = node->CoM[2] - bdr[2];
    double epsilon = 0.25;  // Softening parameter
    double d = sqrt(xrel * xrel + yrel * yrel + zrel * zrel + epsilon * epsilon); 
    double F = G * *node->Mass / (d * d * d);    
    bda[0] += F * xrel;
    bda[1] += F * yrel;
    bda[2] += F * zrel;

}


void Euler(double *bdr, double *bdv, double *bda, int N, double dt) {
    for (int ii = 0; ii < N; ii++) {
        bdv[3*ii] += bda[3*ii] * dt; 
        bdr[3*ii] += bdv[3*ii] * dt;  

        bdv[3*ii + 1] += bda[3*ii + 1] * dt; 
        bdr[3*ii + 1] += bdv[3*ii + 1] * dt;  

        bdv[3*ii + 2] += bda[3*ii + 2] * dt;  
        bdr[3*ii + 2] += bdv[3*ii + 2] * dt; 

    }
}


void PEFRL(double *bdr, double *bdv, double *bda, double *bdm, const int N, double dt, double *rootmin, double *rootmax, double eps, Node *Tree) {
    /*---------------------------------------------------------------------------
    Position Extended Forest-Ruth Like method to calculate position and velocity
    at next time step.
    -----------------------------------------------------------------------------
    Arguments:
      bdr   :   Position of bodies (1D vector).
      bdv   :   Velocity of bodies (1D vector).
      bdm  :   Mass of bodies (1D vector).
      bda   :   Accelerations (1D vector).
      dt    :   Time step.
      N     :   Total number of bodies.
      treeforce :   Function to calculate acceleration.

    ---------------------------------------------------------------------------*/
    // PEFRL Parameters
    double xi = 0.1786178958448091;
    double gamma = -0.2123418310626054;
    double chi = -0.6626458266981849;


    // n arrays

    double *X = (double *)malloc(3 * N * sizeof(double));
    double *V = (double *)malloc(3 * N * sizeof(double));

    // Copy original positions and velocities to temporary arrays
    memcpy(X, bdr, 3 * N * sizeof(double));
    memcpy(V, bdv, 3 * N * sizeof(double));

    
    Tree = BuiltTree(X, bdm, N, rootmin, rootmax, eps);

    // Main loop
    int ii;

    // Update position
    for(ii = 0; ii < N; ii++){
        X[3*ii]   += xi * dt * V[3*ii];
        X[3*ii+1] += xi * dt * V[3*ii+1];
        X[3*ii+2] += xi * dt * V[3*ii+2];
    }
    
    // Calculate new forces again
    freeNode(Tree);
    Tree = BuiltTree(X, bdm, N, rootmin, rootmax, eps);
    treeforce(Tree, X, bda, bdm, N);  


    // update Velocities
    for (ii = 0; ii < N; ii++){
        V[3*ii]   +=  0.5 * (1.-2*gamma) * dt * bda[3*ii];
        V[3*ii+1] +=  0.5 * (1.-2*gamma) * dt * bda[3*ii+1];
        V[3*ii+2] +=  0.5 * (1.-2*gamma) * dt * bda[3*ii+2];
    }
    
    // Update position 
    for (ii = 0; ii < N; ii++){
        X[3*ii]   += chi * dt * V[3*ii];
        X[3*ii+1] += chi * dt * V[3*ii+1];
        X[3*ii+2] += chi * dt * V[3*ii+2];
    }

    memset(bda, 0, 3 * N * sizeof(double));

  
    // Calculate new forces again
    freeNode(Tree);
    Tree = BuiltTree(X, bdm, N, rootmin, rootmax, eps);
    treeforce(Tree, X, bda, bdm, N);


    // update Velocities
    for (ii = 0; ii < N; ii++){
        V[3*ii]   += gamma * dt*bda[3*ii];
        V[3*ii+1] += gamma * dt*bda[3*ii+1];
        V[3*ii+2] += gamma * dt*bda[3*ii+2];
    }

    // Update position
    for (ii = 0; ii < N; ii++){
        X[3*ii]   += (1.-2*(chi+xi)) * dt * V[3*ii];
        X[3*ii+1] += (1.-2*(chi+xi)) * dt * V[3*ii+1];
        X[3*ii+2] += (1.-2*(chi+xi)) * dt * V[3*ii+2];
    }
    memset(bda, 0, 3 * N * sizeof(double));

    // Calculate new forces again
    freeNode(Tree);
    Tree = BuiltTree(X, bdm, N, rootmin, rootmax, eps);
    treeforce(Tree, X, bda, bdm, N);

    // update Velocities
    for (ii = 0; ii < N; ii++){
        V[3*ii]   += gamma * dt * bda[3*ii];
        V[3*ii+1] += gamma * dt * bda[3*ii+1];
        V[3*ii+2] += gamma * dt * bda[3*ii+2];
    }
    // Update position
    for (ii = 0; ii < N; ii++){
        X[3*ii]   += chi * dt * V[3*ii];
        X[3*ii+1] += chi * dt * V[3*ii+1];
        X[3*ii+2] += chi * dt * V[3*ii+2];
    }
    memset(bda, 0, 3 * N * sizeof(double));
  
    // Calculate new forces again
    freeNode(Tree);
    Tree = BuiltTree(X, bdm, N, rootmin, rootmax, eps);
    treeforce(Tree, X, bda, bdm, N);

    // update real velocities
    for (ii = 0; ii < N; ii++){
        bdv[3*ii]   = V[3*ii]   + 0.5 * (1.-2*gamma) * dt * bda[3*ii];
        bdv[3*ii+1] = V[3*ii+1] + 0.5 * (1.-2*gamma) * dt * bda[3*ii+1];
        bdv[3*ii+2] = V[3*ii+2] + 0.5 * (1.-2*gamma) * dt * bda[3*ii+2];
    }
    // update real positions
    for (ii = 0; ii < N; ii++){
        bdr[3*ii]   = X[3*ii]   + xi * dt * bdv[3*ii];
        bdr[3*ii+1] = X[3*ii+1] + xi * dt * bdv[3*ii+1];
        bdr[3*ii+2] = X[3*ii+2] + xi * dt * bdv[3*ii+2];
    }

    //freeNode(Tree);

    free(X);
    free(V);
}


void writeTXT(double *bdr,double *bdm,  int N, const char *filename) {

    FILE* file = fopen(filename, "a"); 
    if (file == NULL) {
        printf("Error to open file...\n");
        return;
    }

    for (int ii = 0; ii < N; ii++) {
        fprintf(file, "%.15lf\t%.15lf\t%.15lf\t%.15lf\n", bdr[3 * ii], bdr[3 * ii + 1], bdr[3 * ii + 2], bdm[ii]); 
    }

    fclose(file); 
}


void read_data(const char *File_address, double *Pos, double *Vel, double *Mass) {
    /*---------------------------------------------------------------------------
    Reads data from bodies: position, velocity and mass.
    -----------------------------------------------------------------------------
    Arguments:
      File_address:   File address from where the data is read.
      Pos         :   Position (1D vector) [xi, yi, zi].
      Vel         :   Velocity (1D vector) [vxi, vyi, vzi].
      Mass        :   Mass (1D vector) [mi].
    ---------------------------------------------------------------------------*/

    FILE *File;
    File = fopen(File_address, "r");
    char line[256];
    int row = 0;
    int ii;
    while (fgets(line, sizeof(line), File)) {
        if (line[0] == '\n' || line[0] == '#')
            continue;
        else {
            char *token;
            token = strtok(line, "\t");
            for (ii = 0; ii < 3; ii++) {
                Pos[ii + 3 * row] = atof(token);
                token = strtok(NULL, "\t");
            }
            for (ii = 0; ii < 3; ii++) {
                Vel[ii + 3 * row] = atof(token);
                token = strtok(NULL, "\t");
            }
            Mass[row] = atof(token);
            row += 1;
        }
    }
    fclose(File);
}


void printtree(Node* node) {
    if (node->type) {
        printf("%d \t %f \t %f \t %f \t %f \t %f \t %f \n", *node->deep, node->min[0], node->min[1], node->min[2],
                                     node->max[0],  node->max[1], node->max[2]);
    } else {
        printtree(node->child);
        printtree(node->sibling);

    }
}

