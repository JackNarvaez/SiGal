#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include "structures.h"

Node* Init_Node(double *min, double *max){

    //Memory allocation
    Node* node = (Node *)malloc(sizeof(Node));
    if (node == NULL) {
        fprintf(stderr, "Error Allocation Memory Init_Node\n");
        exit(EXIT_FAILURE);
    }

    node->CoM       = (double *)malloc(3 * sizeof(double));  // x, y, z
    node->totalMass = (double *)malloc(sizeof(double));      // mass
    node->bbox[0]   = (double *)malloc(3 * sizeof(double));       // min[x,y,z]
    node->bbox[1]   = (double *)malloc(3 * sizeof(double));       // max[x,y,z]

    if (node->CoM == NULL || node->totalMass == NULL || node->bbox[0] == NULL || node->bbox[1] == NULL) {
        fprintf(stderr, "Error Allocation Memory in Init_Node\n");
        exit(EXIT_FAILURE);
    }

    // Initialization
    for(int ii = 0; ii < 3; ++ii) {
        node->bbox[0][ii] = min[ii]; 
        node->bbox[1][ii] = max[ii];
        node->CoM[ii] = 0.0;
    }

    *node->totalMass = 0.0; 
    
    for(int ii = 0; ii < 8; ++ii) {
        node->children[ii] = NULL;
    }
    node->bds = NULL;

    return node;
}

Node* createChildNodeForOctant(Node* node, int octantIndex) {

    // Center of Node
    double xc = 0.5 * (node->bbox[0][0] + node->bbox[1][0]);
    double yc = 0.5 * (node->bbox[0][1] + node->bbox[1][1]);
    double zc = 0.5 * (node->bbox[0][2] + node->bbox[1][2]);

    // Limits of Node
    double min[3] = {node->bbox[0][0], node->bbox[0][1], node->bbox[0][2]};
    double max[3] = {node->bbox[1][0], node->bbox[1][1], node->bbox[1][2]};

    // Adjusting limits based on octant
    if (octantIndex & 4) { 
        min[0] = xc;
    } else {
        max[0] = xc;
    }

    if (octantIndex & 2) { 
        min[1] = yc;
    } else {
        max[1] = yc;
    }

    if (octantIndex & 1) { 
        min[2] = zc;
    } else {
        max[2] = zc;
    }

    return Init_Node(min, max);
}

int OctantIndex(Node* node, double *bd) {

    double xc = 0.5 * (node->bbox[0][0] + node->bbox[1][0]);
    double yc = 0.5 * (node->bbox[0][1] + node->bbox[1][1]);
    double zc = 0.5 * (node->bbox[0][2] + node->bbox[1][2]);

    // Determine Octant
    double x = bd[0];
    double y = bd[1];
    double z = bd[2];
    if (x < xc) {
        if (y < yc) {
            if (z < zc) {
                return 0; // Octant -x, -y, -z
            } else {
                return 1; // Octant -x, -y, +z
            }
        } else {
            if (z < zc) {
                return 2; // Octant -x, +y, -z
            } else {
                return 3; // Octant -x, +y, +z
            }
        }
    } else {
        if (y < yc) {
            if (z < zc) {
                return 4; // Octant +x, -y, -z
            } else {
                return 5; // Octant +x, -y, +z
            }
        } else {
            if (z < zc) {
                return 6; // Octant +x, +y, -z
            } else {
                return 7; // Octant +x, +y, +z
            }
        }
    }
    return NULL;
}

void updateMassCenterOfMass(Node * node, double * bdr, double * bdm){
    //Average CM  
    //(cm * M * + r * m) / (M +m)
    double mass = *bdm;
    double Mmass = *node->totalMass*mass;
    node->CoM[0] = (node->CoM[0] * *node->totalMass + bdr[0] * mass) / Mmass;
    node->CoM[1] = (node->CoM[1] * *node->totalMass + bdr[1] * mass) / Mmass;
    node->CoM[2] = (node->CoM[2] * *node->totalMass + bdr[2] * mass) / Mmass;

    *node->totalMass += mass;
}

void subdivideNode(Node* node) {
    // Create child nodes
    for (int ii = 0; ii < 8; ii++) {
        if (node->children[ii] == NULL) {
            node->children[ii] = createChildNodeForOctant(node, ii);
        }
    }
}

void Insertbd(Node * node, double * bdr, double * bdm, int ii) {

    // Empty Node
    if (node->bds == NULL && *node->totalMass == 0) {
        node->bds = ii;
        updateMassCenterOfMass(node, &bdr[3*ii], &bdm[ii]);
        return;
    }

    //If the node is a leaf but already contains a body, subdivide the node.
    if (node->bds != NULL && node->totalMass != 0) {
        subdivideNode(node);
    }

    // Insert Child
    int newbodyIndex = OctantIndex(node, &bdr[3*ii]);
    if (node->children[newbodyIndex] == NULL) {
        node->children[newbodyIndex] = createChildNodeForOctant(node, newbodyIndex);
    }


    Insertbd(node->children[newbodyIndex], bdr, bdm, ii);

    // Update the parent node Center of Mass
    updateMassCenterOfMass(node, &bdr[3*ii], &bdm[ii]);
}

void BuiltTree(double * bdr, double * bdm, int N, double * rootmin, double * rootmax){
    Node* rootNode = Init_Node(&rootmin, &rootmax);

     // Initial insertion of bodys into the tree
    int ii;
    for (ii = 0; ii < N; ++ii) {
        Insertbd(rootNode, bdr, bdm, ii);
    }
    return 0;
}

void freeNode(Node* node) {
    if (node == NULL) return;
    
    // Free Memory
    if (node->CoM != NULL) free(node->CoM);
    if (node->totalMass != NULL) free(node->totalMass);
    if (node->bbox[0] != NULL) free(node->bbox[0]);
    if (node->bbox[1] != NULL) free(node->bbox[1]);
    
    for (int ii = 0; ii < 8; ii++) {
        if (node->children[ii] != NULL) {
            freeNode(node->children[ii]); // Recursivity
        }
    }
    
    free(node);
}