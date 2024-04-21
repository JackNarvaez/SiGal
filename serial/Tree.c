#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include "structures.h"


Node* RootNode(const double *min, const double *max, const double * bdr, const double * bdm, const int N){
    
    //Memory allocation
    Node* node = (Node *)malloc(sizeof(Node));
    if (node == NULL) {
        fprintf(stderr, "Error Allocation Memory RootNode.\n");
        exit(EXIT_FAILURE);
    }

    node->CoM   = (double *)malloc(3 * sizeof(double)); // Center of mass
    node->min   = (double *)malloc(3 * sizeof(double)); // Low-left corner
    node->max   = (double *)malloc(3 * sizeof(double)); // Upper-right corner
    node->Mass  = (double *)malloc(sizeof(double));     // Mass
    node->deep  = (int *)malloc(sizeof(int));           // deep
    node->slice = (int *)malloc(sizeof(int));           // Slicing addresses
    node->bodies= (int *)malloc((N+1)*sizeof(int));     // Particles in this node
    node->child = (Node *)malloc(sizeof(Node));         // Child Node
    node->sibling= (Node *)malloc(sizeof(Node));        // Sibling Node

    if (node->CoM == NULL || node->Mass == NULL || node->deep == NULL || 
        node->min == NULL || node->max == NULL || node->slice == NULL || node->bodies == NULL ||
        node->child == NULL || node->sibling == NULL) {
        fprintf(stderr, "Error Allocation Memory in RootNode's Features.\n");
        exit(EXIT_FAILURE);
    }

    int ii;

    // Initialize arrays
    for (ii = 0; ii < 3; ++ii) {
        node->min[ii] = min[ii]; 
        node->max[ii] = max[ii];
        node->CoM[ii] = 0.0;
    }
    *node->Mass = 0.0;
    *node->deep = 0;
    node->bodies[0] = N;

    // Calculate Total Mass and CoM
    for (ii = 0; ii < N; ii++) {
        node->bodies[ii+1] = ii;
        *node->Mass  += bdm[ii];
        node->CoM[0] += bdm[ii]*bdr[3*ii];
        node->CoM[1] += bdm[ii]*bdr[3*ii+1];
        node->CoM[2] += bdm[ii]*bdr[3*ii+2];
    }
    double invMass = 1.0 / *node->Mass;
    node->CoM[0]    *= invMass;
    node->CoM[1]    *= invMass;
    node->CoM[2]    *= invMass;

    *node->slice    = 0;
    node->type      = (N>1)? false:true;
    node->next      = NULL;

    return node;
}

void CreateNode(Node* node, Node* father, double *rmin, double *rmax){

    //Memory allocation
    node->CoM   = (double *)malloc(3 * sizeof(double)); // Center of mass
    node->min   = (double *)malloc(3 * sizeof(double)); // Low-left corner
    node->max   = (double *)malloc(3 * sizeof(double)); // Upper-right corner
    node->Mass  = (double *)malloc(sizeof(double));     // Mass
    node->deep  = (int *)malloc(sizeof(int));           // deep
    node->slice = (int *)malloc(sizeof(int));           // Slicing addresses
    node->bodies= (int *)malloc((father->bodies[0]+1)*sizeof(int));         // Particles in this node
    node->child = (Node *)malloc(sizeof(Node));         // Child Node
    node->sibling= (Node *)malloc(sizeof(Node));        // Sibling Node

    if (node->CoM == NULL || node->Mass == NULL || node->deep == NULL || 
        node->min == NULL || node->max == NULL || node->slice == NULL || node->bodies == NULL ||
        node->child == NULL || node->sibling == NULL) {
        fprintf(stderr, "Error Allocation Memory in Node.\n");
        exit(EXIT_FAILURE);
    }

    int ii;

    // Initialize arrays
    for (ii = 0; ii < 3; ++ii) {
        node->min[ii] = rmin[ii]; 
        node->max[ii] = rmax[ii];
        node->CoM[ii] = 0.0;
    }
    *node->Mass     = 0.0;
    *node->deep     = *father->deep+1;
    *node->slice    = 0;
    node->bodies[0] = 0;

    node->type      = true;    
    node->next      = NULL; 
}

void createChildren(Node* node, const int address) {

    *node->slice = address;

    double min[3]   = {node->min[0], node->min[1], node->min[2]};
    double max[3]   = {node->max[0], node->max[1], node->max[2]};

    switch (address)
    {
    case 0:
    {
        // X slicing
        double med[3] = {0.5*(node->max[0]+node->min[0]), node->max[1], node->max[2]};
        CreateNode(node->child, node, min, med);
        med[1] = node->min[1];
        med[2] = node->min[2];
        CreateNode(node->sibling, node, med, max);
        break;
    }
    case 1:
    {

        // Y slicing
        double med[3] = {node->max[0], 0.5*(node->max[1]+node->min[1]), node->max[2]};
        CreateNode(node->child, node, min, med);
        med[0] = node->min[0];
        med[2] = node->min[2];
        CreateNode(node->sibling, node, med, max);
        break;
    }
    case 2:
    {
        // Z slicing
        double med[3] = {node->max[0], node->max[1], 0.5*(node->max[2]+node->min[2])};
        CreateNode(node->child, node, min, med);
        med[0] = node->min[0];
        med[1] = node->min[1];
        CreateNode(node->sibling, node, med, max);
        break;
    }
    default:
        fprintf(stderr, "Error creating children.\n");
        exit(EXIT_FAILURE);
        break;
    }
    (node->child)->next = node->sibling;
    (node->sibling)->next = node;
}

void Insertbd(Node * node, const double * bdr, const double * bdm, const int ii) {
    node->bodies[0] ++;
    node->bodies[node->bodies[0]] = ii;
    *node->Mass += *bdm;
    node->CoM[0] += *bdm*bdr[0]; 
    node->CoM[1] += *bdm*bdr[1]; 
    node->CoM[2] += *bdm*bdr[2]; 
}

void updateCenterOfMass(Node * node){
    double invMass = 1.0 / *node->Mass;
    node->CoM[0] = invMass;
    node->CoM[1] = invMass;
    node->CoM[2] = invMass;
}

void DivideNode(Node* node, const double * bdr, const double * bdm) {
    // Find address for slicing
    double side[3] = {node->max[0]-node->min[0], node->max[1]-node->min[1], node->max[2]-node->min[2]};
    int slicingAddress = 0;
    if (side[1]>side[slicingAddress]) slicingAddress = 1;
    if (side[2]>side[slicingAddress]) slicingAddress = 2;

    // Slice
    createChildren(node, slicingAddress);

    // Update
    int ii;
    bool branch;
    for (ii=0; ii<node->bodies[0]; ii++) {
        branch = (bdr[node->bodies[ii]+*node->slice] < (node->child)->max[*node->slice]) ? true:false;
        if (branch) {
            Insertbd(node->child, &bdr[3*ii], &bdm[ii], ii);
        } else {
            Insertbd(node->sibling, &bdr[3*ii], &bdm[ii], ii);
        }
    }
    updateCenterOfMass(node->child); 
    updateCenterOfMass(node->sibling);

    if ((node->child)->bodies[0]>1) {
        (node->child)->type = false;
    } else {
        free((node->child)->bodies);
        free((node->child)->child);
        free((node->child)->sibling);
    } 
    if ((node->sibling)->bodies[0]>1) {
        (node->sibling)->type = false;
    } else {
        free((node->sibling)->bodies);
        free((node->sibling)->child);
        free((node->sibling)->sibling);
    }
    free(node->bodies);
}

Node *nextnode(Node* node, int sense) {
    if (!node->type) {
        if (sense==0){
            return node->child;
        } else{
            return node->next;
        }
    } else{
        return node->next;
    }
}

Node* BuiltTree(const double * bdr, const double * bdm, const int N, const double * rootmin, const double * rootmax){
    Node* rootNode = RootNode(rootmin, rootmax, bdr, bdm, N);
    DivideNode(rootNode, bdr, bdm);
    int pstnd   = *rootNode->deep;
    Node* Next  = rootNode->child;
    while (Next != rootNode){
        if (!Next->type && pstnd<=*Next->deep){
            DivideNode(Next, bdr, bdm);
        }
        bool sense = pstnd<=*Next->deep;
        pstnd   = *Next->deep;
        if (sense) {
            Next    = nextnode(Next, 0);
        } else {
            Next    = nextnode(Next, 1);
        }
    }
    return rootNode;
}

void freeNode(Node * node){
    int pstnd   = *node->deep;
    Node* Next  = node->child;
    while (Next != node){
        bool sense = pstnd>*Next->deep;
        pstnd   = *Next->deep;
        if (Next->type){
            free(Next->CoM);
            free(Next->min);
            free(Next->max);
            free(Next->Mass);
            free(Next->deep);
            free(Next->slice);
        } else {
            if (sense) {
                free(Next->CoM);
                free(Next->min);
                free(Next->max);
                free(Next->Mass);
                free(Next->deep);
                free(Next->slice);
                free(Next->child);
                free(Next->sibling);
            }
        }
        if (!sense) {
            Next    = nextnode(Next, 0);
        } else {
            Next    = nextnode(Next, 1);
        }
    }
    free(node->CoM);
    free(node->min);
    free(node->max);
    free(node->Mass);
    free(node->deep);
    free(node->slice);
    free(node->child);
    free(node->sibling);
    free(node);
}