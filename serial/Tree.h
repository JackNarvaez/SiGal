#ifndef TREE_H_
#define TREE_H_
#include "structures.h"

Node* RootNode(const double *min, const double *max, const double * bdr, const double * bdm, const int N);

Node* CreateNode(Node* father, double *rmin, double *rmax);

void createChildren(Node* node, const int address);

void Insertbd(Node * node, const double * bdr, const double * bdm, const int ii) ;

void updateCenterOfMass(Node * node);

void DivideNode(Node* node, const double * bdr, const double * bdm);

Node* nextnode(Node* node, int pstnode);

Node* BuiltTree(const double * bdr, const double * bdm, const int N, const double * rootmin, const double * rootmax);

void freeNode(Node* node);

#endif // TREE_H_