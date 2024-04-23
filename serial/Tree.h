#ifndef TREE_H_
#define TREE_H_
#include "structures.h"

Node* RootNode(const double *, const double *, const double *, const double *, const int);

void CreateNode(Node *, Node *, double *, double *);

void createChildren(Node*, const int);

void Insertbd(Node *, const double *, const double *, const int);

void updateCenterOfMass(Node *);

void DivideNode(Node *, const double *, const double *, const double);

Node* nextnode(Node *, int);

Node* BuiltTree(const double *, const double *, const int, const double *, const double *, const double);

void freeNode(Node*);

#endif // TREE_H_