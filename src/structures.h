#ifndef STRUCTURES_H_
#define STRUCTURES_H_

typedef struct {
    double *m;
    double *r;
    double *v;
    double *a;
    double *rtemp;
    double *vtemp;
} body;

typedef struct Node{
    double  *CoM;
    double  *totalMass;
    double  *bbox[2];
    int     *bds;
    Node    *children[8];
} Node;

#endif // STRUCTURES_H_