#ifndef STRUCTURES_H_
#define STRUCTURES_H_
#include <stdbool.h>

typedef struct {
    double *m;
    double *r;
    double *v;
    double *a;
    int    *i;
} body;

typedef struct Node{
    double  *CoM;
    double  *min;
    double  *max;
    double  *Mass;
    int     *deep;
    int     *slice;     // 0: x; 1: y; 2: z;
    bool    type;       // False: Node; True: Leaf
    int     *bodies;    // 0th-position stores the total number of bodies
    struct  Node* child;
    struct  Node* sibling;
    struct  Node* next;
} Node;

#endif // STRUCTURES_H_