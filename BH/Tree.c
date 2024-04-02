#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include <time.h>

const int N = 100;
const double theta = 0.5;
//const double M_PI = 3.14159216;

const double G = 4 * 3.14159216 * 3.14159216; 
const double dt = 0.001; 
const int nSteps = 1;
const double epsilon = 1e-7;
int totalBodiesInserted = 0;


typedef struct body{
    double *m;
    double *r;
    double *v;
    double *a;
} body;

typedef struct Node{
    double *centerOfMass;
    double *totalMass;
    double *bbox[2];
    struct Node *children[8];
    body* bd;
    //int *index;
} Node;


Node* Init_Node(double* min, double* max);
Node* createChildNodeForOctant(Node* node, int octantIndex);

int OctantIndex(Node* node, body* bd);
void Init_bd(body* bd, int N, double* min, double* max);

void updateMassCenterOfMass(Node* nd, body* bd);
void subdivideNode(Node* nd);
void Insertbd(Node* node, body* bd);

void freeNode(Node* node);
void freeBody(body* bd, int N);

double randomDouble(double min, double max);
void writeNodeToFile(Node *node, FILE *file);
void writePositionsToCSV(body* bd, int numbd, const char* filename);


int main() {

    body bd[N];

    double *force = (double *)malloc(3 * N * sizeof(double));
    double rootMin[3] = {-70, -70, -70};
    double rootMax[3] = { 70,  70, 70};
    //double centerg[3] = {0,0,0};

    // Init quantities
    Init_bd(bd, N, rootMin, rootMax);
    Node* rootNode = Init_Node(rootMin, rootMax);

     // Initial insertion of bodys into the tree
    for (int ii = 0; ii < N; ++ii) {
        Insertbd(rootNode, &bd[ii]);
    }

    char filename[nSteps];
    sprintf(filename, "positions_%d.csv",nSteps);
    writePositionsToCSV(bd, N, filename);

    FILE *file = fopen("octants.txt", "w"); // Abre el archivo para escritura
    writeNodeToFile(rootNode, file);
    fclose(file);

    //calcular fuerzas
 
    freeBody(bd, N);
    free(force);
    freeNode(rootNode);

    return 0;
}


Node* Init_Node(double *min, double *max){


    Node* node = (Node *)malloc(sizeof(Node));
    if (node == NULL) {
        fprintf(stderr, "Error Allocation Memory Init_Node\n");
        exit(EXIT_FAILURE);
    }

    //Memory allocation
    node->centerOfMass = (double *)malloc(3 * sizeof(double));  // x, y, z
    node->totalMass    = (double *)malloc(sizeof(double));      // mass
    node->bbox[0] = (double *)malloc(3 * sizeof(double));       // min[x,y,z]
    node->bbox[1] = (double *)malloc(3 * sizeof(double));       // max[x,y,z]

    if (node->centerOfMass == NULL || node->totalMass == NULL || node->bbox[0] == NULL || node->bbox[1] == NULL) {
        fprintf(stderr, "Error Allocation Memory in Init_Node\n");
        exit(EXIT_FAILURE);
    }

    // Initialization
    for(int ii = 0; ii < 3; ++ii) {
        node->bbox[0][ii] = min[ii]; 
        node->bbox[1][ii] = max[ii];
        node->centerOfMass[ii] = 0.0;
    }
    *node->totalMass = 0.0; 
    
    for(int ii = 0; ii < 8; ++ii) {
        node->children[ii] = NULL;
    }
    node->bd = NULL;

    return node;
}


void Init_bd(body *bd, int N, double *min, double *max) {
    // seed
    srand(time(NULL));

    for (int ii = 0; ii < N; ++ii) {
        bd[ii].r = (double *)malloc(3 * sizeof(double)); // x, y, z
        bd[ii].v = (double *)malloc(3 * sizeof(double)); // vx vy vz
        bd[ii].a = (double *)malloc(3 * sizeof(double)); // ax ay az
        bd[ii].m = (double *)malloc(sizeof(double));     // m

        if (bd[ii].r == NULL || bd[ii].v == NULL || bd[ii].a == NULL || bd[ii].m == NULL) {
            fprintf(stderr, "Error Allocation Memory Init_bd\n");
            exit(EXIT_FAILURE);
        }

        // random pos in domain 
        bd[ii].r[0] = randomDouble(min[0] + 20, max[0] - 20);  // x component
        bd[ii].r[1] = randomDouble(min[1] + 20, max[1] - 20);  // y component
        bd[ii].r[2] = randomDouble(min[2] + 20, max[2] - 20);  // z component

        bd[ii].v[0] = 0.0; 
        bd[ii].v[1] = 0.0; 
        bd[ii].v[2] = 0.0; 

        bd[ii].a[0] = 0.0; 
        bd[ii].a[1] = 0.0;
        bd[ii].a[2] = 0.0; 

        // mass 
        *bd[ii].m = 1.0; 
    }
}


void Insertbd(Node* node, body* bd) {

    // Empty Node
    if (node->bd == NULL && *node->totalMass == 0) {
        node->bd = bd;
        updateMassCenterOfMass(node, bd);
        return;
    }

    //If the node is a leaf but already contains a body, subdivide the node.
    if (node->bd != NULL && node->totalMass != 0) {
        subdivideNode(node);  
    }

    // Insert Child
    int newbodyIndex = OctantIndex(node, bd);
    if (node->children[newbodyIndex] == NULL) {
        node->children[newbodyIndex] = createChildNodeForOctant(node, newbodyIndex);

    }

    Insertbd(node->children[newbodyIndex], bd);

    // freeNode() liberar memoria del nodo ya que los hijos son los que tienen la informacin de los cuerpos
    // Update the parent node Center of Mass
    updateMassCenterOfMass(node, bd);
    
}


void updateMassCenterOfMass(Node* node, body* bd){
    if (*node->totalMass == 0){
        //Empty node
        node->centerOfMass[0] = bd->r[0]; 
        node->centerOfMass[1] = bd->r[1]; 
        node->centerOfMass[2] = bd->r[2]; 
    }
    else { //Average CM  
           //(cm * M * + r * m) / (M +m)

        node->centerOfMass[0] = (node->centerOfMass[0] * *node->totalMass + bd->r[0] * *bd->m) / (*node->totalMass + *bd->m);
        node->centerOfMass[1] = (node->centerOfMass[1] * *node->totalMass + bd->r[1] * *bd->m) / (*node->totalMass + *bd->m);
        node->centerOfMass[2] = (node->centerOfMass[2] * *node->totalMass + bd->r[2] * *bd->m) / (*node->totalMass + *bd->m);
    }
    *node->totalMass += *bd->m;
}


void subdivideNode(Node* node) {
    // Create child nodes
    for (int ii = 0; ii < 8; ii++) {
        // Node not Full
        if (node->children[ii] == NULL) {

            node->children[ii] = createChildNodeForOctant(node, ii);
        }
    }
}


int OctantIndex(Node* node, body* bd) {
    // Calculate center of node
    double center[3];

    center[0] = (node->bbox[0][0] + node->bbox[1][0]) / 2.0;
    center[1] = (node->bbox[0][1] + node->bbox[1][1]) / 2.0;
    center[2] = (node->bbox[0][2] + node->bbox[1][2]) / 2.0;

    int index = 0;

    // Determine Octant
    if (bd->r[0] < center[0]) {
        if (bd->r[1] < center[1]) {
            if (bd->r[2] < center[2]) {
                return 0; // Octant -x, -y, -z
            } else {
                return 1; // Octant -x, -y, +z
            }
        } else {
            if (bd->r[2] < center[2]) {
                return 2; // Octant -x, +y, -z
            } else {
                return 3; // Octant -x, +y, +z
            }
        }
    } else {
        if (bd->r[1] < center[1]) {
            if (bd->r[2] < center[2]) {
                return 4; // Octant +x, -y, -z
            } else {
                return 5; // Octant +x, -y, +z
            }
        } else {
            if (bd->r[2] < center[2]) {
                return 6; // Octant +x, +y, -z
            } else {
                return 7; // Octant +x, +y, +z
            }
        }
    }

    
    printf("Índice de octante calculado: %d\n", index);
    return index;
}
 

Node* createChildNodeForOctant(Node* node, int octantIndex) {

    // Center of Node
    double center[3]; 
    center[0] = (node->bbox[0][0] + node->bbox[1][0]) / 2.0;
    center[1] = (node->bbox[0][1] + node->bbox[1][1]) / 2.0;
    center[2] = (node->bbox[0][2] + node->bbox[1][2]) / 2.0;

    // Limits of Node
    double min[3] = {node->bbox[0][0], node->bbox[0][1], node->bbox[0][2]};
    double max[3] = {node->bbox[1][0], node->bbox[1][1], node->bbox[1][2]};

    // Adjusting limits based on octant
    if (octantIndex & 4) { 
        min[0] = center[0];
    } else {
        max[0] = center[0];
    }

    if (octantIndex & 2) { 
        min[1] = center[1];
    } else {
        max[1] = center[1];
    }

    if (octantIndex & 1) { 
        min[2] = center[2];
    } else {
        max[2] = center[2];
    }


    return Init_Node(min, max);
}


void freeBody(body* bd, int N) {

    for (int ii = 0; ii < N; ++ii) {
        if (bd[ii].r != NULL) {
            free(bd[ii].r);
        }
        if (bd[ii].v != NULL) {
            free(bd[ii].v);
        }
        if (bd[ii].a != NULL) {
            free(bd[ii].a);
        }
        if (bd[ii].m != NULL) {
            free(bd[ii].m);
        }
    }
}


void freeNode(Node* node) {
    if (node == NULL) return;
    
    // Free Memory
    if (node->centerOfMass != NULL) free(node->centerOfMass);
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


double randomDouble(double min, double max) {
    return min + (rand() / (RAND_MAX / (max - min)));
}


void writeNodeToFile(Node *node, FILE *file) {
    if (node == NULL) return;
    
    // write limits of nodes
    fprintf(file, "%f,%f,%f,%f,%f,%f\n",
            node->bbox[0][0], node->bbox[0][1], node->bbox[0][2],
            node->bbox[1][0], node->bbox[1][1], node->bbox[1][2]);
    
    // Recursivity for childs
    for (int i = 0; i < 8; i++) {
        writeNodeToFile(node->children[i], file);
    }
}


void writePositionsToCSV(body* bd, int numbd, const char* filename) {
    FILE* file = fopen(filename, "w"); 
    if (file == NULL) {
        printf("Error to open file.\n");
        return;
    }

    fprintf(file, "X,Y,Z\n");

    for (int ii = 0; ii < numbd; ii++) {
        fprintf(file, "%f,%f,%f\n", bd[ii].r[0], bd[ii].r[1], bd[ii].r[2]);
    }

    fclose(file);
}
