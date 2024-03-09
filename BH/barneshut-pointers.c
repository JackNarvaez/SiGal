#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include <time.h>
//#include <mpi.h>


const int N = 100;
const double theta = 0.5;
//const double M_PI = 3.14159216;
int totalBodiesInserted = 0;

const double G = 4 * 3.14159216 * 3.14159216; 
const double dt = 0.001; 
const int nSteps = 100;
const double epsilon = 1e-7;


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
    body *bd;
} Node;


Node* Init_Node(double* min, double* max);
Node* createChildNodeForOctant(Node* node, int octantIndex);

int OctantIndex(Node* node, body* bd);
void Init_bd(body* bd, int N, double* min, double* max);
void updateMassCenterOfMass(Node* nd, body* bd);
void subdivideNode(Node* nd);
void Insertbd(Node* node, body* bd);

double distance(double* ri, double* rj);
void CenterOfMass(Node* node);
void Force(Node* node, body* bd, double* f);
void Euler(body* bd, double* force, double dt);

void freeNode(Node* node);
void freeBody(body* bd, int N);

double randomDouble(double min, double max);
void writePositionsToCSV(body* bd, int numbd, const char* filename);
void Init_bd2(body* bd, int N, double* center);



int main() {

    body bd[N];
    srand(time(NULL));

    double *force = (double *)malloc(3 * N * sizeof(double));
    double rootMin[3] = {-70, -70, -70};
    double rootMax[3] = { 70,  70, 70};
    double centerg[3] = {0,0,0};

    // Init quantities
    //Init_bodys2(bd, N, centerg);
    Init_bd(bd, N, rootMin, rootMax);
    Node* rootNode = Init_Node(rootMin, rootMax);

    // Initial insertion of bodys into the tree
    for (int ii = 0; ii < N; ++ii) {
        Insertbd(rootNode, &bd[ii]);
    }

    char filename[nSteps];

    for (int step = 0; step < nSteps; step++) {
        // Write positions
        if (step % 1 == 0)
        {
            sprintf(filename, "positions_%d.csv", step);
            writePositionsToCSV(bd, N, filename);
        }

        // clean forces
        memset(force, 0, sizeof(force));

        // calculate force of each body
        for (int ii = 0; ii < N; ii++) {

            double f[3] = {0.0, 0.0, 0.0}; 
            Force(rootNode, &bd[ii], f);
            force[3*ii] = f[0];
            force[3*ii + 1] = f[1];
            force[3*ii + 2] = f[2];
        }

        // Update vel and acc for each particle
        for (int ii = 0; ii < N; ii++) {
            Euler(&bd[ii], &force[3*ii], dt);  // Pasamos la dirección del segmento de fuerza correspondiente a la partícula ii
        }


        // Clean last node
        freeNode(rootNode);
        rootNode = Init_Node(rootMin, rootMax);

        for (int ii = 0; ii < N; ii++) {
            Insertbd(rootNode, &bd[ii]); 
        }

        //CenterOfMass(rootNode);

    }

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
    //srand(time(NULL));

    for (int ii = 0; ii < N; ++ii) {
        bd[ii].r = (double *)malloc(3 * sizeof(double)); // x, y, z
        bd[ii].v = (double *)malloc(3 * sizeof(double)); // vx vy vz
        bd[ii].a = (double *)malloc(3 * sizeof(double)); // ax ay az
        bd[ii].m = (double *)malloc(sizeof(double));     // m

        if (bd[ii].r == NULL || bd[ii].v == NULL || bd[ii].a == NULL || bd[ii].m == NULL) {
            fprintf(stderr, "Error Allocation Memory Init_bd\n");
            exit(EXIT_FAILURE);
        }

        // random pos in domain para este cuerpo específico
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


void Init_bd2(body* bd, int N, double *center) {

    double galaxy_radius = 50.0;
    double disk_thickness = 1.5;
    double galaxy_mass = N; 

    for (int ii = 0; ii < N; ++ii) {
        bd[ii].r = (double *)malloc(3 * sizeof(double));
        bd[ii].v = (double *)malloc(3 * sizeof(double));
        bd[ii].a = (double *)malloc(3 * sizeof(double));
        bd[ii].m = (double *)malloc(sizeof(double));

        if (bd[ii].r == NULL || bd[ii].v == NULL || bd[ii].a == NULL || bd[ii].m == NULL) {
            fprintf(stderr, "Error Allocation Memory Init_bd2\n");
            exit(EXIT_FAILURE);
        }
    }

    for (int ii = 0; ii < N; ii++) {
        // Radial distribution

        double radius = ((double)rand() / RAND_MAX) * galaxy_radius;
        double theta = ((double)rand() / RAND_MAX) * 2 * 3.141592; // Ángulo aleatorio en radianes

        // Height to disk
        double height = ((double)rand() / RAND_MAX) * disk_thickness - disk_thickness / 2;

        // coordinates
        bd[ii].r[0] = center[0] + radius * cos(theta); // x
        bd[ii].r[1] = center[1] + radius * sin(theta); // y
        bd[ii].r[2] = center[2] + height;              // z

        // initial velocity: orbital velocity in x-y plane
        double orbital_velocity = sqrt((G * galaxy_mass) / radius); // Fórmula simplificada
        bd[ii].v[0] = -orbital_velocity * sin(theta); // vx
        bd[ii].v[1] = orbital_velocity * cos(theta);  // vy
        bd[ii].v[2] = 0.0;                            // vz

        bd[ii].a[0] = 0.0; 
        bd[ii].a[1] = 0.0; 
        bd[ii].a[2] = 0.0; 

        // Mass
        *bd[ii].m = 1.0; // Asumiendo que cada cuerpo tiene una masa unitaria
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

    // Update the parent node Center of Mass
    updateMassCenterOfMass(node, bd);
}


void updateMassCenterOfMass(Node* node, body* bd){
    if (node->totalMass == 0){
        //Empty node
        node->centerOfMass[0] = bd->r[0]; 
        node->centerOfMass[1] = bd->r[1]; 
        node->centerOfMass[2] = bd->r[2]; 
    }
    else { //Average CM
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
    // Calcula el centro del nodo
    double center[3];

    center[0] = (node->bbox[0][0] + node->bbox[1][0]) / 2.0;
    center[1] = (node->bbox[0][1] + node->bbox[1][1]) / 2.0;
    center[2] = (node->bbox[0][2] + node->bbox[1][2]) / 2.0;

    int index = 0;

    // X-axis
    if (bd->r[0] >= center[0]) {
        index += 4;
    }
    // Y-xis
    if (bd->r[1] >= center[1]) {
        index += 2;
    }
    // Z-axis
    if (bd->r[2] >= center[2]) {
        index += 1;
    }

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

    // Octant index
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


double distance(double* ri, double* rj)
{
    return sqrt( (rj[0] - ri[0]) * (rj[0] - ri[0]) 
               + (rj[1] - ri[1]) * (rj[1] - ri[1]) 
               + (rj[2] - ri[2]) * (rj[2] - ri[2]) 
               );
}


void CenterOfMass(Node* node) {
    double totalMass = 0.0;
    double weightedPositionSum[3] = {0.0, 0.0, 0.0};

    // Si el nodo tiene un cuerpo directamente, sumar su masa y su contribución al centro de masa
    if (node->bd != NULL) {
        totalMass += *node->bd->m;
        weightedPositionSum[0] += node->bd->r[0] * *node->bd->m;
        weightedPositionSum[1] += node->bd->r[1] * *node->bd->m;
        weightedPositionSum[2] += node->bd->r[2] * *node->bd->m;
    }

    // Recorrer los nodos hijos para acumular su masa y contribución al centro de masa
    for (int i = 0; i < 8; ++i) {
        if (node->children[i] != NULL) {
            CenterOfMass(node->children[i]); // Llamada recursiva
            totalMass += *node->children[i]->totalMass;
            weightedPositionSum[0] += node->children[i]->centerOfMass[0] * *node->children[i]->totalMass;
            weightedPositionSum[1] += node->children[i]->centerOfMass[1] * *node->children[i]->totalMass;
            weightedPositionSum[2] += node->children[i]->centerOfMass[2] * *node->children[i]->totalMass;
        }
    }

    // Si hay masa total, calcular el centro de masa del nodo
    if (totalMass > 0) {
        node->centerOfMass[0] = weightedPositionSum[0] / totalMass;
        node->centerOfMass[1] = weightedPositionSum[1] / totalMass;
        node->centerOfMass[2] = weightedPositionSum[2] / totalMass;
    }

    *node->totalMass = totalMass; // Actualizar la masa total del nodo
}


void Force(Node* node, body* bd, double* f ){
    // Emty node
    if (*node->totalMass == 0) return;
    
    double dir[3];
    double r_ij = distance(bd->r, node->centerOfMass);

    // Avoid divergence
    if (r_ij == epsilon) return;

    // size of node
    double nodeSize = node->bbox[1][0] - node->bbox[0][0];
    double h = nodeSize / r_ij;

    // If the node is far enough away or is a leaf, treat the node as a point mass
    if (h < theta || node->children[0] == NULL ){

        dir[0] = node->centerOfMass[0] - bd->r[0];
        dir[1] = node->centerOfMass[1] - bd->r[1];
        dir[2] = node->centerOfMass[2] - bd->r[2];

        double F = -G * *bd->m * *node->totalMass / (r_ij*r_ij*r_ij);
        f[0] += F*dir[0];
        f[1] += F*dir[1];
        f[2] += F*dir[2];
    }
    else{
        // Otherwise calculate force recursively on subnodes
        for (int ii = 0; ii < 8; ii++){
            if (node->children[ii] != NULL){
                Force(node->children[ii], bd, f);
            }    
        }   
    }

}


void Euler(body* bd, double* force, double dt) {
    // Calculate Aceleration

    bd->a[0] = force[0] / *bd->m;
    bd->a[1] = force[1] / *bd->m;
    bd->a[2] = force[2] / *bd->m;

    // Update Velocity
    bd->v[0] += 0.5 * bd->a[0] * dt;
    bd->v[1] += 0.5 * bd->a[1] * dt;
    bd->v[2] += 0.5 * bd->a[2] * dt;

    // Update Positions
    bd->r[0] += 0.5 * bd->v[0] * dt;
    bd->r[1] += 0.5 * bd->v[1] * dt;
    bd->r[2] += 0.5 * bd->v[2] * dt;
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
    
    // Liberar memoria de punteros dentro del nodo
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


void writePositionsToCSV(body* bd, int numbd, const char* filename) {
    FILE* file = fopen(filename, "w"); 
    if (file == NULL) {
        printf("Error al abrir el archivo.\n");
        return;
    }

    fprintf(file, "X,Y,Z\n");

    for (int ii = 0; ii < numbd; ii++) {
        fprintf(file, "%f,%f,%f\n", bd[ii].r[0], bd[ii].r[1], bd[ii].r[2]);
    }

    fclose(file);
}

double randomDouble(double min, double max) {
    return min + (rand() / (RAND_MAX / (max - min)));
}



