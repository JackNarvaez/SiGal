#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include <time.h>
//#include <mpi.h>

#define RAND() ((double)rand()/(double)(RAND_MAX))

const int N = 1000;
const double theta = 0.5;

const double G = 4 * (3.14159216)*(3.14159216); 
const double dt = 0.001; 
const int nSteps = 100000;
const double epsilon = 1e-7;

const double TPI = 2*3.14159265358979323846;
double sqrt2; 
double scaleftr; 
double invsclftr; 
double sqrtscldrt;

typedef struct Vec3 {
    double x, y, z;
} Vec3;

typedef struct Particle{
    Vec3 position;
    Vec3 velocity;
    double mass;
} Particle;

typedef struct Node {
    Vec3 centerOfMass;
    double totalMass;
    Vec3 bbox[2];
    struct Node* children[8];
    Particle* particle;
} Node;


Node* Init_Node(Vec3 min, Vec3 max);
Node* createChildNodeForOctant(Node* parent, int octantIndex);
int OctantIndex(Node* node, Particle* particle);
void Init_Particles(Particle* particles, int N, Vec3 min, Vec3 max);
void updateMassCenterOfMass(Node* node, Particle* particle);
void subdivideNode(Node* node);
void insertParticle(Node* node, Particle* particle);
void CenterOfMass(Node* node);
void Force(Node* node, Particle* particle, Vec3* f);
void Euler(Particle* particle, Vec3 force, double dt);
void freeNode(Node* node);
double randomDouble(double min, double max);
void writePositionsToCSV(Particle* particles, int numParticles, const char* filename);
void writeNodeToFile(Node *node, FILE *file);


void Init_Particles2(Particle* particles, int N);

void spher2cartes(Vec3 * Vec, double r, double sclftr);
double g(double x);
double rndm(double min, double max);
void frm2com(Particle *particles, int N);

void adjust_units(Particle *particles, int N);
void plummer_dist(Particle *particles, int N, double sqrt2, double invsclftr, double sqrtscldrt);

int main() {

    Particle particles[N];
    Vec3 force[N];
    Vec3 rootMin = {-70, -70, -70};
    Vec3 rootMax = {70, 70, 70};

    // Init quantities
    //Init_Particles2(particles, N);
    //Init_Particles(particles, N, rootMin, rootMax);

    sqrt2 = sqrt(2.0);
    scaleftr = 16.0 / (3.0 * TPI);
    invsclftr = 1.0 / scaleftr;
    sqrtscldrt = sqrt(scaleftr);


    plummer_dist(particles, N, sqrt2, invsclftr,sqrtscldrt);
    Node* rootNode = Init_Node(rootMin, rootMax);

    // Initial insertion of particles into the tree
    for (size_t ii = 0; ii < N; ii++) {
        insertParticle(rootNode, &particles[ii]);
    }

    char filename[nSteps];
    int count =0;
    for (int step = 0; step < nSteps; step++) {
        // Write positions
        
        if (step % 100 == 0)
        {
            sprintf(filename, "/home/yo/Documents/NBodySimulations/BH/files/positions_%d.csv", count);
            writePositionsToCSV(particles, N, filename);
            count++;
        }
        

        // clean forces
        memset(force, 0, sizeof(force));

        // calculate force of each particle
        for (size_t ii = 0; ii < N; ii++) {

            Vec3 f = {0, 0, 0}; 
            Force(rootNode, &particles[ii], &f);
            force[ii] = f; 
        }

        // Update vel and acc
        for (size_t ii = 0; ii < N; ii++) {
            Euler(&particles[ii], force[ii], dt);
        }

        // Clean last node
        freeNode(rootNode);
        rootNode = Init_Node(rootMin, rootMax);

        for (size_t ii = 0; ii < N; ii++) {
            insertParticle(rootNode, &particles[ii]); 
        }

        CenterOfMass(rootNode);


    }

    FILE *file = fopen("octants.txt", "w"); 
    writeNodeToFile(rootNode, file);
    fclose(file);


    // CLEAN
    freeNode(rootNode);

    return 0;
}


Node* Init_Node(Vec3 min, Vec3 max){

    Node* node = (Node*)malloc(sizeof(Node));
    if (node == NULL)
    {
        fprintf(stderr, "Error allocation memory Node\n");
        exit(EXIT_FAILURE);
    }

    //Init node parameters
    node->bbox[0] = min; 
    node->bbox[1] = max;
    node->centerOfMass = (Vec3){0,0,0};
    node->totalMass = 0;
    
    for(int ii = 0; ii< 8 ; ++ii){
        node->children[ii]= NULL;
    }

    node->particle = NULL;

    return node;
}


Node* createChildNodeForOctant(Node* parent, int octantIndex) {
    // Calculate center point of the parent node
    Vec3 center = {
        (parent->bbox[0].x + parent->bbox[1].x) / 2,
        (parent->bbox[0].y + parent->bbox[1].y) / 2,
        (parent->bbox[0].z + parent->bbox[1].z) / 2
    };

    Vec3 min = parent->bbox[0];
    Vec3 max = parent->bbox[1];

    if (octantIndex & 1) { // If LSB is set, adjust z for upper half
        min.z = center.z;
    } else { // Lower half
        max.z = center.z;
    }

    if (octantIndex & 2) { // Second bit, adjust y for upper half
        min.y = center.y;
    } else { // Lower half
        max.y = center.y;
    }

    if (octantIndex & 4) { // Third bit, adjust x for upper half
        min.x = center.x;
    } else { // Lower half
        max.x = center.x;
    }

    return Init_Node(min, max);
}


double randomDouble(double min, double max) {
    return min + (rand() / (RAND_MAX / (max - min)));
}


void Init_Particles(Particle* particles, int N, Vec3 min, Vec3 max) {
    // seed
    srand(time(NULL));

    for (size_t ii = 0; ii < N; ii++) {
        // random pos in domain
        particles[ii].position.x = randomDouble(min.x+20, max.x-20);
        particles[ii].position.y = randomDouble(min.y+20, max.y-20);
        particles[ii].position.z = randomDouble(min.z+20, max.z-20);

        particles[ii].velocity.x = 0; 
        particles[ii].velocity.y = 0; 
        particles[ii].velocity.z = 0; 

        // mass
        particles[ii].mass = 1.0; 
    }
}



void Init_Particles2(Particle* particles, int N) {
    srand(time(NULL));

    Vec3 center = {0,0,0};
    double disk_thickness = 1.0;
    double galaxy_radius = 50.0;
    double galaxy_mass = N;

    for (size_t ii = 0; ii < N; ii++) {
        // Radial distribution
        double radius = ((double)rand() / RAND_MAX) * galaxy_radius;
        double theta = ((double)rand() / RAND_MAX) * 2 * 3.14159216; // Ángulo aleatorio en radianes

        // Random height to simulate the thickness of the galaxy's disk
        double height = ((double)rand() / RAND_MAX) * disk_thickness - disk_thickness / 2;

        // changue of coordinates
        particles[ii].position.x = center.x + radius * cos(theta);
        particles[ii].position.y = center.y + radius * sin(theta);
        particles[ii].position.z = center.z + height;

        //Init Velocity: orbital velocity
        double orbital_velocity = sqrt((G * galaxy_mass) / radius); 
        particles[ii].velocity.x = -orbital_velocity * sin(theta);
        particles[ii].velocity.y = orbital_velocity * cos(theta);
        particles[ii].velocity.z = 0; 

        // Mass
        particles[ii].mass = 1.0; 
    }
}


int OctantIndex(Node* node, Particle* particle) {
    Vec3 center = {
        (node->bbox[0].x + node->bbox[1].x) / 2,
        (node->bbox[0].y + node->bbox[1].y) / 2,
        (node->bbox[0].z + node->bbox[1].z) / 2
    };
    
    // Start with a base index of 0
    int index = 0;

    // Check the x dimension
    if (particle->position.x >= center.x) {
        index += 4; // x 
    }

    // Check the y dimension
    if (particle->position.y >= center.y) {
        index += + 2; // y 
    }

    // Check the z dimension
    if (particle->position.z >= center.z) {
        index += 1; // z 
    }

    return index;
}

void updateMassCenterOfMass(Node* node, Particle* particle){
    if (node->totalMass == 0){
        //Empty node
        node->centerOfMass = particle->position;
    }
    else { //Average CM
        node->centerOfMass.x = (node->centerOfMass.x * node->totalMass + particle->position.x * particle->mass) / (node->totalMass + particle->mass);
        node->centerOfMass.y = (node->centerOfMass.y * node->totalMass + particle->position.y * particle->mass) / (node->totalMass + particle->mass);
        node->centerOfMass.z = (node->centerOfMass.z * node->totalMass + particle->position.z * particle->mass) / (node->totalMass + particle->mass);
    }
    node->totalMass += particle->mass;
}


void insertParticle(Node* node, Particle* particle) {
    // Empty Node
    if (node->particle == NULL && node->totalMass == 0) {
        node->particle = particle;
        updateMassCenterOfMass(node, particle);
        return;
    }

    //If the node is a leaf but already contains a particle, subdivide the node.
    if (node->particle != NULL && node->totalMass != 0) {
        subdivideNode(node);  
    }

    // Insert Child
    int newParticleIndex = OctantIndex(node, particle);
    if (node->children[newParticleIndex] == NULL) {
        node->children[newParticleIndex] = createChildNodeForOctant(node, newParticleIndex);
    }

    insertParticle(node->children[newParticleIndex], particle);

    // Update the parent node Center of Mass
    updateMassCenterOfMass(node, particle);
}


void subdivideNode(Node* node) 
{
    Vec3 center = {
        (node->bbox[0].x + node->bbox[1].x) / 2,
        (node->bbox[0].y + node->bbox[1].y) / 2,
        (node->bbox[0].z + node->bbox[1].z) / 2
    };

    // Create child nodes
    for (int ii = 0; ii < 8; ii++) {
        if (node->children[ii] == NULL) {
            node->children[ii] = createChildNodeForOctant(node, ii);
        }
    }
    
}


double distance(Vec3 ri, Vec3 rj)
{
    return sqrt( (rj.x - ri.x)*(rj.x - ri.x) + (rj.y-ri.y)*(rj.y-ri.y) + (rj.z-ri.z)*(rj.z-ri.z) );
}


void CenterOfMass(Node* node) {

    double totalMass = 0.0;
    Vec3 weightedP = (Vec3){0.0, 0.0, 0.0};

// If the node has a body directly, add its mass and its contribution to the center of mass
    if (node->particle != NULL) {
        totalMass += node->particle->mass;
        weightedP.x += node->particle->position.x * node->particle->mass;
        weightedP.y += node->particle->position.y * node->particle->mass;
        weightedP.z += node->particle->position.z * node->particle->mass;
    }

// Loop through the child nodes to accumulate their mass and contribution to the center of mass
    for (int i = 0; i < 8; ++i) {
        if (node->children[i] != NULL) {

            CenterOfMass(node->children[i]); // Llamada recursiva

            totalMass += node->children[i]->totalMass;
            weightedP.x += node->children[i]->centerOfMass.x * node->children[i]->totalMass;
            weightedP.y += node->children[i]->centerOfMass.y * node->children[i]->totalMass;
            weightedP.z += node->children[i]->centerOfMass.z * node->children[i]->totalMass;
        }
    }

    //If there is total mass, calculate the center of mass of the node
    if (totalMass > 0) {
        node->centerOfMass.x = weightedP.x / totalMass;
        node->centerOfMass.y = weightedP.y / totalMass;
        node->centerOfMass.z = weightedP.z / totalMass;
    }

    node->totalMass = totalMass; 
}



void Force(Node* node, Particle* particle, Vec3* f ){
    // Emty node
    if (node->totalMass == 0) return;
    
    Vec3 dir;
    double r_ij = distance(particle->position, node->centerOfMass);

    // Avoid divergence
    if (r_ij <= epsilon) return;

    // size of node
    double nodeSize = node->bbox[1].x - node->bbox[0].x;
    double h = nodeSize / r_ij;

    // If the node is far enough away or is a leaf, treat the node as a point mass
    if (h < theta || node->children[0] == NULL ){

        dir.x = node->centerOfMass.x - particle->position.x;
        dir.y = node->centerOfMass.y - particle->position.y;
        dir.z = node->centerOfMass.z - particle->position.z;

        double F = G * particle->mass * node->totalMass / (r_ij*r_ij*r_ij);
        f->x += F*dir.x;
        f->y += F*dir.y;
        f->z += F*dir.z;
    }
    else{
        // Otherwise calculate force recursively on subnodes
        for (size_t ii = 0; ii < 8; ii++){
            if (node->children[ii] != NULL){
                Force(node->children[ii], particle, f);
            }    
        }   
    }

}


void Euler(Particle* particle, Vec3 force, double dt) {
    // Calculate Aceleration
    Vec3 acceleration = {force.x / particle->mass, 
                         force.y / particle->mass, 
                         force.z / particle->mass};

    // Update Velocity
    particle->velocity.x += acceleration.x * dt;
    particle->velocity.y += acceleration.y * dt;
    particle->velocity.z += acceleration.z * dt;

    // Update Positions
    particle->position.x += particle->velocity.x * dt;
    particle->position.y += particle->velocity.y * dt;
    particle->position.z += particle->velocity.z * dt;
}


void freeNode(Node* node) {
    // Empty node
    if (node == NULL) {
        return;
    }
    
    // Run over nodes
    for (size_t ii = 0; ii < 8; ii++) {
        if (node->children[ii] != NULL) {
            freeNode(node->children[ii]);
            node->children[ii] = NULL; 
        }
    }

    free(node);
}


void writePositionsToCSV(Particle* particles, int numParticles, const char* filename) {
    FILE* file = fopen(filename, "w"); 
    if (file == NULL) {
        printf("Error to open file.\n");
        return;
    }

    fprintf(file, "X,Y,Z\n");

    for (int i = 0; i < numParticles; i++) {
        fprintf(file, "%f,%f,%f\n", particles[i].position.x, particles[i].position.y, particles[i].position.z);
    }

    fclose(file);
}



void writeNodeToFile(Node *node, FILE *file) {
    if (node == NULL) return;
    
    // write limits of nodes
    fprintf(file, "%f,%f,%f,%f,%f,%f\n",
            node->bbox[0].x, node->bbox[0].y, node->bbox[0].z,
            node->bbox[1].x, node->bbox[1].y, node->bbox[1].z);
    
    // Recursivity for childs
    for (int ii = 0; ii < 8; ii++) {
        writeNodeToFile(node->children[ii], file);
    }
}





void spher2cartes(Vec3 *vec, double r, double sclftr) {
    double X1 = RAND();
    double X2 = RAND();
    double aux1 = (1.-2.*X1)*r;
    double aux2 = sqrt(r*r - aux1*aux1)*sclftr;
    vec->x = aux2 * cos(TPI*X2);
    vec->y = aux2 * sin(TPI*X2);
    vec->z = aux1*sclftr;
}

double g(double x) {
    double x2 = x*x;
    return x2 * pow(1 - x2, 3.5);
}

double rndm(double min, double max) {
    return min + (max - min) * RAND();
}

void frm2com(Particle *particles, int N) {
    Vec3 com_pos = {0.0, 0.0, 0.0}; // Center of mass position
    Vec3 com_vel = {0.0, 0.0, 0.0}; // Center of mass velocity
    for (int i = 0; i < N; i++) {
        com_pos.x += particles[i].mass * particles[i].position.x;
        com_pos.y += particles[i].mass * particles[i].position.y;
        com_pos.z += particles[i].mass * particles[i].position.z;
        com_vel.x += particles[i].mass * particles[i].velocity.x;
        com_vel.y += particles[i].mass * particles[i].velocity.y;
        com_vel.z += particles[i].mass * particles[i].velocity.z;
    }

    for (int i = 0; i < N; i++) {
        particles[i].position.x -= com_pos.x;
        particles[i].position.y -= com_pos.y;
        particles[i].position.z -= com_pos.z;
        particles[i].velocity.x -= com_vel.x;
        particles[i].velocity.y -= com_vel.y;
        particles[i].velocity.z -= com_vel.z;
    }
}

void adjust_units(Particle *particles, int N) {
    double epot = 0.0;
    double ekin = 0.0;
    double dist = 0.0;

    for (int ii = 0; ii < N; ii++) {

        Vec3 vi = particles[ii].velocity;
        ekin += particles[ii].mass * (vi.x*vi.x + vi.y*vi.y + vi.z*vi.z) / 2;

        for (int jj = ii + 1; jj < N; jj++) {
            Vec3 r_ij = {particles[ii].position.x - particles[jj].position.x,
                    particles[ii].position.y - particles[jj].position.y,
                    particles[ii].position.z - particles[jj].position.z};

            dist = sqrt(r_ij.x*r_ij.x + r_ij.y*r_ij.y + r_ij.z*r_ij.z);
            epot -= particles[ii].mass * particles[jj].mass / dist;
        }
    }

    double alpha = -2 * epot;
    double beta = 0.5/sqrt(ekin);
    for (int ii = 0; ii < N; ii++) {
        particles[ii].position.x *= alpha;
        particles[ii].position.y *= alpha;
        particles[ii].position.z *= alpha;
        particles[ii].velocity.x *= beta;
        particles[ii].velocity.y *= beta;
        particles[ii].velocity.z *= beta;
    }
}

void plummer_dist(Particle *particles, int N, double sqrt2, double invsclftr, double sqrtscldrt) {

    double X1, X2;
    double cum_mass, r, Ve;
    double m = 1.0/N;
    double cum_mass_min = 0.0;
    double cum_mass_max = m;

    for (int ii = 0; ii < N; ii++) {
        particles[ii].mass = m;

        // Position
        cum_mass = rndm(cum_mass_min, cum_mass_max);
        cum_mass_min = cum_mass_max;
        cum_mass_max += m;

        r = 1.0 / sqrt(pow(cum_mass, -2.0/3.0) - 1.0);
        Vec3 pos;
        spher2cartes(&pos, r, invsclftr);
        particles[ii].position = pos;

        // Velocity
        do {
            X1 = RAND();
            X2 = RAND();
        } while (0.1*X2 >= g(X1));

        Ve = sqrt2 * pow(1.0 + r*r, -0.25) * X1;
        Vec3 vel;
        spher2cartes(&vel, Ve, sqrtscldrt);
        particles[ii].velocity = vel;
    }

    frm2com(particles, N);
    adjust_units(particles, N);

}