#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "utils.h"
#include "plummer.h"
#include "Kuzmin.h"


#define G 1.0
#define M_c 1.0         // Mass central body for Kuzmin galaxy
#define Q 1.5           // Toomre Parameter for Kuzmin galaxy
#define TPI (2*M_PI)
#define RAND() ((double)rand()/(double)(RAND_MAX))


// Function to simulate the galaxy collision setup
void galaxy_collision(double *Pos, double *Vel, double *Mass, const int N,const double R, const double M, const double seed) {
    srand(seed);

    // Parameters for the collision setup
    const int N_kuzmin = N - 200;
    const int N_plummer = N - 800;
    const double R_kuzmin = R; // Scale radius for Kuzmin galaxy
    const double M_kuzmin = M; // Total mass for Kuzmin galaxy
    const double R_plummer = R*0.1; // Scale radius for Plummer galaxy
    const double M_plummer = M*0.07; // Total mass for Plummer galaxy
    const double d_x = 0.0;    // Initial separation in the x-direction
    const double d_z = 5.0;   // Initial separation in the z-direction
    const double d_y = 0.0;   // Initial separation in the y-direction
    const double v_comp = 1.3; // Initial velocity of the companion galaxy in the x-direction

    
    // Initialize Kuzmin disk
    kuzmin_disk(Pos, Vel, Mass, N_kuzmin, R_kuzmin, M_kuzmin, seed);

    // Initialize Plummer sphere
    plummer_dist(Pos + 3 * N_kuzmin, Vel + 3 * N_kuzmin, Mass + N_kuzmin, N_plummer, R_plummer, M_plummer, seed);

    // Adjust Plummer galaxy to initial collision conditions
    for (int ii = 0; ii < N_plummer; ii++) {
        Pos[3 * (N_kuzmin + ii)]     += d_x;  // shift in x
        Pos[3 * (N_kuzmin + ii) + 1] += d_y;  // shift in y
        Pos[3 * (N_kuzmin + ii) + 2] += d_z;  // shift in z
        Vel[3 * (N_kuzmin + ii) + 2]     -= v_comp; // initial velocity in x
    }

    // Adjust both galaxies to center of mass frame
    frm2com(Pos, Vel, Mass, N_kuzmin + N_plummer);
}