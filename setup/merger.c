#include <math.h>
#include <iostream>
#include "utils.h"
#include "plummer.h"
#include "kuzmin.h"


void galaxy_collision(double *Pos, double *Vel, double *Mass, int *i, const int Nl, const int N, const double R, const double M, const int I) {
    // N = 1200 
    // Parameters for each galaxy (Total N, mass and R)
    const int N_1  = N-4000;  // 700
    const int N_2  = N-20000;  // 500
    const int Nl_1 = (Nl * N_1) / (N_1 + N_2);
    const int Nl_2 = Nl - Nl_1;                

    const double R_1 = R * 12.0;
    const double R_2 = R * 2.0;
    const double M_1 = 10.;
    const double M_2 = .6;

    // GALAXY 1
    const int   N_D_T_1   = N_1* (3./7.);
    const int   N_B_T_1   = N_1* (1./7.); 
    const int   N_H_T_1   = N_1* (3./7.);
    // Disk
    const double M_D_1    = M_1;
    const double R_D_1    = R_1;
    //const int N_D_1 = Nl_1 * (3.0 / 7.0);
    const int N_D_1 = static_cast<int>(std::round(Nl_1 * (3.0 / 7.0)));


    kuzmin_disk(Pos, Vel, Mass, i, N_D_1, N_D_T_1, R_D_1, M_D_1, I);
    // Bulge
    const double M_B_1    = 0.1 * M_D_1;
    const double R_B_1    = 0.5 * R_D_1;
    //const int N_B_1 = Nl_1 * (1.0 / 7.0);
    const int N_B_1 = static_cast<int>(std::round(Nl_1 * (1.0 / 7.0)));


    plummer_dist(Pos + 3 * N_D_1, Vel + 3 * N_D_1, Mass + N_D_1, i + N_D_1, N_B_1, N_B_T_1, R_B_1, M_B_1, I + 1);
    // Dark Halo
    const double M_H_1    = 3. * M_D_1;
    const double R_H_1    = 3. * R_D_1;
    const int N_H_1 = Nl_1 - N_D_1 - N_B_1;  // Asigna lo que falta al halo

    plummer_dist(Pos + 3 * (N_D_1 + N_B_1), Vel + 3 * (N_D_1 + N_B_1), Mass + N_D_1 + N_B_1, i + N_D_1 + N_B_1, N_H_1, N_H_T_1, R_H_1, M_H_1, I + 2);

    frm2com(Pos, Vel, Mass, Nl_1);

    // GALAXY 2
    const int   N_D_T_2   = N_2 * (2./5.);
    const int   N_B_T_2   = N_2 * (1./5.);
    const int   N_H_T_2   = N_2 * (2./5.);
    // Disk
    const double M_D_2    = M_2;
    const double R_D_2    = R_2;
    //const int N_D_2 = Nl_2 * (2.0 / 5.0);
    const int N_D_2 = static_cast<int>(std::round(Nl_2 * (2.0 / 5.0)));


    kuzmin_disk(Pos + 3 * (Nl_1), Vel + 3 * (Nl_1), Mass + Nl_1, i + Nl_1, N_D_2, N_D_T_2, R_D_2, M_D_2, I + 3);
    // Bulge
    const double M_B_2    = 0.1 * M_D_2;
    const double R_B_2    = 0.5 * R_D_2;
    //const int N_B_2 = Nl_2 * (1.0 / 5.0);
    const int N_B_2 = static_cast<int>(std::round(Nl_2 * (1.0 / 5.0)));

    plummer_dist(Pos + 3 * (Nl_1 + N_D_2), Vel + 3 * (Nl_1 + N_D_2), Mass + Nl_1 + N_D_2, i + Nl_1 + N_D_2, N_B_2, N_B_T_2, R_B_2, M_B_2, I + 4);   
    // Dark Halo
    const double M_H_2    = 3 * M_D_2;
    const double R_H_2    = 3*R_D_2;
    const int N_H_2 = Nl_2 - N_B_2 - N_D_2;
    plummer_dist(Pos + 3 * (Nl_1 + N_D_2 + N_B_2), Vel + 3 * (Nl_1 + N_D_2 + N_B_2), Mass + Nl_1 + N_D_2 + N_B_2, i + Nl_1 + N_D_2 + N_B_2, N_H_2, N_H_T_2, R_H_2, M_H_2, I + 5);
    
    frm2com(Pos + 3 * Nl_1, Vel + 3 * Nl_1, Mass + Nl_1, Nl_2);

    // Separation of galaxy
    const double d_x = 12.0;   
    const double d_y = 10.0; 
    const double d_z = 0.0;   

    // Compute the direction from Plummer galaxy to Kuzmin galaxy center
    double dir_x = -d_x;  
    double dir_y = -d_y;
    double dir_z = -d_z;
    double norm = sqrt(dir_x * dir_x + dir_y * dir_y + dir_z * dir_z);

    // Normalize the velocity vector towards the center of galaxy
    double v_comp_x = dir_x / norm;
    double v_comp_y = dir_y / norm;
    double v_comp_z = dir_z / norm;

    for (int ii = 0; ii < Nl_2; ii++) { 
        Pos[3 * (Nl_1 + ii)]     += d_x; 
        Pos[3 * (Nl_1 + ii) + 1] += d_y; 
        Pos[3 * (Nl_1 + ii) + 2] += d_z; 
        Vel[3 * (Nl_1 + ii)]     += 0.6 * v_comp_x;  
        Vel[3 * (Nl_1 + ii) + 1] += 1.2 * v_comp_y; 
        Vel[3 * (Nl_1 + ii) + 2] += 0.0 * v_comp_z;  
    }
}