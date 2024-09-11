#include <math.h>
#include "utils.h"
#include "plummer.h"
#include "kuzmin.h"
#include "isothermal.h"

void mwg2(double *Pos, double *Vel, double *Mass, int *i, const int Nl, const int N,const double R, const double M, const int I) {

    const int   N_D_T   = 3*N/7;
    const int   N_B_T   = N_D_T/3;
    const int   N_H_T   = N_D_T;

    // Disk
    const double M_D    = M;
    const double R_D    = R;
    const int    N_D    = 3*Nl/7;
    kuzmin_disk(Pos, Vel, Mass, i, N_D, N_D_T, R_D, M_D, I);

    // Bulge
    const double M_B    = M_D/9;
    const double R_B    = 0.125*R_D;
    const int    N_B    = N_D/3;
    plummer_dist(Pos+3*N_D, Vel+3*N_D, Mass+N_D, i+N_D, N_B, N_B_T, R_B, M_B, I+1);

    // Hallo
    const double M_H    = 20/0.9*M_D;
    const double R_H    = 250./12*R_D;
    const int    N_H    = N_D;
    isothermal_dist(Pos+3*(N_D+N_B), Vel+3*(N_D+N_B), Mass+N_D+N_B, i+N_D+N_B, N_H, N_H_T, R_H, M_H, I+2);

    // Adjust both galaxies to center of mass frame
    frm2com(Pos, Vel, Mass, Nl);
}