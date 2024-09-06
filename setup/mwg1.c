#include <math.h>
#include "utils.h"
#include "plummer.h"
#include "kuzmin.h"

void mwg1(double *Pos, double *Vel, double *Mass, const int Nl, const int N,const double R, const double M) {

    const int   N_D_T   = 4*N/7;
    const int   N_B_T   = N_D_T/4;
    const int   N_H_T   = N_D_T/2;

    // Disk
    const double M_D    = M;
    const double R_D    = R;
    const int    N_D    = 4*Nl/7;
    kuzmin_disk(Pos, Vel, Mass, N_D, N_D_T, R_D, M_D);

    // Bulge
    const double M_B    = 0.1*M_D;
    const double R_B    = 0.5*R_D;
    const int    N_B    = N_D/4;
    plummer_dist(Pos+3*N_D, Vel+3*N_D, Mass+N_D, N_B, N_B_T, R_B, M_B);

    // Hallo
    const double M_H    = 0.3*M_D;
    const double R_H    = R_D;
    const int    N_H    = N_D/2;
    plummer_dist(Pos+3*(N_D+N_B), Vel+3*(N_D+N_B), Mass+N_D+N_B, N_H, N_H_T, R_H, M_H);

    // Adjust both galaxies to center of mass frame
    frm2com(Pos, Vel, Mass, Nl);
}