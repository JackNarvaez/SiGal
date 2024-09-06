// --------------------------------------------------------------- //
// The Kepler’s disk with a massive central object.                //
// --------------------------------------------------------------- //

#include <math.h>
#include "utils.h"
   
double keplerian_velocity(double r, double M) {
    return sqrt(M / r);
}

void keplerian_disk(double *Pos, double *Vel, double *Mass, const int Nl, const int N, const double R, const double M) {
    double r, phi, th, z, vtan;
    double Md = 0.1*M;
    double dm = Md / N;
    int ii;
    for (ii = 0; ii < Nl; ii++) {
        Mass[ii]    = dm;
        r   = R*sqrt(rndm(0.005, 1.0));
        do {
            z   = rand_normal(0, 1);
        } while (fabs(z)>0.01*R);
        phi = rndm(0, 2*M_PI);
        Pos[3*ii]   = r * cos(phi);
        Pos[3*ii+1] = r * sin(phi);
        Pos[3*ii+2] = z;
        r   = sqrt(r*r + z*z);
        th  = atan2(r, z);
        vtan= keplerian_velocity(r, M);
        Vel[3*ii]   =   -vtan*sin(th)*sin(phi);
        Vel[3*ii+1] =   vtan*sin(th)*cos(phi);
        Vel[3*ii+2] =   vtan*cos(th);
    }
    frm2com(Pos, Vel, Mass, Nl);
}