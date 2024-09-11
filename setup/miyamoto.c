// --------------------------------------------------------------- //
// The Miyamoto and Nagai model for a flat disk.                   //
// Miyamoto & Nagai (1975)                                         //
// --------------------------------------------------------------- //

#include <math.h>
#include "utils.h"

double rotational_velocity_miy(double r, double z, double M, double a, double b) {
    double f1 = sqrt(z*z+b*b);
    return M*r*(f1+a*z*z/r*r)/(f1*pow(a*(2*f1+a)+r*r+b*b, 1.5));
}

double density_miy(double R, double z, double M, double a, double b) {
    double f1 = z*z + b*b;
    double f2 = sqrt(f1);
    double f3 = (a+f2)*(a+f2);
    double f4 = a*R*R + (a+3*f2)*f3/(pow(R*R+f3, 2.5)*pow(f1, 1.5));
    return (b*b*M/(4*M_PI))*f4;
}

void miyamoto_disk(double *Pos, double *Vel, double *Mass, int *i, const int Nl, const int N, const double R, const double M, const int I) {
    double x, y, z, r, rho, phi, sigma_rad, sigma_z;
    double vR, vx, vy, vz, S, Om;
    double dm   = M/N;
    double a    = 1.0;
    double b    = 0.05;
    int ii;
    double rho_max  = density_miy(a, 0, M, a, b);
    double rho_0    = density_miy(0, 0, M, a, b);
    for (ii = 0; ii < Nl; ii++) {
        Mass[ii] = dm;
        i[ii]    = I;
        // Position
        do {
            x   = rndm(-a, a);
            y   = rndm(-a, a);
            z   = 0.0; 
            r  = sqrt(x*x+y*y);   
            rho = rndm(rho_0, rho_max);
        } while (rho > density_miy(r, z, M, a, b) || r < 0.001*R);

        Pos[3 * ii] = x*R;
        Pos[3 * ii + 1] = y*R;
        Pos[3 * ii + 2] = z;

        vR = 0.1;              
        sigma_rad = vR/10;
        sigma_z = sqrt(b*b*density_miy(a, 0, M, a, b));

        vz  = rand_normal(0, sigma_z);
        vR  = rand_normal(0, sigma_rad);
    
        phi = atan2(y, x);
        vx  = vR*cos(phi);
        vy  = vR*sin(phi);
        
        S   = sqrt(r*r+pow(a+sqrt(z*z+b*b),2));
        Om  = pow(S,-1.5);

        vx += y*Om;
        vy += -x*Om;

        Vel[3 * ii] = vx;
        Vel[3 * ii + 1] = vy;
        Vel[3 * ii + 2] = vz;
    }
    frm2com(Pos, Vel, Mass, Nl);
}
