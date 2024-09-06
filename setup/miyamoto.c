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

void miyamoto_disk(double *Pos, double *Vel, double *Mass, const int Nl, const int N, const double R, const double M) {
    double a = 1.0;
    double b = 0.1;
    double r, z, rho, phi, sigma_rad, sigma_tan, sigma_z;
    double vx, vy, vz, v, sigma;
    double X1;
    double dm = M/N;
    int ii;
    for (ii = 0; ii < Nl; ii++) {
        Mass[ii] = dm;
        // Position
        do {
            r   = R*rndm(0.005, 1.0);
            z   = 0.1*R*rndm(0.005, 1.0);
            X1  = RAND();
            rho = density_miy(r, z, M, a, b);
        } while (X1 < rho);

        phi = rndm(0, 2.0*M_PI);
        Pos[3 * ii] = r*cos(phi);
        Pos[3 * ii + 1] = r*sin(phi);
        Pos[3 * ii + 2] = z;

        r = sqrt(r*r+z*z);
        
        /*
          Warning!!!! 
          Velocity profile is not stable yet.
        */
        v       = rotational_velocity_miy(r, z, M, a , b);
        sigma   = density_miy(r, z, M, a , b);

        vx  = -v*sin(phi);
		vy  = v*cos(phi);
		vz  = 0;

        sigma_rad   = radial_velocity_dispersion(r, v, sigma);
        sigma_tan   = tangential_velocity_dispersion(r, v, sigma);
        sigma_z     = vertical_velocity_dispersion(r, v, sigma);

		double sigmaRadX = sigma_rad * sin( phi );
		double sigmaRadY = sigma_rad * cos( phi );

        vx = rand_normal(vx, sigmaRadX);
        vy = rand_normal(vy, sigmaRadY);

		double sigmaTanX = sigma_tan * sin( phi );
		double sigmaTanY = sigma_tan * cos( phi );

		vx = rand_normal(vx, sigmaTanX);
        vy = rand_normal(vy, sigmaTanY);

		vz = rand_normal(vz, sigma_z);

        Vel[3 * ii] = vx;
        Vel[3 * ii + 1] = vy;
        Vel[3 * ii + 2] = vz;
    }
    frm2com(Pos, Vel, Mass, Nl);
}
