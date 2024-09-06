// --------------------------------------------------------------- //
// The Kuzmin’s disk model. Kuzmin (1956)                          //
// --------------------------------------------------------------- //

#include <math.h>
#include "utils.h"

double rotational_velocity_kuzmin(double r, double M) {
    return sqrt(M*r*r*pow(r*r+1., -1.5));
}

double surface_density_kuzmin(double r, double M) {
    return (M/(2*M_PI))*pow(1+r*r, -1.5);
}

void kuzmin_disk(double *Pos, double *Vel, double *Mass, const int Nl, const int N, const double R, const double M) {
    double r, phi, sigma_rad, sigma_tan, sigma_z;
    double vx, vy, vz, v, sigma;
    double X1;
    double dm = M/N;
    int ii;
    for (ii = 0; ii < Nl; ii++) {
        Mass[ii] = dm;
        // Position
        do {
            X1  = 1.- RAND();
            r   = sqrt(1./(X1*X1)-1.);
        } while (r > R || r < 0.005*R);

        phi = rndm(0, 2.0*M_PI);
        Pos[3 * ii] = r*cos(phi);
        Pos[3 * ii + 1] = r*sin(phi);
        Pos[3 * ii + 2] = 0.0;

        v       = rotational_velocity_kuzmin(r, M);
        sigma   = surface_density_kuzmin(r, M);

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