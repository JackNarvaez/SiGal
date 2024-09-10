// --------------------------------------------------------------- //
// The exponential disk model.                                     //
// --------------------------------------------------------------- //

#include <math.h>
#include <gsl/gsl_sf_bessel.h>
#include "utils.h"  

double cum_mass_exp(double r) {
    return (1.-(r+1.)*exp(-r));
}

double rotational_velocity_exp(double r, double M) {
    double rm   = r/2;
    double fac  = gsl_sf_bessel_I0(rm)*gsl_sf_bessel_K0(rm) - gsl_sf_bessel_I1(rm)*gsl_sf_bessel_K1(rm);
    return rm*sqrt(2*M*fac);
}

double surface_density_exp(double r, double M) {
    return (M/(2*M_PI))*exp(-r);
}

double f_root_exp(double r, double x) {
    return x-cum_mass_exp(r);
}

void exponential_disk(double *Pos, double *Vel, double *Mass, int *i, const int Nl, const int N, const double R, const double M, const int I) {
    double r, phi, sigma_rad, sigma_tan, sigma_z;
    double vx, vy, vz, v, sigma;
    double X1;
    double dm = M/N;
    int ii, a;
    for (ii = 0; ii < Nl; ii++) {
        Mass[ii]= dm;
        i[ii]   = I;
        // Position
        do {
            a   = 1;
            X1  = RAND();
            r   = ZerosBisection(X1, &a, f_root_exp);
        } while (a==0 || r < 0.001*R || r > R);

        phi = rndm(0, 2*M_PI);
        Pos[3 * ii]     = r * cos(phi);
        Pos[3 * ii + 1] = r * sin(phi);
        Pos[3 * ii + 2] = 0.0;

        v       = rotational_velocity_exp(r, M);
        sigma   = surface_density_exp(r, M);

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