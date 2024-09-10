// --------------------------------------------------------------- //
// Isothermal sphere model.                                        //
// --------------------------------------------------------------- //
#include <math.h>
#include "utils.h"

void isothermal_dist(double *Pos, double *Vel, double *Mass, int *i, const int Nl, const int N, const double R, const double M, const int I) {
    double r, Ve;
    double dm = M/N;
    int ii;
    for (ii=0; ii<Nl; ii++) {
        Mass[ii]= dm;
        i[ii]   = I;

        // Position
        r  = R*rndm(0.001, 1.0);
        spher2cartes(Pos + 3*ii, r);

        Ve = sqrt(M/R);
        spher2cartes(Vel + 3*ii, Ve);
    }
    frm2com(Pos, Vel, Mass, Nl);
}