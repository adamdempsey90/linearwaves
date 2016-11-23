#include "linearwaves.h"

void init_params(void) {


    // Ints
    params.n = 8192;
    params.nphi = 2000;
    params.nrhs =3;

    // doubles
    params.h = .05;
    params.mu = -.5;
    params.delta = -1;
    params.nuindx = 0.5;
    params.eta = 2./3;
    params.alpha = 0.001;
    params.omf = 1;
    params.f = 0;
    params.sig0 = 1;
    params.rmin = .05;
    params.rmax = 30;
    params.ieps = 0;

    //Bools
    params.iso = TRUE;
    params.pcorrect = FALSE;
    params.zero_inner_bc = FALSE;
    params.zero_outer_bc = FALSE;
    params.simple_visc = FALSE;
    params.fromfile = FALSE;


    planet.mp = 1;
    planet.eps = .03;
    planet.eps2 = planet.eps*planet.eps;
    planet.a = 1;
    planet.indirect = TRUE;

    return;
}
