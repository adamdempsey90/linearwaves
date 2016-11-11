#include "linearwaves.h"

void init_params(void) {


    params.n = 32000;
    params.nrhs =3;
    params.h = .05;
    params.mu = -.5;
    params.delta = -1;
    params.nuindx = 0;
    params.eta = 2./3;
    params.alpha = 1e-3;
    params.omf = 1;
    params.f = 0;
    params.sig0 = 1;
    params.rmin = .01;
    params.rmax = 30;
    params.iso = TRUE;
    params.pcorrect = FALSE;
    params.ieps = 0;
    params.zero_inner_bc = FALSE;
    params.zero_outer_bc = FALSE;


    planet.mp = 1;
    planet.eps = .6*.05;
    planet.eps2 = planet.eps*planet.eps;
    planet.a = 1;
    planet.indirect = TRUE;

    return;
}
