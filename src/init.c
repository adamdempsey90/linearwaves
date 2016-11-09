#include "linearwaves.h"

void init_params(void) {


    params.n = 512;
    params.nrhs = 1;
    params.h = .05;
    params.mu = -1;
    params.delta = 0;
    params.nuindx = .5;
    params.eta = 0;
    params.alpha = 0;
    params.omf = 1;
    params.f = 0;
    params.sig0 = 1;
    params.rmin = .1;
    params.rmax = 5.;
    params.iso = TRUE;


    planet.mp = 1e-5;
    planet.eps = .6*.05;
    planet.eps2 = planet.eps*planet.eps;
    planet.a = 1;
    planet.indirect = TRUE;

    return;
}
