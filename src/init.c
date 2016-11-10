#include "linearwaves.h"

void init_params(void) {


    params.n = 8192;
    params.nrhs =3;
    params.h = .0316;
    params.mu = 0;
    params.delta = 0;
    params.nuindx = .5;
    params.eta = 0;
    params.alpha = 0.002;
    params.omf = 1;
    params.f = 0;
    params.sig0 = 1;
    params.rmin = .5;
    params.rmax = 1.8;
    params.iso = FALSE;
    params.ieps = 1e-2;


    planet.mp = 1;
    planet.eps = 1e-3;
    planet.eps2 = planet.eps*planet.eps;
    planet.a = 1;
    planet.indirect = TRUE;

    return;
}
