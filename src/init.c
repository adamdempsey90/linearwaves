#include "linearwaves.h"

void init_params(void) {


    params.n = 32768;
    params.nrhs =3;
    params.h = .0316;
    params.mu = 0;
    params.delta = 0;
    params.nuindx = .5;
    params.eta = 0;
    params.alpha = 1e-4;
    params.omf = 1;
    params.f = 0;
    params.sig0 = 1;
    params.rmin = .05;
    params.rmax = 30;
    params.iso = FALSE;
    params.pcorrect = FALSE;
    params.ieps = 0;


    planet.mp = 1;
    planet.eps = .01;
    planet.eps2 = planet.eps*planet.eps;
    planet.a = 1;
    planet.indirect = TRUE;

    return;
}
