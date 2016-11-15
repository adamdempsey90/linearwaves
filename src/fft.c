#include "linearwaves.h"
#include <fftw3.h>

void init_planet(void) {

    planet.mp = 1;
    planet.eps = .6*.05;
    planet.eps2 = planet.eps*planet.eps;
    planet.a = 1;
    planet.indirect = TRUE;
    return;

}

double potential(double phi, double x) {
    double res = -pow(x*x + planet.a*planet.a + planet.eps2 - 2*planet.a*x*cos(phi),-.5);
    if (planet.indirect) {
        res += cos(phi)*planet.a/(x*x);
    }
    return res;
}

double dr_potential(double phi, double x) {
    double res = ( x - planet.a*cos(phi)) * pow(x*x + planet.a*planet.a + planet.eps2 - 2*planet.a*x*cos(phi),-1.5);

    if (planet.indirect) {
        res -= 2* cos(phi)*planet.a/(x*x*x);
    }
    return res;
}



void fft_potential(double *r, double *pot1, double *pot2, int num_modes) {
    int nphi = params.nphi;
    int rank = 1; 
    int n[] = {nphi}; 
    int howmany = params.n;
    int idist = nphi;
    int odist = nphi;
    int istride = 1;
    int ostride = 1; 
    int *inembed = n, *onembed = n;
    int i,j;
    fftw_r2r_kind kind[] = {FFTW_REDFT00};

    int norm = 2*(nphi-1);
    double *in = (double *)malloc(sizeof(double)*params.n*nphi);

    fftw_plan pr2r = fftw_plan_many_r2r(rank, n, howmany,in,inembed, istride, idist, in, onembed,ostride,odist, kind, FFTW_ESTIMATE);


    for(i=0;i<params.n;i++) {
        for(j=0;j<nphi;j++) {
            in[j + i*nphi] = potential(M_PI*j/(nphi-1),r[i]);
        }
    }
    fftw_execute(pr2r);
    for(i=0;i<params.n;i++) {
        for(j=1;j<num_modes+1;j++) {
            pot1[i + (j-1)*params.n] = in[j + i*nphi]/norm;
        }
    }
    for(i=0;i<params.n;i++) {
        for(j=0;j<nphi;j++) {
            in[j + i*nphi] = dr_potential(M_PI*j/(nphi-1),r[i]);
        }
    }
    fftw_execute(pr2r);
    for(i=0;i<params.n;i++) {
        for(j=1;j<num_modes+1;j++) {
            pot2[i + (j-1)*params.n] = in[j + i*nphi]/norm;
        }
    }
    fftw_destroy_plan(pr2r);
    free(in);
    return;
}

