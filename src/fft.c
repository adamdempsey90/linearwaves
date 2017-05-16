#include "linearwaves.h"
#include <fftw3.h>


double potential(double phi, double x, double eps2, int indirect) {
    double res = -pow(x*x + 1 + eps2 - 2*x*cos(phi),-.5);
    if (indirect) {
        //res += cos(phi)*1./(x*x);
        res += cos(phi)*x;
    }
    return res;
}

double dr_potential(double phi, double x, double eps2, int indirect) {
    double res = ( x - cos(phi)) * pow(x*x + 1 + eps2 - 2*x*cos(phi),-1.5);

    if (indirect) {
        res += cos(phi);
    }
    return res;
}



void fft_potential(double *r, double *pot1, double *pot2, int num_modes, int nphi, int howmany, double eps2, int indirect) {
    int rank = 1; 
    int n[] = {nphi}; 
    int idist = nphi;
    int odist = nphi;
    int istride = 1;
    int ostride = 1; 
    int *inembed = n, *onembed = n;
    int i,j;
    fftw_r2r_kind kind[] = {FFTW_REDFT00};

    int norm = 2*(nphi-1);
    double *in = (double *)malloc(sizeof(double)*howmany*nphi);

    fftw_plan pr2r = fftw_plan_many_r2r(rank, n, howmany,in,inembed, istride, idist, in, onembed,ostride,odist, kind, FFTW_ESTIMATE);


    for(i=0;i<howmany;i++) {
        for(j=0;j<nphi;j++) {
            in[j + i*nphi] = potential(M_PI*j/(nphi-1),r[i],eps2,indirect);
        }
    }
    fftw_execute(pr2r);
    for(i=0;i<howmany;i++) {
        for(j=1;j<num_modes+1;j++) {
            pot1[i + (j-1)*howmany] = in[j + i*nphi]/norm;
        }
    }
    for(i=0;i<howmany;i++) {
        for(j=0;j<nphi;j++) {
            in[j + i*nphi] = dr_potential(M_PI*j/(nphi-1),r[i],eps2,indirect);
        }
    }
    fftw_execute(pr2r);
    for(i=0;i<howmany;i++) {
        for(j=1;j<num_modes+1;j++) {
            pot2[i + (j-1)*howmany] = in[j + i*nphi]/norm;
        }
    }
    fftw_destroy_plan(pr2r);
    SAFE_FREE(in);
    return;
}

