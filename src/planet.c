#include "linearwaves.h"
#include <gsl/gsl_integration.h>
#include <fftw3.h>
typedef struct GSL_params {
    double x;
    int m;
} GSL_params;

fftw_plan pr2r,pr2r1d;

void init_planet(void) {

    planet.mp = 1e-5;
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

double dp_potential(double phi, double x) {
    double res = planet.a*x*sin(phi)* pow(x*x + planet.a*planet.a + planet.eps2 - 2*planet.a*x*cos(phi),-1.5);

    if (planet.indirect) {
        res +=  sin(phi)*planet.a/(x*x);
    }
    return res;
}

double pfunc(double phi, void *params) {
    GSL_params p = *(GSL_params *)params;
    double x  = p.x;
    int m = p.m;
    return cos(m*phi) * potential(phi,x);
}
double dpfunc(double phi, void *params) {
    GSL_params p = *(GSL_params *)params;
    double x  = p.x;
    int m = p.m;
    return cos(m*phi) * dr_potential(phi,x);
}

void force(double x, int m, double complex *res) {
    int nspace = 1000;
    int gsl_order = 5;
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(nspace);
    
    gsl_function F1, F2;
    F1.function = &pfunc;
    F2.function = &dpfunc;
    GSL_params p;
    p.x = x;
    p.m = m;
    F1.params = &p;
    F2.params = &p;
    
    double tol = 1e-7;
    double error;

    double res_r, dr_res_r;
    gsl_integration_qag(&F1,0,M_PI, tol, tol, nspace, gsl_order , w, &res_r, &error);
    gsl_integration_qag(&F2,0,M_PI, tol, tol, nspace, gsl_order , w, &dr_res_r, &error);
    res_r *= planet.mp/M_PI;
    dr_res_r *= planet.mp/M_PI;

    res[0] = -dr_res_r ;
    res[1] = res_r ;//* -I * m/x;
    res[2] = 0;
    gsl_integration_workspace_free(w); 
    return;
}

void fft_init(double *in, double *out,int nphi) {
    int rank = 1; 
    int n[] = {nphi}; /* 1d transforms of length 10 */
    int howmany = params.n;
    int idist = nphi;
    int odist = nphi;
    int istride = 1;
    int ostride = 1; /* distance between two elements in 
                                                               the same column */
    int *inembed = n, *onembed = n;
    fftw_r2r_kind kind[] = {FFTW_REDFT00};

    pr2r = fftw_plan_many_r2r(rank, n, howmany,in,inembed, istride, idist, out, onembed,ostride,odist, kind, FFTW_ESTIMATE);

    pr2r1d = fftw_plan_r2r_1d(nphi, in, out,FFTW_REDFT00, FFTW_ESTIMATE);
}

void do_fft(double *in, double *out) {
    fftw_execute(pr2r1d);
    return;
}
void free_fft(double *in, double *out) {
    free(in);
    free(out);
    fftw_destroy_plan(pr2r);
    fftw_destroy_plan(pr2r1d);
}


int main(int argc, char *argv[]) {

    planet.mp = 1;
    planet.a = 1.0;
    planet.eps = .6*.05;
    planet.eps2 = planet.eps*planet.eps;
    planet.indirect = TRUE;

    if (planet.indirect) {
        printf("Using indirect potential\n");
    }

    int NM = 100;
    int NR = 10;
    double rmin = .1;
    double rmax = 1.5;
    params.n = NR;
    
    double complex res[3];
    int m;
    double x = .8;
    int j,indx;
    int nphi = 200;
    double phi;
    double *out = (double *)malloc(sizeof(double)*nphi);
    double *in = (double *)malloc(sizeof(double)*nphi);
    force(x,2,res);
    printf("%lg\n",creal(res[1]));
    for(j=0;j<nphi;j++) {
        phi = 2*M_PI/(nphi-1) *j;
        in[j] = potential(phi,x);
    }




    fft_init(in,out,nphi);
            
    do_fft(in,out);
    printf("%lg\n",creal(out[1]));


    free_fft(in,out);



    return 1;
}
