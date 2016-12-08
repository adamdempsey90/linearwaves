#include "linearwaves.h"
#include <gsl/gsl_integration.h>
typedef struct GSL_params {
    double x;
    int m;
} GSL_params;


void init_params(void) {

    params.mp = 1e-5;
    params.softening = .6*.05;
    params.softening2 = params.softening*params.softening;
    params.a = 1;
    params.indirect = TRUE;
    return;

}

double potential(double phi, double x) {
    double res = -pow(x*x + params.a*params.a + params.softening2 - 2*params.a*x*cos(phi),-.5);
    if (params.indirect) {
        res += cos(phi)*params.a/(x*x);
    }
    return res;
}

double dr_potential(double phi, double x) {
    double res = ( x - params.a*cos(phi)) * pow(x*x + params.a*params.a + params.softening2 - 2*params.a*x*cos(phi),-1.5);

    if (params.indirect) {
        res -= 2* cos(phi)*params.a/(x*x*x);
    }
    return res;
}

double dp_potential(double phi, double x) {
    double res = params.a*x*sin(phi)* pow(x*x + params.a*params.a + params.softening2 - 2*params.a*x*cos(phi),-1.5);

    if (params.indirect) {
        res +=  sin(phi)*params.a/(x*x);
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
    
    double tol = params.tol;
    double error;

    double res_r, dr_res_r;
    gsl_integration_qag(&F1,0,M_PI, tol, tol, nspace, gsl_order , w, &res_r, &error);
    gsl_integration_qag(&F2,0,M_PI, tol, tol, nspace, gsl_order , w, &dr_res_r, &error);
    res_r *= params.mp/M_PI;
    dr_res_r *= params.mp/M_PI;

    res[0] = dr_res_r;// *-1;
    res[1] = res_r ;//* -I * m/x;
    res[2] = 0;
    gsl_integration_workspace_free(w); 
    return;
}
