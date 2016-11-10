#include "linearwaves.h"
#include <gsl/gsl_integration.h>
typedef struct GSL_params {
    double x;
    int m;
} GSL_params;


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
    res[1] = res_r * -I * m/x;
    res[2] = 0;
    gsl_integration_workspace_free(w); 
    return;
}

/*
int main(int argc, char *argv[]) {

    planet.mp = 1e-5;
    planet.a = 1.0;
    planet.eps = .6*.05;
    planet.eps2 = planet.eps*planet.eps;
    planet.indirect = TRUE;

    if (planet.indirect) {
        printf("Using indirect potential\n");
    }

   
    int NM = atoi(argv[1]);
    int NR = atoi(argv[2]);


    double complex *res = (double complex *)malloc(sizeof(double complex)*NM*NR);
    double complex *dr_res = (double complex *)malloc(sizeof(double complex)*NM*NR);
    
    int m;
    double x;
    int j,indx;
    for(m=1;m<=NM;m++) {
        for(j=0;j<NR;j++) {
            x = .6 + j*(2.-.6)/NR;
            indx = j + (m-1)*NR;
            force(x,m,&res[indx],&dr_res[indx]);
        }
    }


    FILE *f = fopen("pot_test.dat","w");

    fwrite(res,sizeof(double complex),NM*NR,f);
    fwrite(dr_res,sizeof(double complex),NM*NR,f);
    fclose(f);



    return 1;
}
*/
