#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>

#define TRUE 1
#define FALSE 0

typedef struct Params {
    
    int n,nrhs;

    double h, mu, delta, nuindx, eta, alpha,omf,f,sig0;
    double dlr,rmin,rmax;
    double ieps;
    int iso;


}    Params;
typedef struct Planet {
    
    double mp, eps,eps2, a;

    int indirect;

}    Planet;


Planet planet;
Params params;
double *r;


double scaleH(double x) ;
double omegaK(double x) ;
double nu(double x) ;
double cs(double x) ;
double cs2(double x) ;
double omega2(double x) ;
double kappa2(double x) ;
double kappa(double x) ;
double omega(double x) ;
double k2om(double x);
double sigma(double x) ;
double pres(double x) ;
double dsdr(double x) ;
double dpdr(double x) ;
double dc2dr(double x) ;
double Dfunc(double x, double omp, int m) ;
double complex sig(double x, double complex s);

void force(double x, int m, double complex *res);
void construct_matrix(double *r, double complex *ld, double complex *md, double complex *ud, double complex *fd, int m);



void cthomas_alg_block(double complex *a, double complex *b, double complex *c, double complex *d, int n, int m);
void cconstruct_total_matrix(double complex *ld, double complex *md, double complex *ud, double complex *mat, int n, int m);
void output(double *r, double complex *sol, double *lamex, double *fw, char *fname);
void init_params(void);
void init_planet(void);
void output_matrix(double complex *ld, double complex *md, double complex *ud, double complex *fd);
