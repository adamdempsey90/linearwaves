#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>

#define TRUE 1
#define FALSE 0

typedef struct Params {
    
    int n,nrhs,nm,nphi;

    double h, mu, delta, nuindx, eta, alpha,omf,f,sig0;
    double dlr,rmin,rmax,tol;
    double ieps;
    int iso,pcorrect,zero_inner_bc, zero_outer_bc;
    int simple_visc;


}    Params;
typedef struct Planet {
    
    double mp, eps,eps2, a;

    int indirect;

}    Planet;
typedef struct Grid {
    int n, nm;
    int *mvals;
    double *r;
//    double complex *ld,*md,*ud,*fd;
    double *lamdep, *lamex, *fw, *drfw;
    double *TL, *TR;
    double complex *u, *v, *s;
    double  *dppot, *drpot;

} Grid;

typedef struct Disk {
    double *sigma;
    double *omega;
    double *drsigma;
    double *dr2sigma;
    double *dromega;
    double *kappa2;

} Disk;


Planet planet;
Params params;
double *r;
int np, rank;


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
void construct_matrix(double *r, double complex *ld, double complex *md, double complex *ud, double complex *fd, double *dppot, double *drpot, int m);



void cthomas_alg_block(double complex *a, double complex *b, double complex *c, double complex *d, int n, int m);
void cconstruct_total_matrix(double complex *ld, double complex *md, double complex *ud, double complex *mat, int n, int m);
void output(double *r, double complex *sol, double *lamex, double *lamdep, double *drfw, double *fw,char *fname);
void init_params(void);
void init_planet(void);
void output_matrix(double complex *ld, double complex *md, double complex *ud, double complex *fd);
void calc_torques(double *r, double *fw, double *drfw, double *lamex, double *lamdep, double complex *sol, double *dppot, double *TL, double *TR, int m) ;
void viscosity_coeffs_u(double r, double complex *res, int m);
void viscosity_dcoeffs_u(double r, double complex *res, int m, double invdlr);
void viscosity_d2coeffs_u(double r, double complex *res, int m, double invdlr2);
void viscosity_coeffs_v(double r, double complex *res, int m);
void viscosity_dcoeffs_v(double r, double complex *res, int m, double invdlr);
void viscosity_d2coeffs_v(double r, double complex *res, int m, double invdlr2);
void init_grid(int num_modes, Grid *grid) ;
void free_grid(Grid *grid) ;
void linearwaves(int i, Grid *grid) ;
void output_torques(char *fname, Grid *grid);
void fft_potential(double *r, double *pot1, double *pot2, int num_modes);
