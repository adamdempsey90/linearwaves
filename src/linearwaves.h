#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>

#define TRUE 1
#define FALSE 0

#define SAFE_FREE(ptr) free(ptr); ptr = NULL;

typedef struct Params {
    
    int n,nrhs,nm,nphi;
    int mstart,mend;

    double h, mu, delta, nuindx, eta, alpha,omf,f,sig0;
    double dlr,rmin,rmax,tol;
    double ieps;
    int iso,pcorrect,zero_inner_bc, zero_outer_bc;
    int simple_visc;
    int indirect;
    int fromfile;
    char outputname[512];
    char diskfile[512];
    double eps;


}    Params;
typedef struct Planet {
    
    double mp, eps,eps2, a;

    int indirect;

}    Planet;
typedef struct Grid {
    int n, nm;
    int *mvals;
    double *r, *lr;
//    double complex *ld,*md,*ud,*fd;
    double *lamdep, *lamex, *fw, *drfw;
    double *TL, *TR;
    double complex *u, *v, *s;
    double  *dppot, *drpot;

} Grid;

typedef struct Disk {
    double *c2;
    double *sigma;
    double *omega;
    double *dlsdlr;
    double *d2lsdlr;
    double *dlomdlr;
    double *dlnudlr;
    double *dlTdlr;
    double *d2lTdlr;
    double *kappa2;
    double *nu;
    double *pres;
    double *dpdr;
    double *lambda;
} Disk;

int np, rank;


double scaleH(double x, Params params, Disk disk) ;
double omegaK(double x, Params params, Disk disk) ;
double nu_func(double x, Params params, Disk disk) ;
double cs_func(double x, Params params, Disk disk) ;
double cs2_func(double x, Params params, Disk disk) ;
double omega2_func(double x, Params params, Disk disk) ;
double kappa2_func(double x, Params params, Disk disk) ;
double kappa_func(double x, Params params, Disk disk) ;
double omega_func(double x, Params params, Disk disk) ;
double k2om_func(double x, Params params, Disk disk);
double sigma_func(double x, Params params, Disk disk) ;
double pres_func(double x, Params params, Disk disk) ;
double dsdr_func(double x, Params params, Disk disk) ;
double dpdr_func(double x, Params params, Disk disk) ;
double dc2dr_func(double x, Params params, Disk disk) ;
double Dfunc(int indx, double omp, int m, Params params, Disk disk) ;
double complex sig(int i, double complex s, Params params, Disk disk);

void force(double x, int m, double complex *res);
void construct_matrix(double *r, double complex *ld, double complex *md, double complex *ud, double complex *fd, double *dppot, double *drpot, int m,Params params, Disk disk);



void cthomas_alg_block(double complex *a, double complex *b, double complex *c, double complex *d, int n, int m);
void cconstruct_total_matrix(double complex *ld, double complex *md, double complex *ud, double complex *mat, int n, int m);
void output(double *r, double complex *sol, double *lamex, double *lamdep, double *drfw, double *fw,char *fname);
void init_params(void);
void init_planet(void);
void output_matrix(double complex *ld, double complex *md, double complex *ud, double complex *fd);
void calc_torques(double *r, double *fw, double *drfw, double *lamex, double *lamdep, double complex *sol, double *dppot, double *TL, double *TR, int m) ;
void viscosity_coeffs_u(int indx, double r, double complex *res, int m, double fac);
void viscosity_dcoeffs_u(int indx, double r, double complex *res, int m, double fac);
void viscosity_d2coeffs_u(int indx, double r, double complex *res, double fac);
void viscosity_coeffs_v(int indx, double r, double complex *res, int m, double fac);
void viscosity_dcoeffs_v(int indx, double r, double complex *res, int m, double fac);
void viscosity_d2coeffs_v(int indx, double r, double complex *res, double fac);
void init_grid(int mstar,int mend, Grid *grid) ;
void free_grid(Grid *grid) ;
void linearwaves(int i, Grid *grid) ;
void output_torques(char *fname, Grid *grid);
void fft_potential(double *r, double *pot1, double *pot2, int num_modes);
void read_sigma(char *fname, double *lr, double *sigma, double *dlsdlr, double *d2lsdlr, int n);
void init_disk(char *fname, double *r);
void free_disk(void);
void output_disk(char *fname);
void read_param_file(char *fname, Params params, int argc, char *argv[]);
