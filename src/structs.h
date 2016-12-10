typedef struct Proc {
    int np,rank;
} Proc;
typedef struct Params {
    int np,rank; 
    int n,nrhs,nm,nphi;
    int mstart,mend;

    double h, mu, delta, nuindx, eta, alpha,omf,f,sig0;
    double dlr,rmin,rmax,tol;
    double ieps;
    int iso,pcorrect,zero_inner_bc, zero_outer_bc;
    int simple_visc;
    int fromfile;
    char outputname[512];
    char diskfile[512];
    double eps;
    double eps2;
    int indirect;
    double softening, softening2;
    double a;


}    Params;
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
    double *vrbar;
    double *dlvrbar;
} Disk;

Proc proc;
