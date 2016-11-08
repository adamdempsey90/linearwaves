#include "linearwaves.h"


int main(int argc, char *argv[]) {
    int i;

    int m = atoi(argv[1]); 

    double *r = (double *)malloc(sizeof(double)*params.n);
    double complex *md =  (double complex *)malloc(sizeof(double complex)*params.n * params.nrhs * params.nrhs);
    double complex *fd =  (double complex *)malloc(sizeof(double complex)*params.n * params.nrhs);
    double complex *ud =  (double complex *)malloc(sizeof(double complex)*(params.n -1)* params.nrhs * params.nrhs);
    double complex *ld =  (double complex *)malloc(sizeof(double complex)*(params.n -1)* params.nrhs * params.nrhs);
    double complex *sol = (double complex *)malloc(sizeof(double complex)*params.n * params.nrhs);


    params.dlr = log(params.rmax/params.rmin) / (double)params.n;
    for(i=0;i<params.n;i++) {
        r[i] = params.rmin * exp(i*params.dlr);
    }

    construct_matrix(ld,md,ud,fd,m);
    thomas_alg_block(ld,md,ud,fd,sol,params.n,params.nrhs);
    output(sol);




    free(r); free(md); free(fd); free(ud); free(ld); free(sol);
    return 1;
}
