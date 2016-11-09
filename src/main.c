#include "linearwaves.h"


int main(int argc, char *argv[]) {
    int i;

    int m = atoi(argv[1]); 

    printf("Using m = %d\nSaving to %s\n",m,argv[2]);

    printf("Initializing\n");
    init_planet();
    init_params();
    printf("Allocating\n");
    double *r = (double *)malloc(sizeof(double)*params.n);
    double complex *md =  (double complex *)malloc(sizeof(double complex)*params.n * params.nrhs * params.nrhs);
    double complex *fd =  (double complex *)malloc(sizeof(double complex)*params.n * params.nrhs);
    double complex *ud =  (double complex *)malloc(sizeof(double complex)*(params.n -1)* params.nrhs * params.nrhs);
    double complex *ld =  (double complex *)malloc(sizeof(double complex)*(params.n -1)* params.nrhs * params.nrhs);


    printf("Setting grid\n");
    params.dlr = log(params.rmax/params.rmin) / (double)params.n;
    for(i=0;i<params.n;i++) {
        r[i] = params.rmin * exp(i*params.dlr);
    }


    printf("Constructing matrix\n");
    construct_matrix(r,ld,md,ud,fd,m);

    printf("Solving\n");
    cthomas_alg_block(ld,md,ud,fd,params.n,params.nrhs);

    printf("Outputting to %s\n",argv[2]);
    output(r,fd,argv[2]);

    printf("Done\n");


    free(r); free(md); free(fd); free(ud); free(ld);
    return 1;
}
