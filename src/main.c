#include "linearwaves.h"
#include <time.h>


int main(int argc, char *argv[]) {
    int i;
    clock_t tic, toc;
    double dt;
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
    double *lamex = (double *)malloc(sizeof(double)*params.n);
    double *fw= (double *)malloc(sizeof(double)*params.n);


    printf("Setting grid\n");
    params.dlr = log(params.rmax/params.rmin) / (double)params.n;
    for(i=0;i<params.n;i++) {
        r[i] = params.rmin * exp(i*params.dlr);
    }

    printf("Constructing matrix\n");
    tic = clock();
    construct_matrix(r,ld,md,ud,fd,m);
    toc = clock();
    dt = (double)(toc - tic)/CLOCKS_PER_SEC;
    printf("Construct Matrix took %lg\n seconds",dt);
    output_matrix(ld,md,ud,fd);
    printf("Solving\n");
    tic = clock();
    cthomas_alg_block(ld,md,ud,fd,params.n,params.nrhs);
    toc = clock();
    dt = (double)(toc - tic)/CLOCKS_PER_SEC;
    printf("Solve took %lg seconds\n",dt);

    double complex u,v,s;
    double complex res_f[3];
    for(i=0;i<params.n;i++) {
        u = fd[i*3];
        v = fd[i*3+1];
        s = sig(r[i],fd[i*3+2]);
        fw[i] = r[i]*r[i]*sigma(r[i])*2*creal(u * conj(v));
        force(r[i],m,res_f);
        lamex[i] = 2*M_PI*r[i]*2*creal(r[i]*conj(s)*res_f[1]);
    }

    printf("Outputting to %s\n",argv[2]);
    output(r,fd,lamex,fw,argv[2]);

    printf("Done\n");


    free(r); free(md); free(fd); free(ud); free(ld);
    free(fw); free(lamex);
    return 1;
}
