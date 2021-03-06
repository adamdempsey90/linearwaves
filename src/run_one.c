#include "linearwaves.h"
#include <time.h>


int run_one(int m, char *fname) {
    int i;
    clock_t tic, toc;
    double timers[3];
    const char *timer_names[3] = {"Matrix: ", "Solve:  ", "Torque: "};



    printf("Using m = %d\nSaving to %s\n",m,fname);

    init_planet();
    init_params();
    double *r = (double *)malloc(sizeof(double)*params.n);
    double complex *md =  (double complex *)malloc(sizeof(double complex)*params.n * params.nrhs * params.nrhs);
    double complex *fd =  (double complex *)malloc(sizeof(double complex)*params.n * params.nrhs);
    double complex *ud =  (double complex *)malloc(sizeof(double complex)*(params.n -1)* params.nrhs * params.nrhs);
    double complex *ld =  (double complex *)malloc(sizeof(double complex)*(params.n -1)* params.nrhs * params.nrhs);
    double *lamex = (double *)malloc(sizeof(double)*params.n);
    double *lamdep = (double *)malloc(sizeof(double)*params.n);
    double *fw= (double *)malloc(sizeof(double)*params.n);
    double *drfw= (double *)malloc(sizeof(double)*params.n);


    params.dlr = log(params.rmax/params.rmin) / (double)params.n;
    for(i=0;i<params.n;i++) {
        r[i] = params.rmin * exp(i*params.dlr);
    }

    tic = clock();
    construct_matrix(r,ld,md,ud,fd,m);
    toc = clock();
    timers[0] = (double)(toc - tic)/CLOCKS_PER_SEC;


    output_matrix(ld,md,ud,fd);

    tic = clock();
    cthomas_alg_block(ld,md,ud,fd,params.n,params.nrhs);
    toc = clock();
    timers[1] = (double)(toc - tic)/CLOCKS_PER_SEC;

    
    tic = clock();
    calc_torques(r,fw,drfw,lamex,lamdep,fd,m);
    toc = clock();
    timers[2] = (double)(toc - tic)/CLOCKS_PER_SEC;

    output(r,fd,lamex,lamdep,drfw,fw,fname);


    for(i=0;i<3;i++) printf("%s\t%.3e s\n",timer_names[i],timers[i]);


    SAFE_FREE(r); SAFE_FREE(md); SAFE_FREE(fd); SAFE_FREE(ud); SAFE_FREE(ld);
    SAFE_FREE(fw); SAFE_FREE(lamex);
    return 1;
}
