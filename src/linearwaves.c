#include "linearwaves.h"

void linearwaves(int i, Grid *grid) {
    
    int m = grid->mvals[i];
    double *r = grid->r;
    double complex *md =  (double complex *)malloc(sizeof(double complex)*params.n * params.nrhs * params.nrhs);
    double complex *fd =  (double complex *)malloc(sizeof(double complex)*params.n * params.nrhs);
    double complex *ud =  (double complex *)malloc(sizeof(double complex)*(params.n -1)* params.nrhs * params.nrhs);
    double complex *ld =  (double complex *)malloc(sizeof(double complex)*(params.n -1)* params.nrhs * params.nrhs);
    double *lamdep = &grid->lamdep[(i)*grid->n];
    double *lamex = &grid->lamex[(i)*grid->n];
    //double *ilamex = &grid->ilamex[(i)*grid->n];
    double *fw = &grid->fw[(i)*grid->n];
    double *drfw = &grid->drfw[(i)*grid->n];
    double *TL = &grid->TL[i];
    double *TR = &grid->TR[i];
    double complex *u = &grid->u[(i)*grid->n];
    double complex *v = &grid->v[(i)*grid->n];
    double complex *s = &grid->s[(i)*grid->n];
    double *dppot = &grid->dppot[(i)*grid->n];
    double *drpot = &grid->drpot[(i)*grid->n];
    FILE *f;
    int j;
    if (m==10) {
        f = fopen("potential.dat","w");
        for(j=0;j<params.n;j++) fprintf(f,"%lg\n",dppot[j]);
        fclose(f);
    }
    //printf("Working on m=%d\n",m);
    construct_matrix(r,ld,md,ud,fd,dppot,drpot,m);
    cthomas_alg_block(ld,md,ud,fd,params.n,params.nrhs);

    for(i=0;i<grid->n;i++) {
        u[i] = fd[i*params.nrhs];
        v[i] = fd[i*params.nrhs+1];
        s[i] = sig(i,fd[i*params.nrhs+2]);
    }

    calc_torques(r,fw,drfw,lamex,lamdep,fd,dppot,TL,TR,m);
    
    free(md);
    free(fd);
    free(ud);
    free(ld);
    return;
}


void init_grid(int num_modes, Grid *grid) {
    grid->mvals = (int *)malloc(sizeof(int)*num_modes);
    grid->r = (double *)malloc(sizeof(double)*params.n);
    grid->lr = (double *)malloc(sizeof(double)*params.n);
    grid->lamex = (double *)malloc(sizeof(double)*params.n*num_modes);
    grid->lamdep = (double *)malloc(sizeof(double)*params.n*num_modes);
    grid->fw= (double *)malloc(sizeof(double)*params.n*num_modes);
    grid->drfw= (double *)malloc(sizeof(double)*params.n*num_modes);
    grid->TL = (double *)malloc(sizeof(double)*num_modes);
    grid->TR = (double *)malloc(sizeof(double)*num_modes);
    grid->u = (double complex *)malloc(sizeof(double)*num_modes*params.n);
    grid->v = (double complex *)malloc(sizeof(double)*num_modes*params.n);
    grid->s = (double complex *)malloc(sizeof(double)*num_modes*params.n);
    grid->dppot = (double *)malloc(sizeof(double)*num_modes*params.n);
    grid->drpot = (double *)malloc(sizeof(double)*num_modes*params.n);

    return;
}
void free_grid(Grid *grid) {
    printf("1\n");
    free(grid->mvals);
    printf("2\n");
    free(grid->r );
    printf("3\n");
    free(grid->lamex );
    printf("4\n");
    free(grid->lamdep );
    printf("5\n");
    free(grid->fw);
    printf("6\n");
    free(grid->drfw);
    printf("7\n");
    free(grid->TL);
    printf("8\n");
    free(grid->TR);
    printf("9\n");
    free(grid->u);
    printf("10\n");
    free(grid->v);
    printf("11\n");
    free(grid->s);
    printf("12\n");
    free(grid->dppot);
    printf("13\n");
    free(grid->drpot);
    printf("14\n");
    return;
}
