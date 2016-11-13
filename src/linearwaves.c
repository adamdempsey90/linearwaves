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


    construct_matrix(r,ld,md,ud,fd,m);
    cthomas_alg_block(ld,md,ud,fd,params.n,params.nrhs);
    calc_torques(r,fw,drfw,lamex,lamdep,fd,TL,TR,m);
    
    free(md);
    free(fd);
    free(ud);
    free(ld);
    return;
}


void init_grid(int num_modes, Grid *grid) {
    grid->mvals = (int *)malloc(sizeof(int)*num_modes);
    grid->r = (double *)malloc(sizeof(double)*params.n);
    grid->lamex = (double *)malloc(sizeof(double)*params.n*num_modes);
    grid->lamdep = (double *)malloc(sizeof(double)*params.n*num_modes);
    grid->fw= (double *)malloc(sizeof(double)*params.n*num_modes);
    grid->drfw= (double *)malloc(sizeof(double)*params.n*num_modes);
    grid->TL = (double *)malloc(sizeof(double)*num_modes);
    grid->TR = (double *)malloc(sizeof(double)*num_modes);

    return;
}
void free_grid(Grid *grid) {
    free(grid->mvals);
    free(grid->r );
    free(grid->lamex );
    free(grid->lamdep );
    free(grid->fw);
    free(grid->drfw);
    free(grid->TL);
    free(grid->TR);
    return;
}
