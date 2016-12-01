#include "linearwaves.h"

void excited_torques(Grid *grid, double *TL, double *TR, double *sol, int m) {

    double lamex;
    int i;
    double r,s,dp, rfac;

    for(i=0;i<params.n;i++) {
        r = grid->r[i];
        rfac = r*r* 2*M_PI*params.dlr;
        s = sig(i,sol[i*params.nrhs+2]);
        dp = grid->dppot[i];    
        lamex = -2*m*creal(cimag(s)*dp);
        if (r >= planet.a) *TR += lamex * rfac;
        if (r <= planet.a) *TL -= lamex * rfac;
    }
    return;
}
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
    
    SAFE_FREE(md);
    SAFE_FREE(fd);
    SAFE_FREE(ud);
    SAFE_FREE(ld);
    return;
}


void init_grid(int mstart, int mend, Grid *grid) {
    int num_modes = (mend-mstart+1)/np;
    
    

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
    grid->n = params.n;
    params.dlr = log(params.rmax/params.rmin) / (double)params.n;
    int i;
    for(i=0;i<params.n;i++) {
        grid->lr[i] = log(params.rmin) + params.dlr*i;
        grid->r[i] = exp(grid->lr[i]);
    }
    init_disk(params.diskfile,grid->lr);

    grid->n = params.n;
    grid->nm = num_modes;
    for(i=0;i<num_modes;i++) {
        grid->mvals[i] = rank*num_modes + mstart + i;
    }
    double *dp_pot_full = (double *)malloc(sizeof(double)*params.n*(params.nm));
    double *dr_pot_full = (double *)malloc(sizeof(double)*params.n*(params.nm));
/* Do FFT of potential on root */
    if (rank == 0) {
        fft_potential(grid->r,dp_pot_full,dr_pot_full,params.nm);
    }
#ifdef _MPI
    MPI_Scatter(dp_pot_full,grid->n*num_modes,MPI_DOUBLE,grid->dppot,grid->n*num_modes,MPI_DOUBLE,0, MPI_COMM_WORLD);

    MPI_Scatter(dr_pot_full,grid->n*num_modes,MPI_DOUBLE,grid->drpot,grid->n*num_modes,MPI_DOUBLE,0, MPI_COMM_WORLD);
#else
    memcpy(grid->dppot,dp_pot_full,sizeof(double)*grid->n*num_modes);
    memcpy(grid->drpot,dr_pot_full,sizeof(double)*grid->n*num_modes);
#endif
    SAFE_FREE(dp_pot_full);
    SAFE_FREE(dr_pot_full);

    return;
}
void free_grid(Grid *grid) {
    printf("1\n");
    SAFE_FREE(grid->mvals);
    printf("2\n");
    SAFE_FREE(grid->r );
    printf("3\n");
    SAFE_FREE(grid->lamex );
    printf("4\n");
    SAFE_FREE(grid->lamdep );
    printf("5\n");
    SAFE_FREE(grid->fw);
    printf("6\n");
    SAFE_FREE(grid->drfw);
    printf("7\n");
    SAFE_FREE(grid->TL);
    printf("8\n");
    SAFE_FREE(grid->TR);
    printf("9\n");
    SAFE_FREE(grid->u);
    printf("10\n");
    SAFE_FREE(grid->v);
    printf("11\n");
    SAFE_FREE(grid->s);
    printf("12\n");
    SAFE_FREE(grid->dppot);
    printf("13\n");
    SAFE_FREE(grid->drpot);
    printf("14\n");
    return;
}


void get_excited_torques(int mstart, int mend, double *TL, double *TR) {
    Grid *grid = (Grid *)malloc(sizeof(Grid));
    init_grid(mstart,mend, grid);
    int i;
    int num_modes = mend-mstart+1;
    for(i=0;i<num_modes;i++) {
        linearwaves(i, grid);
        excited_torques(grid,TL,TR,grid->mvals[i]);
    }

    free_grid(grid);

    return;

}
