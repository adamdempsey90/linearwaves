#include "linearwaves.h"
/*
void excited_torques(double *r, double complex *sp, double *dppot, double *TL, double *TR, int m, Params params, Planet planet, Disk *disk) {

    double lamex;
    int i;
    double rval,rfac;

    for(i=0;i<params.n;i++) {
        rval = r[i];
        rfac = rval*rval* 2*M_PI*params.dlr;
        lamex = -2*m*creal(cimag(sp[i])*dppot[i]);
        if (rval >= planet.a) *TR += lamex * rfac;
        if (rval <= planet.a) *TL -= lamex * rfac;
    }
    return;
}
*/

void linearwaves(int i, Grid *grid, Params params, Disk *disk, int silent,int second_order) {
    
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
    double complex *u;
    double complex *v;
    double complex *s;
    if (second_order) {
        u = &grid->u2[(i)*grid->n];
        v = &grid->v2[(i)*grid->n];
        s = &grid->s2[(i)*grid->n];
    }
    else {
        u = &grid->u[(i)*grid->n];
        v = &grid->v[(i)*grid->n];
        s = &grid->s[(i)*grid->n];
    }
    double *dppot = &grid->dppot[(i)*grid->n];
    double *drpot = &grid->drpot[(i)*grid->n];
    double complex *Ru = &grid->Ru[i*(grid->n)];
    double complex *Rv = &grid->Rv[i*(grid->n)];
    double complex *Rs = &grid->Rs[i*(grid->n)];
    FILE *f;
    int j;
    /*
   if (m==2) {
       f = fopen("potential.dat","w");
       for(j=0;j<params.n;j++) fprintf(f,"%lg\n",dppot[j]);
       fclose(f);
   }
   */
    if (second_order) {
        construct_matrix_second(r,ld,md,ud,fd,Ru,Rv,Rs,m,params,disk);
    }
    else {
        construct_matrix(r,ld,md,ud,fd,drpot,dppot,m,params,disk);
    }
    /*
    if (m==2) {
    f = fopen("matrix.dat","w");
    for (j=0;j<params.n;j++) {
        fprintf(f,"%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n",r[j],
                creal(md[3*j]),cimag(md[3*j]),
                creal(md[3*j + 1]),cimag(md[3*j + 1]),
                creal(md[3*j + 2]),cimag(md[3*j + 2]));
    }
    fclose(f);
    }
    */
    cthomas_alg_block(ld,md,ud,fd,params.n,params.nrhs);

    for(i=0;i<grid->n;i++) {
        u[i] = fd[i*params.nrhs];
        v[i] = fd[i*params.nrhs+1];
        s[i] = sig(i,fd[i*params.nrhs+2],params,disk);
    }

    if (!second_order) {
        calc_torques(r,u,v,s,fw,drfw,lamex,lamdep,fd,dppot,TL,TR,m,params,disk,silent);
    }
    
    SAFE_FREE(md);
    SAFE_FREE(fd);
    SAFE_FREE(ud);
    SAFE_FREE(ld);
    return;
}


void init_grid(int mstart, int mend, Grid *grid, Params params, Disk *disk) {
    Proc localproc;
#ifndef _MPI
    localproc.rank = 0;
    localproc.np = 1;
#else
    localproc.rank = proc.rank;
    localproc.np = proc.np;
#endif
    int num_modes = (mend-mstart+1)/localproc.np;

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
    grid->u2 = (double complex *)malloc(sizeof(double)*num_modes*params.n);
    grid->v2 = (double complex *)malloc(sizeof(double)*num_modes*params.n);
    grid->s2 = (double complex *)malloc(sizeof(double)*num_modes*params.n);
    grid->Ru = (double complex *)malloc(sizeof(double)*num_modes*params.n);
    grid->Rv = (double complex *)malloc(sizeof(double)*num_modes*params.n);
    grid->Rs = (double complex *)malloc(sizeof(double)*num_modes*params.n);
    grid->dppot = (double *)malloc(sizeof(double)*num_modes*params.n);
    grid->drpot = (double *)malloc(sizeof(double)*num_modes*params.n);
    grid->n = params.n;
    int i;
    for(i=0;i<params.n;i++) {
        grid->lr[i] = log(params.rmin) + params.dlr*i;
        grid->r[i] = exp(grid->lr[i]);
    }
    init_disk(params.diskfile,grid->lr,disk,params);

    grid->n = params.n;
    grid->nm = num_modes;
    for(i=0;i<num_modes;i++) {
        grid->mvals[i] = localproc.rank*num_modes + mstart + i;
    }
    double *dp_pot_full = (double *)malloc(sizeof(double)*params.n*(params.nm));
    double *dr_pot_full = (double *)malloc(sizeof(double)*params.n*(params.nm));
/* Do FFT of potential on root */
    if (localproc.rank == 0) {
        fft_potential(grid->r,dp_pot_full,dr_pot_full,params.nm,params.nphi,params.n,params.eps2,params.indirect);
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
void free_linear_grid(Grid *grid) {
    SAFE_FREE(grid->mvals);
    SAFE_FREE(grid->r );
    SAFE_FREE(grid->lamex );
    SAFE_FREE(grid->lamdep );
    SAFE_FREE(grid->fw);
    SAFE_FREE(grid->drfw);
    SAFE_FREE(grid->TL);
    SAFE_FREE(grid->TR);
    SAFE_FREE(grid->u);
    SAFE_FREE(grid->v);
    SAFE_FREE(grid->s);
    SAFE_FREE(grid->u2);
    SAFE_FREE(grid->v2);
    SAFE_FREE(grid->s2);
    SAFE_FREE(grid->Ru);
    SAFE_FREE(grid->Rv);
    SAFE_FREE(grid->Rs);
    SAFE_FREE(grid->dppot);
    SAFE_FREE(grid->drpot);
    return;
}


void get_excited_torques(int mstart, int mend, double *TL, double *TR, Grid *grid, Params params, Disk *disk) {
    //Grid *grid = (Grid *)malloc(sizeof(Grid));
    //init_grid(mstart,mend, grid,params,planet,disk);
    int i;
    int num_modes = mend-mstart+1;
    for(i=0;i<num_modes;i++) {
        linearwaves(i, grid, params,disk,FALSE,FALSE);
        *TL += grid->TL[i];
        *TR += grid->TR[i];
    }
    //free_linear_grid(grid);

    return;

}
