#include "linearwaves.h"
#include <mpi.h>


int main(int argc, char *argv[]) {
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&np);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    init_planet();
    init_params();
    int mstart = atoi(argv[1]);
    int mend = atoi(argv[2]);
    params.nm = mend-mstart+1;
    int num_modes = params.nm/np;

    Grid *grid = (Grid *)malloc(sizeof(Grid));
    init_grid(num_modes, grid);
    grid->n = params.n;
    params.dlr = log(params.rmax/params.rmin) / (double)params.n;
    int i;
    for(i=0;i<params.n;i++) {
        grid->lr[i] = log(params.rmin) + params.dlr*i;
        grid->r[i] = exp(grid->lr[i]);
    }
    init_disk(argv[3],grid->lr);

    grid->n = params.n;
    grid->nm = num_modes;
    for(i=0;i<num_modes;i++) {
        grid->mvals[i] = rank*num_modes + mstart + i;
    }
    double *dp_pot_full = (double *)malloc(sizeof(double)*params.n*(params.nm));
    double *dr_pot_full = (double *)malloc(sizeof(double)*params.n*(params.nm));
/* Do FFT of potentail on root */
    if (rank == 0) {
        fft_potential(grid->r,dp_pot_full,dr_pot_full,params.nm);
    }
    MPI_Scatter(dp_pot_full,grid->n*num_modes,MPI_DOUBLE,grid->dppot,grid->n*num_modes,MPI_DOUBLE,0, MPI_COMM_WORLD);

    MPI_Scatter(dr_pot_full,grid->n*num_modes,MPI_DOUBLE,grid->drpot,grid->n*num_modes,MPI_DOUBLE,0, MPI_COMM_WORLD);
    free(dp_pot_full);
    free(dr_pot_full);
    for(i=0;i<num_modes;i++) {
        linearwaves(i, grid);
    }

    char fname[256];
    sprintf(fname,"%s.%d",argv[4],rank);
    output_torques(fname,grid);
    if (rank == 0) {
        sprintf(fname,"%s.disk",argv[4]);
        output_disk(fname);
    }

    free_grid(grid);
    free(grid);
    free_disk();

    int mpi_status  =  MPI_Finalize();

    return mpi_status;
}  
