#include "linearwaves.h"
#include <mpi.h>

int np, rank;

int main(int argc, char *argv[]) {
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&np);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    init_planet();
    init_params();
    int mstart = atoi(argv[1]);
    int mend = atoi(argv[2]);
    int num_modes = (mend-mstart)/np;

    Grid *grid = (Grid *)malloc(sizeof(Grid));
    init_grid(num_modes, grid);
    int i;

    grid->n = params.n;
    grid->nm = num_modes;
    for(i=0;i<num_modes;i++) {
        grid->mvals[i] = rank*num_modes + mstart + i;
    }
    params.dlr = log(params.rmax/params.rmin) / (double)params.n;
    for(i=0;i<params.n;i++) {
        grid->r[i] = params.rmin * exp(i*params.dlr);
    }

    for(i=0;i<num_modes;i++) {
        //printf("%d\t%d\n",rank,i);
        linearwaves(i, grid);
    }

    char fname[256];
    sprintf(fname,"%s.%d",argv[3],rank);
    output_torques(fname,grid);

    free_grid(grid);
    free(grid);

    int mpi_status  =  MPI_Finalize();

    return mpi_status;
}  
