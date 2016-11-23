#include "linearwaves.h"
#ifdef _MPI
#include <mpi.h>
#endif

int main(int argc, char *argv[]) {
#ifdef _MPI
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&np);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);    
#else
    np = 1;
    rank = 0;
#endif
    char parfile[256];
    if (argc < 2) {
        strcpy(parfile,"in/in.par");
    }
    else {
        strcpy(parfile,argv[1]);
    }
    printf("Reading parameters from %s...\n",parfile);
    if (argc > 2) {
        read_param_file(parfile,argc-2,&argv[2]);
    }
    else {
        read_param_file(parfile,0,NULL);
    }
    printf("Read file\n");
    int mstart = params.mstart;
    int mend = params.mend;
    params.nm = mend-mstart+1;
    int num_modes = params.nm/np;
    printf("%lg\t%lg\n",params.rmin,params.rmax);

    Grid *grid = (Grid *)malloc(sizeof(Grid));
    init_grid(num_modes, grid);
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
    free(dp_pot_full);
    free(dr_pot_full);
    for(i=0;i<num_modes;i++) {
        linearwaves(i, grid);
    }

    char fname[256];
    sprintf(fname,"%s.%d",params.outputname,rank);
    output_torques(fname,grid);
    if (rank == 0) {
        sprintf(fname,"%s.disk",params.outputname);
        output_disk(fname);
        printf("Finished\n");
    }
    /*
    printf("free grid\n");
    free_grid(grid);
    printf("free grid\n");
    free(grid);
    printf("free grid\n");
    free_disk();
    printf("free grid\n");
*/
#ifdef _MPI
    int mpi_status  =  MPI_Finalize();
#else
    int mpi_status = 0;
#endif
    return mpi_status;
}  
