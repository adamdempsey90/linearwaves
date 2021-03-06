#include "linearwaves.h"

int main(int argc, char *argv[]) {
#ifdef _MPI
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&proc.np);
    MPI_Comm_rank(MPI_COMM_WORLD,&proc.rank);    
#else
    proc.np = 1;
    proc.rank = 0;
#endif
    char parfile[256];
    if (argc < 2) {
        strcpy(parfile,"in/in.par");
    }
    else {
        strcpy(parfile,argv[1]);
    }
    printf("Reading parameters from %s...\n",parfile);


    Params params;
    Disk *disk = (Disk *)malloc(sizeof(Disk));
    if (argc > 2) {
        read_param_file(parfile,argc-2,&argv[2],&params);
    }
    else {
        read_param_file(parfile,0,NULL,&params);
    }
    printf("Read file\n");
    int mstart = params.mstart;
    int mend = params.mend;
    params.nm = mend-mstart+1;
    int num_modes = params.nm/proc.np;
    printf("%lg\t%lg\n",params.rmin,params.rmax);

    
    Grid *grid = (Grid *)malloc(sizeof(Grid));
    init_grid(mstart,mend, grid,params,disk);
/*
    double TL, TR;
    get_excited_torques(mstart,mend,&TL, &TR, grid,params,disk);
    printf("Final %lg\t%lg\n",TL,TR);
*/

    int i,j;

    for(i=0;i<num_modes;i++) {
        linearwaves(i, grid,params,disk,FALSE,FALSE);
    }



    char fname[256];
    sprintf(fname,"%s.%d",params.outputname,proc.rank);
    output_torques(fname,params,grid);

    if (params.second_order) {
        printf("Calculating second order forcing...\n");
        calc_forcing(grid->r,grid->u,grid->v,grid->s,grid->Ru,grid->Rv,grid->Rs,params,disk);
        printf("Calculating second order solution...\n");
        for(i=0;i<num_modes;i++) {
            linearwaves(i, grid,params,disk,FALSE,TRUE);
        }
    }
    if (proc.rank == 0) {
        sprintf(fname,"%s.disk",params.outputname);
        output_disk(fname,params,disk);
        printf("Finished\n");
    }
    /*
    printf("free grid\n");
    free_linear_grid(grid);
    printf("free grid\n");
    SAFE_FREE(grid);
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
