#include "defines.h"
#include "structs.h"
void interpolate_onto_grid(double *xd, double *yd, int ndata, double *lr, double *sigma, double *dlsdlr, double *d2lsdlr, int n, int add_flag);
void linearwaves(int i, Grid *grid, Params params, Disk *disk, int silent, int second_order) ;
void init_grid(int mstart, int mend, Grid *grid, Params params, Disk *disk) ;
void get_excited_torques(int mstart, int mend, double *TL, double *TR, Grid *grid, Params params, Disk *disk) ;
void free_linear_grid(Grid *grid) ;
void calc_forcing(double *r, double complex *u, double complex *v, double complex *s, double complex *Rr, double complex *Rp, double complex *Rs, Params params, Disk *disk);
