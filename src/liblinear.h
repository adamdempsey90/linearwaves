#include "defines.h"
#include "structs.h"
void interpolate_onto_grid(double *xd, double *yd, int ndata, double *lr, double *sigma, double *dlsdlr, double *d2lsdlr, int n);
void linearwaves(int i, Grid *grid, Params params, Disk *disk, int silent) ;
void init_grid(int mstart, int mend, Grid *grid, Params params, Disk *disk) ;
void get_excited_torques(int mstart, int mend, double *TL, double *TR, Grid *grid, Params params, Disk *disk) ;
void free_linear_grid(Grid *grid) ;
