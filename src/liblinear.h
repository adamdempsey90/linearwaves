#include "defines.h"
#include "structs.h"
void linearwaves(int i, Grid *grid, Params params, Planet planet, Disk *disk) ;
void init_grid(int mstart, int mend, Grid *grid, Params params, Planet planet, Disk *disk) ;
void get_excited_torques(int mstart, int mend, double *TL, double *TR, Grid *grid, Params params, Planet planet, Disk *disk) ;
void free_grid(Grid *grid) ;
