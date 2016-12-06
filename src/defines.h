#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#ifdef _MPI
#include <mpi.h>
#endif

#define TRUE 1
#define FALSE 0

#define SAFE_FREE(ptr) free(ptr); ptr = NULL;
