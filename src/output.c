#include "linearwaves.h"


void output(double *r, double complex *sol, char *fname) {

    FILE *f = fopen(fname,"w");

    int n = params.n;
    int m = 3;
    int i;
    double nout[1] = {(double)n};
    fwrite(nout,sizeof(double),1,f);
    fwrite(r,sizeof(double),n,f);

    double rval;
    double ival;
    for(i=0;i<n;i++) {
        rval = creal(sol[i*m]);
        fwrite(&rval,sizeof(double),1,f);
    }
    for(i=0;i<n;i++) {
        ival = cimag(sol[i*m]);
        fwrite(&ival,sizeof(double),1,f);
    }
    for(i=0;i<n;i++) {
        rval = creal(sol[i*m+1]);
        fwrite(&rval,sizeof(double),1,f);
    }
    for(i=0;i<n;i++) {
        ival = cimag(sol[i*m+1]);
        fwrite(&ival,sizeof(double),1,f);
    }
    for(i=0;i<n;i++) {
        rval = creal(sol[i*m+2]);
        fwrite(&rval,sizeof(double),1,f);
    }
    for(i=0;i<n;i++) {
        ival = cimag(sol[i*m+2]);
        fwrite(&ival,sizeof(double),1,f);
    }
    return;
}
