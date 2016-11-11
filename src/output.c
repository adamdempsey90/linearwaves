#include "linearwaves.h"


void output(double *r, double complex *sol, double *lamex, double *lamdep, double *drfw, double *fw, char *fname) {
 
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

    fwrite(lamex,sizeof(double),params.n,f);
    fwrite(lamdep,sizeof(double),params.n,f);
    fwrite(drfw,sizeof(double),params.n,f);
    fwrite(fw,sizeof(double),params.n,f);
    fclose(f);
    return;
}

void output_matrix(double complex *ld, double complex *md, double complex *ud, double complex *fd) {

    FILE *f = fopen("outputs/matrix.dat","w");

    int i,j;
    double rpart, ipart;
    for(i=0;i<params.n;i++) {
        for(j=0;j<9;j++) {
            rpart = creal(md[i*9 + j]);
            ipart = cimag(md[i*9 + j]);
            fwrite(&rpart,sizeof(double),1,f);
            fwrite(&ipart,sizeof(double),1,f);
        }
    }

    //fwrite((double *)ld,sizeof(double),(params.n-1)*params.nrhs*params.nrhs*2,f);
    //fwrite((double *)ud,sizeof(double),(params.n-1)*params.nrhs*params.nrhs*2,f);
    //fwrite((double *)md,sizeof(double),params.n*params.nrhs*params.nrhs*2,f);
    //fwrite((double *)fd,sizeof(double),params.n*params.nrhs*2,f);
    fclose(f);
    return;
}
