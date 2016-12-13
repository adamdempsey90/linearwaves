#include "linearwaves.h"


void output_disk(char *fname,Params params, Disk *disk) {
    FILE *f = fopen(fname,"w");
    double temp;
    temp = (int)params.n;
    fwrite(&temp,sizeof(double),1,f);
    fwrite(disk->sigma, sizeof(double),params.n,f);
    fwrite(disk->dlsdlr, sizeof(double),params.n,f);
    fwrite(disk->d2lsdlr, sizeof(double),params.n,f);
    fwrite(disk->omega, sizeof(double),params.n,f);
    fwrite(disk->dlomdlr, sizeof(double),params.n,f);
    fwrite(disk->kappa2, sizeof(double),params.n,f);
    fwrite(disk->c2, sizeof(double),params.n,f);
    fwrite(disk->dlnudlr, sizeof(double),params.n,f);
    fwrite(disk->dlTdlr, sizeof(double),params.n,f);
    fwrite(disk->d2lTdlr, sizeof(double),params.n,f);
    fwrite(disk->nu, sizeof(double),params.n,f);
    fwrite(disk->pres, sizeof(double),params.n,f);
    fwrite(disk->dpdr, sizeof(double),params.n,f);

    fclose(f);
}

void output(double *r, double complex *sol, double *lamex, double *lamdep, double *drfw, double *fw, char *fname, Params params) {
 
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
void output_torques(char *fname, Params params,Grid *grid) {
    FILE *f = fopen(fname,"w");

    printf("Outputting to %s\n",fname);
    double n = (double)params.n;
    double m = (double)grid->nm;
    double mstart = (double)grid->mvals[0];
    double mend = (double)grid->mvals[grid->nm-1];
    fwrite(&n,sizeof(double),1,f);
    fwrite(&m,sizeof(double),1,f);
    fwrite(&mstart,sizeof(double),1,f);
    fwrite(&mend,sizeof(double),1,f);
    fwrite(grid->TL,sizeof(double),grid->nm,f);
    fwrite(grid->TR,sizeof(double),grid->nm,f);
    fwrite(grid->r,sizeof(double),params.n,f);
    fwrite(grid->lamex,sizeof(double),params.n*grid->nm,f);
    fwrite(grid->lamdep,sizeof(double),params.n*grid->nm,f);
    fwrite(grid->drfw,sizeof(double),params.n*grid->nm,f);
    fwrite(grid->fw,sizeof(double),params.n*grid->nm,f);
    fwrite(grid->dppot,sizeof(double),params.n*grid->nm,f);

    fclose(f);


}

void output_matrix(double complex *ld, double complex *md, double complex *ud, double complex *fd, Params params) {

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
    for(i=0;i<params.n-1;i++) {
        for(j=0;j<9;j++) {
            rpart = creal(ld[i*9 + j]);
            ipart = cimag(ld[i*9 + j]);
            fwrite(&rpart,sizeof(double),1,f);
            fwrite(&ipart,sizeof(double),1,f);
        }
    }
    for(i=0;i<params.n-1;i++) {
        for(j=0;j<9;j++) {
            rpart = creal(ud[i*9 + j]);
            ipart = cimag(ud[i*9 + j]);
            fwrite(&rpart,sizeof(double),1,f);
            fwrite(&ipart,sizeof(double),1,f);
        }
    }
    for(i=0;i<params.n;i++) {
        for(j=0;j<3;j++) {
            rpart = creal(fd[i*3 + j]);
            ipart = cimag(fd[i*3 + j]);
            fwrite(&rpart,sizeof(double),1,f);
            fwrite(&ipart,sizeof(double),1,f);
        }
    }

    fclose(f);
    return;
}
