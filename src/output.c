#include "linearwaves.h"


void output_disk(char *fname) {
    FILE *f = fopen(fname,"w");
    double temp;
    temp = (int)params.n;
    fwrite(&temp,sizeof(double),1,f);
    fwrite(disk.sigma, sizeof(double),params.n,f);
    fwrite(disk.dlsdlr, sizeof(double),params.n,f);
    fwrite(disk.d2lsdlr, sizeof(double),params.n,f);
    fwrite(disk.c2, sizeof(double),params.n,f);
    fwrite(disk.omega, sizeof(double),params.n,f);
    fwrite(disk.dlomdlr, sizeof(double),params.n,f);
    fwrite(disk.dlnudlr, sizeof(double),params.n,f);
    fwrite(disk.dlTdlr, sizeof(double),params.n,f);
    fwrite(disk.d2lTdlr, sizeof(double),params.n,f);
    fwrite(disk.kappa2, sizeof(double),params.n,f);
    fwrite(disk.nu, sizeof(double),params.n,f);
    fwrite(disk.pres, sizeof(double),params.n,f);
    fwrite(disk.dpdr, sizeof(double),params.n,f);

    fclose(f);
}

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
void output_torques(char *fname, Grid *grid) {
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
    //fwrite((double *)grid->mvals,sizeof(double),grid->nm,f);
    fwrite(grid->r,sizeof(double),params.n,f);
    fwrite(grid->lamex,sizeof(double),params.n*grid->nm,f);
    fwrite(grid->lamdep,sizeof(double),params.n*grid->nm,f);
    fwrite(grid->drfw,sizeof(double),params.n*grid->nm,f);
    fwrite(grid->fw,sizeof(double),params.n*grid->nm,f);
    fwrite(grid->dppot,sizeof(double),params.n*grid->nm,f);
    double rval;
    double ival;
    int j,i;
    /*
    for(j=0;j<grid->nm;j++) {
        for(i=0;i<grid->n;i++) {
            rval = creal(grid->u[i + j*grid->n]);
            fwrite(&rval,sizeof(double),1,f);
        }
        for(i=0;i<grid->n;i++) {
            ival = cimag(grid->u[i + j*grid->n]);
            fwrite(&ival,sizeof(double),1,f);
        }
    }
    for(j=0;j<grid->nm;j++) {
        for(i=0;i<grid->n;i++) {
            rval = creal(grid->v[i + j*grid->n]);
            fwrite(&rval,sizeof(double),1,f);
        }
        for(i=0;i<grid->n;i++) {
            ival = cimag(grid->v[i + j*grid->n]);
            fwrite(&ival,sizeof(double),1,f);
        }
    }
    for(j=0;j<grid->nm;j++) {
        for(i=0;i<grid->n;i++) {
            rval = creal(grid->s[i + j*grid->n]);
            fwrite(&rval,sizeof(double),1,f);
        }
        for(i=0;i<grid->n;i++) {
            ival = cimag(grid->s[i + j*grid->n]);
            fwrite(&ival,sizeof(double),1,f);
        }
    }
    */
    fclose(f);



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
