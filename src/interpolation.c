#include "linearwaves.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>


#define FREAD_CHECK(i,j) ((i) >= (j)) ? (void) 0 : printf("fread returned read fewer items than expected (%zu < %zu) at line %d in %s of %s\n",i,j,__LINE__,__func__,__FILE__)

void interpolation(double *xd, double *yd, int nd, double *x, double *y, double *dy, double *d2y, int n) {
    gsl_interp_accel *acc = gsl_interp_accel_alloc ();
    gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, nd);

    gsl_spline_init (spline, xd, yd, nd);

    int i;
    for(i=0;i<n;i++) {
        y[i] = gsl_spline_eval(spline,x[i],acc);
        dy[i] = gsl_spline_eval_deriv(spline,x[i],acc);
        d2y[i] = gsl_spline_eval_deriv2(spline,x[i],acc);
    }
    gsl_spline_free (spline);
    gsl_interp_accel_free (acc);

    return;
}

void read_sigma(char *fname, double *lr, double *sigma, double *dlsdlr, double *d2lsdlr, int n) {
    double temp;
    size_t ndata;
    double *xd,*yd;
    int i;
    size_t fres;

    FILE *f = fopen(fname,"r");
    ndata = 1;
    fres = fread(&temp,sizeof(double),ndata,f);
    FREAD_CHECK(fres,ndata);
    ndata=(size_t)temp;

    xd = (double *)malloc(sizeof(double)*ndata);
    yd = (double *)malloc(sizeof(double)*ndata);

    fres = fread(xd,sizeof(double),ndata,f);
    FREAD_CHECK(fres,ndata);
    fres =  fread(yd,sizeof(double),ndata,f);
    FREAD_CHECK(fres,ndata);
    fclose(f);


    int jstart = 0;
    int jend = n-1;

    for(i=0;i<n;i++) {
        if (lr[i] >= xd[0]) {
            jstart = i;
            break;
        }
    }
    for(i=jstart;i<n;i++) {
        if (lr[i] >= xd[ndata-1]) {
            if (lr[i] == xd[ndata-1]) {
                jend = i;
            }
            else {
                jend = i-1;
            }
            break;
        }
    }
    int n_interp = jend-jstart +1;

    interpolation(xd,yd,ndata,&lr[jstart],&sigma[jstart],&dlsdlr[jstart],&d2lsdlr[jstart],n_interp);

    for(i=0;i<jstart;i++) {
        sigma[i] = sigma[jstart] + dlsdlr[jstart]*(lr[i]-lr[jstart]);
        dlsdlr[i] = dlsdlr[jstart];
        d2lsdlr[i] = 0;
    }
    for(i=jend+1;i<n;i++) {
        sigma[i] = sigma[jend] + dlsdlr[jend]*(lr[i]-lr[jend]);
        dlsdlr[i] = dlsdlr[jend];
        d2lsdlr[i] = 0;
    }
    SAFE_FREE(xd);
    SAFE_FREE(yd);
    return;
}
/*
int main(void) {
    double *x,*y,*dy,*d2y;
    int i;

    double rmin = .01;
    double rmax = 10;
    int n = 1000;
    x = (double *)malloc(sizeof(double)*n);
    y = (double *)malloc(sizeof(double)*n);
    dy = (double *)malloc(sizeof(double)*n);
    d2y = (double *)malloc(sizeof(double)*n);

    for(i=0;i<n;i++) {
        x[i] = log(rmin) + i*log(rmax/rmin)/n;
        y[i]  = 0;
        dy[i] = 0;
        d2y[i] = 0;
    }

    read_sigma("../profiles/avg_profs/prof_q1e-3_a1e-2.bin",x,y,dy,d2y,n);

    FILE *f = fopen("interp_res.dat","w");

    double temp = (double)n;
    fwrite(&temp,sizeof(double),1,f);
    fwrite(x,sizeof(double),n,f);
    fwrite(y,sizeof(double),n,f);
    fwrite(dy,sizeof(double),n,f);
    fwrite(d2y,sizeof(double),n,f);
    SAFE_FREE(x);
    SAFE_FREE(y);
    SAFE_FREE(dy);
    SAFE_FREE(d2y);
}
*/
