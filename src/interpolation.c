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

void interpolate_onto_grid(double *xd, double *yd, int ndata, double *lr, double *ynew, double *dlydlr, double *d2lydlr, int n,int add_flag) {
    int i;
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

    interpolation(xd,yd,ndata,&lr[jstart],&ynew[jstart],&dlydlr[jstart],&d2lydlr[jstart],n_interp);


    
    for(i=0;i<jstart;i++) {
        ynew[i] = ynew[jstart] + dlydlr[jstart]*(lr[i]-lr[jstart]);
        dlydlr[i] = dlydlr[jstart];
        d2lydlr[i] = 0;
    }
    for(i=jend+1;i<n;i++) {
        ynew[i] = ynew[jend] + dlydlr[jend]*(lr[i]-lr[jend]);
        dlydlr[i] = dlydlr[jend];
        d2lydlr[i] = 0;
    }

// Add in inner truncation and outer truncation
    if (add_flag) {
        double ri = exp(lr[0])*.99;
        double ro = exp(lr[n-1])/2.0;
        double x;
        for(i=0;i<n;i++) {
            x = exp(lr[i]);
            ynew[i] += log(1 - sqrt(ri/x)) - pow(x/ro,2);
            dlydlr[i] += .5*sqrt(ri/x)*pow(1-sqrt(ri/x),-1) - 2*pow(x/ro,2);
            d2lydlr[i] += -.25*sqrt(ri/x)*pow(1-sqrt(ri/x),-2) - 4*pow(x/ro,2);
        }  
    }
    return;

}

void read_sigma(char *fname, double *lr, double *sigma, double *dlsdlr, double *d2lsdlr, double *omega, double *dlomdlr, double *d2lomdlr, int n,int readomega) {
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

    interpolate_onto_grid(xd,yd,ndata,lr,sigma,dlsdlr,d2lsdlr,n,TRUE);


    if (readomega) { 
        fres =  fread(yd,sizeof(double),ndata,f);
        FREAD_CHECK(fres,ndata);
        for(i=0;i<(int)ndata;i++) {
            yd[i] = log( exp(yd[i])/pow(exp(xd[i]),2));
        }
        interpolate_onto_grid(xd,yd,ndata,lr,omega,dlomdlr,d2lomdlr,n,FALSE);
    }


    fclose(f);
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
