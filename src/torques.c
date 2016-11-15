#include "linearwaves.h"


void calc_torques(double *r, double *fw, double *drfw, double *lamex, double *lamdep, double complex *sol, double *dppot, double *TL, double *TR,int m) {
    double complex lf[3], uf[3], mf[3];
    double invdlr, invdlr2;
    double om;
    double k2om;
    double complex res, u, v, s;
    double complex du, dv;
    int i,j;
    int size = params.nrhs;
    FILE *f;
    char fname[256];
    sprintf(fname,"outputs/sol%d.dat.%d",m,rank);
    f = fopen(fname,"w");
    for(j=0;j<params.n;j++) {
       fprintf(f,"%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\n",creal(sol[j*3]),cimag(sol[j*3]),
               creal(sol[j*3+1]),cimag(sol[j*3+1]),
               creal(sol[j*3+2]),cimag(sol[j*3+2]));
    }
    invdlr = 1./(params.dlr);
    invdlr2 = invdlr*invdlr;
    *TL = 0;
    *TR = 0;
    for(i=1;i<params.n-1;i++) {
        for (j=0;j<params.nrhs;j++) {
            mf[j] = 0;
            lf[j] = 0;
            uf[j] = 0;
        }
        lamdep[i] = 0;

        u = sol[i*params.nrhs];
        v = sol[i*params.nrhs+1];
        dv = sol[(i+1)*params.nrhs+1] - sol[(i-1)*params.nrhs+1];
        du = sol[(i+1)*params.nrhs] - sol[(i-1)*params.nrhs];
        s = sig(r[i],sol[i*params.nrhs+2]);



        om = omega(r[i]);
        k2om = kappa2(r[i])/(2*om);

        viscosity_coeffs_v(r[i],mf,m);
        viscosity_d2coeffs_v(r[i],mf,m,-2*invdlr2);
        viscosity_dcoeffs_v(r[i],lf,m,-.5*invdlr);
        viscosity_d2coeffs_v(r[i],lf,m,invdlr2);
        viscosity_dcoeffs_v(r[i],uf,m,.5*invdlr);
        viscosity_d2coeffs_v(r[i],uf,m,invdlr2);

        res = 0;
        for(j=0;j<params.nrhs;j++) {
            res += lf[j]*sol[(i-1)*size+j] + mf[j]*sol[i*size+j] + uf[j]*sol[(i+1)*size+j]; 
        }
        res *= r[i];
        lamdep[i] += 2*creal(conj(s) * res);
        lamdep[i] += 2*creal(conj(s)*u)*r[i]*k2om;
        lamdep[i] -= sigma(r[i])*2*creal(u*conj(dv*.5*invdlr + v));
        lamdep[i] += params.ieps*2*creal(u*conj(s))*r[i];

        fw[i] = r[i]*r[i]*sigma(r[i])*2*creal(u * conj(v));
        drfw[i] = (2*sigma(r[i]) + r[i]*dsdr(r[i]))*2*creal(u*conj(v)) + .5*invdlr*sigma(r[i])*2*creal(du*conj(v)+u*conj(dv));
        for(j=0;j<params.nrhs;j++) mf[j] = 0;
        //lamex[i] = -2*m*creal(cimag(s)*dppot[i]);
        lamex[i] = r[i]*2*creal(conj(s)*-dppot[i]*I*m/r[i]);
        if (r[i] >= planet.a) *TR += lamex[i]*2*M_PI*r[i]*r[i]*params.dlr;
        if (r[i] <= planet.a) *TL -= lamex[i]*2*M_PI*r[i]*r[i]*params.dlr;
    }

    lamdep[0] = 0;
    lamdep[params.n-1] = 0;

    fw[0] = 0;
    drfw[0] = 0;
    lamex[0] = 0;
    fw[params.n-1] = 0;
    drfw[params.n-1] = 0;
    lamex[params.n-1] = 0;


    return;

}
