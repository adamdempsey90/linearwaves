#include "linearwaves.h"

void calc_forcing(double *r, double complex *u, double complex *v, double complex *s, double complex *Rr, double complex *Rp, double complex *Rs, Params params, Disk *disk) {
    int i,j,k;
    int nm, nr;
    double invdlr = .5/(params.dlr);
    /* Three sums */
    
    double x,dbar;
    
    nr = params.n;
    nm = params.nm;
    


    for(k=1;k<nr-1;k++) {
        x = r[k];
        dbar = disk->sigma[k];
        for(i=0;i<nm;i++) {
            Rr[k+i*nr] = 0;
            Rp[k+i*nr] = 0;
            Rs[k+i*nr] = 0;
        /* j = 1 - i */
            for(j=0;j<=i;j++) {
                Rr[k+i*nr] -= (.5/x)*( u[k+j*nr]*invdlr*(u[k+1+(i-j)*nr]-u[k-1+(i-j)*nr]) + I*(i-j)*v[k+j*nr]*u[k+(i-j)*nr] - v[k+j*nr]*v[k+(i-j)*nr]);
                Rp[k+i*nr] -= (.5/x)*( u[k+j*nr]*invdlr*(v[k+1+(i-j)*nr]-v[k-1+(i-j)*nr]) + I*(i-j)*v[k+j*nr]*v[k+(i-j)*nr] + u[k+j*nr]*v[k+(i-j)*nr]);
                Rs[k+i*nr] -= (.5/(dbar*x))*( s[k+j*nr]*u[k+(i-j)*nr] + (u[k+(i-j)*nr]*(s[k+1+j*nr]-s[k-1+j*nr]) + (u[k+1+(i-j)*nr]-u[k-1+(i-j)*nr])*s[k+j*nr])*invdlr + I*(i+1)*s[k+j*nr]*v[k+(i-j)*nr]);
            }
        /* j = i+1 - infty */
            for(j=i+1;j<nm;j++) {
                Rr[k+i*nr] -= (.5/x)*( u[k+j*nr]*conj((u[k+1+(j-i)*nr]-u[k-1+(j-i)*nr]))*invdlr + I*(j-i)*v[k+j*nr]*conj(u[k+(j-i)*nr]) - v[k+j*nr]*conj(v[k+(j-i)*nr]));
                Rp[k+i*nr] -= (.5/x)*( u[k+j*nr]*conj((v[k+1+(j-i)*nr]-v[k-1+(j-i)*nr]))*invdlr + I*(j-i)*v[k+j*nr]*conj(v[k+(j-i)*nr]) + u[k+j*nr]*conj(v[k+(j-i)*nr]));
                Rs[k+i*nr] -= (.5/(dbar*x))*( s[k+j*nr]*conj(u[k+(j-i)*nr]) + (conj(u[k+(j-i)*nr])*(s[k+1+j*nr]-s[k-1+j*nr]) + conj((u[k+1+(j-i)*nr]-u[k-1+(j-i)*nr]))*s[k+j*nr])*invdlr - I*(i+1)*s[k+j*nr]*conj(v[k+(j-i)*nr]));


            }
        /* j = 1 - infy */
            for(j=0;j<nm;j++) {
                Rr[k+i*nr] -= (.5/x)*( conj(u[k+j*nr])*(u[k+1+(i+j)*nr]-u[k-1+(i+j)*nr])*invdlr + I*(i+j+2)*conj(v[k+j*nr])*u[k+(i+j)*nr] - conj(v[k+j*nr])*v[k+(i+j)*nr]);
                Rp[k+i*nr] -= (.5/x)*( conj(u[k+j*nr])*(v[k+1+(i+j)*nr]-v[k-1+(i+j)*nr])*invdlr + I*(i+j+2)*conj(v[k+j*nr])*v[k+(i+j)*nr] - conj(u[k+j*nr])*v[k+(i+j)*nr]);
                Rs[k+i*nr] -= (.5/(dbar*x))*( conj(s[k+j*nr])*u[k+(i+j)*nr] + (u[k+(i+j)*nr]*conj((s[k+1+j*nr]-s[k-1+j*nr])) + (u[k+1+(i+j)*nr]-u[k-1+(i+j)*nr])*conj(s[k+j*nr]))*invdlr - I*(i+1)*conj(s[k+j*nr])*v[k+(i+j)*nr]);

            }
        }
    }

    k = 0;
    x = r[k];
    dbar = disk->sigma[k];
    double c2 = disk->c2[k];
    double rb = x * exp(.5*params.dlr);
    double kr;
    double complex *fac =(double complex *)malloc(sizeof(double complex)*nm);
    for(i=0;i<nm;i++) {
        kr = rb*pow( fabs(Dfunc(k,params.omf,i+1,params,disk)/c2) ,.5);   
        fac[i] = (1 - .5*I*kr*params.dlr)/(1 + .5*I*kr*params.dlr);
    }
    for(i=0;i<nm;i++) {
        Rr[k+i*nr] = 0;
        Rp[k+i*nr] = 0;
        Rs[k+i*nr] = 0;
    /* j = 1 - i */
        for(j=0;j<=i;j++) {
            Rr[k+i*nr] -= (.5/x)*( u[k+j*nr]*invdlr*(u[k+1+(i-j)*nr]-fac[i-j]*u[k+(i-j)*nr]) + I*(i-j)*v[k+j*nr]*u[k+(i-j)*nr] - v[k+j*nr]*v[k+(i-j)*nr]);
            Rp[k+i*nr] -= (.5/x)*( u[k+j*nr]*invdlr*(v[k+1+(i-j)*nr]-fac[i-j]*v[k+(i-j)*nr]) + I*(i-j)*v[k+j*nr]*v[k+(i-j)*nr] + u[k+j*nr]*v[k+(i-j)*nr]);
            Rs[k+i*nr] -= (.5/(dbar*x))*( s[k+j*nr]*u[k+(i-j)*nr] + (u[k+(i-j)*nr]*(s[k+1+j*nr]-fac[j]*s[k+j*nr]) + (u[k+1+(i-j)*nr]-fac[i-j]*u[k+(i-j)*nr])*s[k+j*nr])*invdlr + I*(i+1)*s[k+j*nr]*v[k+(i-j)*nr]);
        }
    /* j = i+1 - infty */
        for(j=i+1;j<nm;j++) {
            Rr[k+i*nr] -= (.5/x)*( u[k+j*nr]*conj((u[k+1+(j-i)*nr]-fac[j-i]*u[k+(j-i)*nr]))*invdlr + I*(j-i)*v[k+j*nr]*conj(u[k+(j-i)*nr]) - v[k+j*nr]*conj(v[k+(j-i)*nr]));
            Rp[k+i*nr] -= (.5/x)*( u[k+j*nr]*conj((v[k+1+(j-i)*nr]-fac[j-i]*v[k+(j-i)*nr]))*invdlr + I*(j-i)*v[k+j*nr]*conj(v[k+(j-i)*nr]) + u[k+j*nr]*conj(v[k+(j-i)*nr]));
            Rs[k+i*nr] -= (.5/(dbar*x))*( s[k+j*nr]*conj(u[k+(j-i)*nr]) + (conj(u[k+(j-i)*nr])*(s[k+1+j*nr]-fac[j]*s[k+j*nr]) + conj((u[k+1+(j-i)*nr]-fac[j-i]*u[k+(j-i)*nr]))*s[k+j*nr])*invdlr - I*(i+1)*s[k+j*nr]*conj(v[k+(j-i)*nr]));


        }
    /* j = 1 - infy */
        for(j=0;j<nm;j++) {
            Rr[k+i*nr] -= (.5/x)*( conj(u[k+j*nr])*(u[k+1+(i+j)*nr]-fac[i+j]*u[k+(i+j)*nr])*invdlr + I*(i+j+2)*conj(v[k+j*nr])*u[k+(i+j)*nr] - conj(v[k+j*nr])*v[k+(i+j)*nr]);
            Rp[k+i*nr] -= (.5/x)*( conj(u[k+j*nr])*(v[k+1+(i+j)*nr]-fac[i+j]*v[k+(i+j)*nr])*invdlr + I*(i+j+2)*conj(v[k+j*nr])*v[k+(i+j)*nr] - conj(u[k+j*nr])*v[k+(i+j)*nr]);
            Rs[k+i*nr] -= (.5/(dbar*x))*( conj(s[k+j*nr])*u[k+(i+j)*nr] + (u[k+(i+j)*nr]*conj((s[k+1+j*nr]-fac[j]*s[k+j*nr])) + (u[k+1+(i+j)*nr]-fac[i+j]*u[k+(i+j)*nr])*conj(s[k+j*nr]))*invdlr - I*(i+1)*conj(s[k+j*nr])*v[k+(i+j)*nr]);

        }
    }
    k = nr-1;
    x = r[k];
    dbar = disk->sigma[k];
    c2 = disk->c2[k];
    rb = x * exp(.5*params.dlr);
    for(i=0;i<nm;i++) {
        kr = rb*pow( fabs(Dfunc(k,params.omf,i+1,params,disk)/c2) ,.5);   
        fac[i] = (1 + .5*I*kr*params.dlr)/(1 - .5*I*kr*params.dlr);
    }
    for(i=0;i<nm;i++) {
        Rr[k+i*nr] = 0;
        Rp[k+i*nr] = 0;
        Rs[k+i*nr] = 0;
    /* j = 1 - i */
        for(j=0;j<=i;j++) {
            Rr[k+i*nr] -= (.5/x)*( u[k+j*nr]*invdlr*(fac[i-j]*u[k+(i-j)*nr]-u[k-1+(i-j)*nr]) + I*(i-j)*v[k+j*nr]*u[k+(i-j)*nr] - v[k+j*nr]*v[k+(i-j)*nr]);
            Rp[k+i*nr] -= (.5/x)*( u[k+j*nr]*invdlr*(fac[i-j]*v[k+(i-j)*nr]-v[k-1+(i-j)*nr]) + I*(i-j)*v[k+j*nr]*v[k+(i-j)*nr] + u[k+j*nr]*v[k+(i-j)*nr]);
            Rs[k+i*nr] -= (.5/(dbar*x))*( s[k+j*nr]*u[k+(i-j)*nr] + (u[k+(i-j)*nr]*(fac[j]*s[k+j*nr]-s[k-1+j*nr]) + (fac[i-j]*u[k+(i-j)*nr]-u[k-1+(i-j)*nr])*s[k+j*nr])*invdlr + I*(i+1)*s[k+j*nr]*v[k+(i-j)*nr]);
        }
    /* j = i+1 - infty */
        for(j=i+1;j<nm;j++) {
            Rr[k+i*nr] -= (.5/x)*( u[k+j*nr]*conj((fac[j-i]*u[k+(j-i)*nr]-u[k-1+(j-i)*nr]))*invdlr + I*(j-i)*v[k+j*nr]*conj(u[k+(j-i)*nr]) - v[k+j*nr]*conj(v[k+(j-i)*nr]));
            Rp[k+i*nr] -= (.5/x)*( u[k+j*nr]*conj((fac[j-i]*v[k+(j-i)*nr]-v[k-1+(j-i)*nr]))*invdlr + I*(j-i)*v[k+j*nr]*conj(v[k+(j-i)*nr]) + u[k+j*nr]*conj(v[k+(j-i)*nr]));
            Rs[k+i*nr] -= (.5/(dbar*x))*( s[k+j*nr]*conj(u[k+(j-i)*nr]) + (conj(u[k+(j-i)*nr])*(fac[j]*s[k+j*nr]-s[k-1+j*nr]) + conj((fac[j-i]*u[k+(j-i)*nr]-u[k-1+(j-i)*nr]))*s[k+j*nr])*invdlr - I*(i+1)*s[k+j*nr]*conj(v[k+(j-i)*nr]));


        }
    /* j = 1 - infy */
        for(j=0;j<nm;j++) {
            Rr[k+i*nr] -= (.5/x)*( conj(u[k+j*nr])*(fac[i+j]*u[k+(i+j)*nr]-u[k-1+(i+j)*nr])*invdlr + I*(i+j+2)*conj(v[k+j*nr])*u[k+(i+j)*nr] - conj(v[k+j*nr])*v[k+(i+j)*nr]);
            Rp[k+i*nr] -= (.5/x)*( conj(u[k+j*nr])*(fac[i+j]*v[k+(i+j)*nr]-v[k-1+(i+j)*nr])*invdlr + I*(i+j+2)*conj(v[k+j*nr])*v[k+(i+j)*nr] - conj(u[k+j*nr])*v[k+(i+j)*nr]);
            Rs[k+i*nr] -= (.5/(dbar*x))*( conj(s[k+j*nr])*u[k+(i+j)*nr] + (u[k+(i+j)*nr]*conj((fac[j]*s[k+j*nr]-s[k-1+j*nr])) + (fac[i+j]*u[k+(i+j)*nr]-u[k-1+(i+j)*nr])*conj(s[k+j*nr]))*invdlr - I*(i+1)*conj(s[k+j*nr])*v[k+(i+j)*nr]);

        }
    }
}
