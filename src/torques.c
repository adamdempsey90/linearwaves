#include "linearwaves.h"


void calc_deposited_torque(double *r, double *lamdep, double complex *sol, int m) {
    double complex lf[3], uf[3], mf[3];
    double visc, ir2;
    double invdlr, invdlr2;
    double om;
    double k2om;
    double c2;
    double complex res, u, v, s;
    double norm;
    int i,j;
    int size = params.nrhs;
    for(i=1;i<params.n-1;i++) {
        for (j=0;j<params.nrhs;j++) {
            mf[j] = 0;
            lf[j] = 0;
            uf[j] = 0;
        }
        lamdep[i] = 0;
        om = omega(r[i]);
        k2om = kappa2(r[i])/(2*om);
        c2 = cs2(r[i]);
        invdlr = -.5/params.dlr;
        invdlr2 = -2./(params.dlr*params.dlr);
        visc = nu(r[i]);
        ir2 = 1./(r[i]*r[i]);
        norm = (params.iso) ? 1. : c2;

        u = sol[i*size];
        v = sol[i*size+1];
        s = sol[i*size+2];

        mf[0] += visc*I*m*(3 - params.eta + params.nuindx)*ir2;
        mf[1] -= visc*(1 + params.nuindx + m*m*(2 - params.eta))*ir2;
        mf[1] += visc*ir2*invdlr2;
        if (params.iso) {
            mf[2] += visc*I*m*(k2om-2*om)/r[i]*invdlr;
        }
        else {
             mf[2] += visc*I*m*(k2om-2*om)/(r[i]*c2)*invdlr;
        }
    /* uf */
        invdlr = .5/params.dlr;
        invdlr2 = 1./(params.dlr*params.dlr);
        uf[0] += visc*(1-params.eta)*I*m*ir2*invdlr;
        uf[1] += visc*params.nuindx*ir2*invdlr;
        uf[1] += visc*ir2*invdlr2;
        if (params.iso) {
            uf[2] += visc*I*m*(k2om-2*om)/r[i]*invdlr;
        }
        else {
             uf[2] += visc*I*m*(k2om-2*om)/(r[i]*c2)*invdlr;
        }

    /* lf */
        invdlr = -.5/params.dlr;
        invdlr2 = 1./(params.dlr*params.dlr);
        lf[0] += visc*(1-params.eta)*I*m*ir2*invdlr;
        lf[1] += visc*params.nuindx*ir2*invdlr;
        lf[1] += visc*ir2*invdlr2;
        if (params.iso) {
            lf[2] += visc*I*m*(k2om-2*om)/r[i] *invdlr;
        }
        else {
             lf[2] += visc*I*m*(k2om-2*om)/(r[i]*c2) *invdlr;
        }
        res = 0;
        for(j=0;j<params.nrhs;j++) {
            res += lf[j]*sol[(i-1)*size+j] + mf[j]*sol[i*size+j] + uf[j]*sol[(i+1)*size+j]; 
        }
        res *= r[i];
        lamdep[i] -= 2*creal(conj(s/norm) * res);

        
        lamdep[i] += 2*creal(conj(s/norm)*u)*r[i]*k2om;
        res = (sol[(i+1)*size+1]*r[i+1] - r[i]*sol[(i-1)*size+1])/(2*params.dlr*r[i]);
        lamdep[i] -= 2*creal(conj(u)*res);
        lamdep[i] *= 2*M_PI*r[i];
    }

    lamdep[0] = 0;
    lamdep[params.n-1] = 0;


    return;

}
