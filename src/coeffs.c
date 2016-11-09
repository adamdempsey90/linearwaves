#include "linearwaves.h"


void main_diag(double r, int m, double complex *res) {
    double om = omega(r);
    double c2 = cs2(r);
    double k2om = kappa2(r)/(2*omega(r));
    double invdlr2 = -2./(params.dlr*params.dlr);


    res[0] = I*m*(om - params.omf);     
    res[1] = -2*om;                     
    res[2] = 0.0;                     
    res[3] = k2om;                     
    res[4] = I*m*(om - params.omf);     

    if (params.iso) {
        res[5] = I*m*c2/r ;
        res[6] = (params.mu+1)/r ;
        res[7] = I*m/r ;
    }
    else {
        res[5] = I*m*r;                     
        res[6] = (params.mu + 1)*c2/r;                     
        res[7] = I*m*c2/r;                     
    }
    res[8] = I*m*(om - params.omf);    

    double visc = -nu(r);
    double ir2 = 1./(r*r);

    res[0] -= visc*(2 + m*m + (1 - params.eta)*params.nuindx)*ir2;
    res[1] -= visc*I*m*(3 + (1 + params.nuindx)*params.eta)*ir2;
    res[3] += visc*I*m*(3 - params.eta + params.nuindx)*ir2;
    res[4] -= visc*(1 + params.nuindx + m*m*(2 - params.eta))*ir2;
    res[0] += visc*(2+params.eta)*ir2*invdlr2;
    res[4] += visc*ir2*invdlr2;
    if (params.iso) {
        res[2] += visc*1j*m*(k2om - 2*om)/r;
    }
    else {
        res[2] += visc*1j*m*(k2om - 2*om)/(r*c2) ;
        res[5] -= visc*1j*m*(k2om - 2*om)/(r*c2) * params.delta;
    }
    return;
}

void upper_diag(double r, int m, double complex *res) {

   double c2 = cs2(r);
   double om = omega(r);
   double k2om = kappa2(r)/(2*om);
   double invdlr = .5/params.dlr;
   double invdlr2 = 1./(params.dlr*params.dlr);
    res[0] = 0;
    res[1] = 0;
    res[3] = 0;
    res[4] = 0;
    res[5] = 0;
    res[7] = 0;
    res[8] = 0;
    if (params.iso) {
        res[2] = c2/r * invdlr;
        res[6] = invdlr / r;
    }
    else {
        res[2] = 1./r * invdlr;
        res[6] = c2*invdlr / r;
    }
    double visc = -nu(r);
    double ir2 = 1./(r*r);
    res[0] += visc*(2*params.nuindx+1)*ir2*invdlr;
    res[1] += visc*I*m*(1+params.eta)*ir2*invdlr;
    res[3] += visc*(1-params.eta)*I*m*ir2*invdlr;
    res[4] += visc*params.nuindx*ir2*invdlr;
    res[0] += visc*(2+params.eta)*ir2*invdlr2;
    res[4] += visc*ir2*invdlr2;
    if (params.iso) {
        res[5] += visc*I*m*(k2om-2*om)/r*invdlr;
    }
    else {
         res[5] += visc*I*m*(k2om-2*om)/(r*c2)*invdlr;
    }

    return;
}

void lower_diag(double r, int m, double complex *res) {
    double om = omega(r);
    double k2om = kappa2(r)/(2*om);
    double c2 = cs2(r);
    double invdlr = -.5/params.dlr;
    double invdlr2 = 1./(params.dlr*params.dlr);

    res[0] = 0;
    res[1] = 0;
    res[3] = 0;
    res[4] = 0;
    res[5] = 0;
    res[7] = 0;
    res[8] = 0;
    if (params.iso) {
        res[2] = c2/r * invdlr;
        res[6] = invdlr / r;
    }
    else {
        res[2] = 1./r * invdlr;
        res[6] = c2*invdlr / r;
    }
    
    double visc = -nu(r);
    double ir2 = 1./(r*r);

    res[0] += visc*(2*params.nuindx+1)*ir2*invdlr;
    res[1] += visc*I*m*(1+params.eta)*ir2*invdlr;
    res[3] += visc*(1-params.eta)*I*m*ir2*invdlr;
    res[4] += visc*params.nuindx*ir2*invdlr;
    res[0] += visc*(2+params.eta)*ir2*invdlr2;
    res[4] += visc*ir2*invdlr2;
    if (params.iso) {
        res[5] += visc*I*m*(k2om-2*om)/r *invdlr;
    }
    else {
         res[5] += visc*I*m*(k2om-2*om)/(r*c2) *invdlr;
    }
    


    return;
}

void add_force(double r, int m, double complex *res) {

    force(r,m,&res[1],&res[0]);

    return;
}
void lw_inner_bc(double r0, int m, int eps, double complex *md0, double complex *ud0, double complex *ld0) {
    main_diag(r0,m,md0);
    upper_diag(r0,m,ud0);
    lower_diag(r0,m,ld0);
    
    double rb = r0 * exp(-.5*params.dlr);
    double kr = rb*pow( fabs(Dfunc(rb,params.omf,m)/cs2(rb)) ,.5);
    
    double fac = (1 - .5*eps*I*kr*params.dlr)/(1 + .5*eps*I*kr*params.dlr);

    int i;
    for(i=0;i<9;i++) {
        md0[i] += fac*ld0[i];
    }
    return;
} 
void lw_outer_bc(double rn, int m, int eps, double complex *mdn, double complex *udn, double complex *ldn) {
    main_diag(rn,m,mdn);
    upper_diag(rn,m,udn);
    lower_diag(rn,m,ldn);
    
    double rb = rn * exp(.5*params.dlr);
    double kr = rb*pow( fabs(Dfunc(rb,params.omf,m)/cs2(rb)) ,.5);
    
    double fac = (1 + .5*eps*I*kr*params.dlr)/(1 - .5*eps*I*kr*params.dlr);

    int i;
    for(i=0;i<9;i++) {
        mdn[i] += fac*udn[i];
    }
    return;
} 

void construct_matrix(double *r, double complex *ld, double complex *md, double complex *ud, double complex *fd, int m) {
    int i;
    int size = params.nrhs*params.nrhs;
    int n = params.n;


    lw_inner_bc(r[0], m, 1, &md[0], &ud[0], &fd[0]);
    add_force(r[0], m, &fd[0]);

    for(i=1;i<n-1;i++) {
        main_diag(r[i],m, &md[i*size]);
        upper_diag(r[i],m, &ud[i*size]);
        lower_diag(r[i],m, &ld[(i-1)*size]);
        add_force(r[i], m, &fd[i*params.nrhs]);
    }

    i = n-1;
    lw_outer_bc(r[i], m, 1, &md[i*size], &fd[i*params.nrhs], &ld[(i-1)*size]);
    add_force(r[i], m, &fd[i*params.nrhs]);

    return;
}

