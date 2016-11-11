#include "linearwaves.h"

#define INDX(i,j) j + 3*i

void main_diag(double r, int m, double complex *res) {
    double om = omega(r);
    double c2 = cs2(r);
    double k2om = kappa2(r)/(2*om);
    double invdlr2 = -2./(params.dlr*params.dlr);
    double norm = (params.iso) ? 1.0 : c2;
    int i,j;
   for(i=0;i<9;i++) res[i] = 0;

    res[0] = I*(m*(om - params.omf) - params.ieps*I);     
    res[1] = -2*om;                     
    res[2] = 0.0;                     
    res[3] = k2om;                     
    res[4] = I*(m*(om - params.omf) -params.ieps*I);     
    res[5] = I*m*c2/(norm*r);                     
    res[6] = (params.mu + 1)*norm/r;                     
    res[7] = I*m*norm/r;                     
    res[8] = I*(m*(om - params.omf)-params.ieps*I);    
/*
    printf("%lg\n",r);
    for(i=0;i<3;i++) {
        for(j=0;j<3;j++) {
            printf("%.03f+%.03fI\t",creal(res[INDX(i,j)]),cimag(res[INDX(i,j)]));
        }
        printf("\n");
    }
    */
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
   double norm = (params.iso) ? 1.0 : c2;
   int i;
   for(i=0;i<9;i++) res[i] = 0;
   res[2] = c2/(norm*r) * invdlr;
   res[6] = norm*invdlr / r;
    /*
    printf("%lg\n",r);
    int j;
    for(i=0;i<3;i++) {
        for(j=0;j<3;j++) {
            printf("%.03f+%.03fI\t",creal(res[INDX(i,j)]),cimag(res[INDX(i,j)]));
        }
        printf("\n");
    }
    */
    double visc = -nu(r);
    double ir2 = 1./(r*r);
    res[0] += visc*(2*params.nuindx+1)*ir2*invdlr;
    res[1] += visc*I*m*(1+params.eta)*ir2*invdlr;
    res[3] += visc*(1-params.eta)*I*m*ir2*invdlr;
    res[4] += visc*params.nuindx*ir2*invdlr;
    res[0] += visc*(2+params.eta)*ir2*invdlr2;
    res[4] += visc*ir2*invdlr2;
    res[5] += visc*I*m*(k2om-2*om)/(r*norm)*invdlr;
    return;
}

void lower_diag(double r, int m, double complex *res) {
    double om = omega(r);
    double k2om = kappa2(r)/(2*om);
    double c2 = cs2(r);
    double invdlr = -.5/params.dlr;
    double invdlr2 = 1./(params.dlr*params.dlr);
    double norm = (params.iso) ? 1.0 : c2;
    int i;
    for(i=0;i<9;i++) res[i] = 0;
    res[2] = c2/(norm*r) * invdlr;
    res[6] = norm*invdlr / r;
    double visc = -nu(r);
    double ir2 = 1./(r*r);

    res[0] += visc*(2*params.nuindx+1)*ir2*invdlr;
    res[1] += visc*I*m*(1+params.eta)*ir2*invdlr;
    res[3] += visc*(1-params.eta)*I*m*ir2*invdlr;
    res[4] += visc*params.nuindx*ir2*invdlr;
    res[0] += visc*(2+params.eta)*ir2*invdlr2;
    res[4] += visc*ir2*invdlr2;
    res[5] += visc*I*m*(k2om-2*om)/(r*norm) *invdlr;


    return;
}

void add_force(double r, int m, double complex *res) {

    force(r,m,res);

    return;
}

void zero_inner_bc(double complex *md0, double complex *ud0) {

    int i;
    for(i=0;i<9;i++) {
        ud0[i] = 0;
        md0[i] = 0;
    }
    md0[0] = 1;
    md0[4] = 1;
    md0[8] = 1;

    return;

}
void zero_outer_bc(double complex *mdn, double complex *ldn) {

    int i;
    for(i=0;i<9;i++) {
        ldn[i] = 0;
        mdn[i] = 0;
    }

    mdn[0] = 1;
    mdn[4] = 1;
    mdn[8] = 1;
    return;

}

void lw_inner_bc(double r0, int m, int eps, double complex *md0, double complex *ud0) {
    main_diag(r0,m,md0);
    lower_diag(r0,m,ud0);
    
    double rb = r0 * exp(-.5*params.dlr);
    double kr = rb*pow( fabs(Dfunc(rb,params.omf,m)/cs2(rb)) ,.5);
    
    double complex fac = (1 - .5*eps*I*kr*params.dlr)/(1 + .5*eps*I*kr*params.dlr);

//    printf("%lg\t%lg\t%lg+I%lg\n",rb,kr,creal(fac),cimag(fac));
    int i;
    for(i=0;i<9;i++) {
        md0[i] += fac*ud0[i];
    }
    upper_diag(r0,m,ud0);
    return;
} 
void lw_outer_bc(double rn, int m, int eps, double complex *mdn,  double complex *ldn) {
    main_diag(rn,m,mdn);
    upper_diag(rn,m,ldn);
    
    double rb = rn * exp(.5*params.dlr);
    double kr = rb*pow( fabs(Dfunc(rb,params.omf,m)/cs2(rb)) ,.5);
    
    double complex fac = (1 + .5*eps*I*kr*params.dlr)/(1 - .5*eps*I*kr*params.dlr);

    int i;
    for(i=0;i<9;i++) {
        mdn[i] += fac*ldn[i];
    }
    lower_diag(rn,m,ldn);
    return;
} 

void construct_matrix(double *r, double complex *ld, double complex *md, double complex *ud, double complex *fd, int m) {
    int i;
    int size = params.nrhs*params.nrhs;
    int n = params.n;

    if (params.zero_inner_bc) {
        zero_inner_bc(&md[0], &ud[0]);
        fd[0] = 0; fd[1] = 0; fd[2] = 0;
    }
    else {
        printf("Outgoing inner boundary\n");
        lw_inner_bc(r[0], m, 1, &md[0], &ud[0]);
        add_force(r[0], m, &fd[0]);
    }
#ifdef _OPENMP
#pragma omp parallel for private(i)
#endif
    for(i=1;i<n-1;i++) {
        main_diag(r[i],m, &md[i*size]);
        upper_diag(r[i],m, &ud[i*size]);
        lower_diag(r[i],m, &ld[(i-1)*size]);
        add_force(r[i], m, &fd[i*params.nrhs]);
    }

    i = n-1;
    if (params.zero_outer_bc) {
        zero_outer_bc(&md[i*size],  &ld[(i-1)*size]);
        fd[i*params.nrhs] = 0; fd[i*params.nrhs+1] = 0; fd[i*params.nrhs+1] = 0;
    }
    else {
        printf("Outgoing outer boundary\n");
        lw_outer_bc(r[i], m, 1, &md[i*size],  &ld[(i-1)*size]);
        add_force(r[i], m, &fd[i*params.nrhs]);
    }
    return;
}

