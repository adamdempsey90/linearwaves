#include "linearwaves.h"


void main_diag(int indx, double r, int m, double complex *res, Params params, Disk disk) {
    double om = disk.omega[indx];
    double c2 = disk.c2[indx];
    double mu = disk.dlsdlr[indx];
    double k2om = disk.kappa2[indx]/(2*om);
    double phat = m*(om - params.omf);
    double invdlr2 = -2./(params.dlr*params.dlr);
    double norm = (params.iso) ? 1.0 : c2;
    int i;
   for(i=0;i<9;i++) res[i] = 0;

    res[0] = I*(phat - params.ieps*I); 
    res[1] = -2*om;                     
    res[2] = 0.0;                     
    res[3] = k2om;                     
    res[4] = I*(phat -params.ieps*I);     
    res[5] = I*m*c2/(norm*r);                     
    res[6] = (mu + 1)*norm/r ;                     
    res[7] = I*m*norm/r  ;                     
    res[8] = I*(phat - params.ieps*I)  ;    

    viscosity_coeffs_u(indx,r,res,m,1.0);
    viscosity_d2coeffs_u(indx,r,res,invdlr2);
    viscosity_coeffs_v(indx,r,&res[params.nrhs],m,1.0);
    viscosity_d2coeffs_v(indx,r,&res[params.nrhs],invdlr2);
    return;
}

void upper_diag(int indx, double r, int m, double complex *res, Params params, Disk disk) {

   double c2 = disk.c2[indx];
   double invdlr = .5/params.dlr;
   double invdlr2 = 1./(params.dlr*params.dlr);
   double norm = (params.iso) ? 1.0 : c2;
   int i;
   for(i=0;i<9;i++) res[i] = 0;

   res[2] = c2/(norm*r) * invdlr;
   res[6] = norm/ r * invdlr;

   viscosity_dcoeffs_u(indx,r,res,m,invdlr);
   viscosity_dcoeffs_v(indx,r,&res[params.nrhs],m,invdlr);
    viscosity_d2coeffs_u(indx,r,res,invdlr2);
    viscosity_d2coeffs_v(indx,r,&res[params.nrhs],invdlr2);
    return;
}


void lower_diag(int indx, double r, int m, double complex *res, Params params, Disk disk) {
    double c2 = disk.c2[indx];
    double invdlr = -.5/params.dlr;
    double invdlr2 = 1./(params.dlr*params.dlr);
    double norm = (params.iso) ? 1.0 : c2;
    int i;
    for(i=0;i<9;i++) res[i] = 0;
    res[2] = c2/(norm*r) * invdlr;
    res[6] = norm/ r * invdlr;

   viscosity_dcoeffs_u(indx,r,res,m,invdlr);
   viscosity_dcoeffs_v(indx,r,&res[params.nrhs],m,invdlr);
    viscosity_d2coeffs_u(indx,r,res,invdlr2);
    viscosity_d2coeffs_v(indx,r,&res[params.nrhs],invdlr2);


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

void lw_inner_bc(int indx, double r0, int m, int eps, double complex *md0, double complex *ud0, Params params, Disk disk) {
    main_diag(indx,r0,m,md0,params,disk);
    lower_diag(indx,r0,m,ud0,params,disk);
        
    double c2 = disk.c2[indx];
    double rb = r0 * exp(-.5*params.dlr);
    double kr = rb*pow( fabs(Dfunc(indx,params.omf,m,params,disk)/c2) ,.5);
    
    double complex fac = (1 - .5*eps*I*kr*params.dlr)/(1 + .5*eps*I*kr*params.dlr);

//    printf("%lg\t%lg\t%lg+I%lg\n",rb,kr,creal(fac),cimag(fac));
    int i;
    for(i=0;i<9;i++) {
        md0[i] += fac*ud0[i];
    }
    upper_diag(indx,r0,m,ud0,params,disk);
    return;
} 
void lw_outer_bc(int indx, double rn, int m, int eps, double complex *mdn,  double complex *ldn, Params params, Disk disk) {
    main_diag(indx,rn,m,mdn,params,disk);
    upper_diag(indx,rn,m,ldn,params,disk);
    
    double c2 = disk.c2[indx];
    double rb = rn * exp(.5*params.dlr);
    double kr = rb*pow( fabs(Dfunc(indx,params.omf,m,params,disk)/c2) ,.5);
    
    double complex fac = (1 + .5*eps*I*kr*params.dlr)/(1 - .5*eps*I*kr*params.dlr);

    int i;
    for(i=0;i<9;i++) {
        mdn[i] += fac*ldn[i];
    }
    lower_diag(indx,rn,m,ldn,params,disk);
    return;
} 

void construct_matrix(double *r, double complex *ld, double complex *md, double complex *ud, double complex *fd, double *dppot, double *drpot, int m, Params params, Disk disk) {
    int i;
    int size = params.nrhs*params.nrhs;
    int n = params.n;

    if ((params.zero_inner_bc) || (m==1)) {
        zero_inner_bc(&md[0], &ud[0]);
        fd[0] = 0; fd[1] = 0; fd[2] = 0;
    }
    else {
        lw_inner_bc(0,r[0], m, 1, &md[0], &ud[0],params,disk);
        fd[0] = 0; fd[1] = 0; fd[2] = 0;
    }
    for(i=1;i<n-1;i++) {
        main_diag(i,r[i],m, &md[i*size],params,disk);
        upper_diag(i,r[i],m, &ud[i*size],params,disk);
        lower_diag(i,r[i],m, &ld[(i-1)*size],params,disk);
        fd[i*params.nrhs] = -drpot[i]; 
        fd[i*params.nrhs+1] = -dppot[i]*I*m/r[i];
        fd[i*params.nrhs+2] = 0;
    }

    i = n-1;
    if (params.zero_outer_bc) {
        zero_outer_bc(&md[i*size],  &ld[(i-1)*size]);
        fd[i*params.nrhs] = 0; fd[i*params.nrhs+1] = 0; fd[i*params.nrhs+1] = 0;
    }
    else {
        lw_outer_bc(i,r[i], m, 1, &md[i*size],  &ld[(i-1)*size],params,disk);
        fd[i*params.nrhs] = 0; fd[i*params.nrhs+1] = 0; fd[i*params.nrhs+1] = 0;
    }
    return;
}

