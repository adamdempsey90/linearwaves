#include "linearwaves.h"

void viscosity_coeffs_u(int indx, double r, double complex *res, int m) {
    
    double visc = -disk.nu[indx];
    double eta = params.eta;
    double om = disk.omega[indx];
    double gamma = disk.dlnudlr[indx] + disk.dlsdlr[indx];
    double c2 = disk.c2[indx];
    double norm = (params.iso) ? 1.0 : c2;
    double drom = disk.kappa2[indx]/(2*om) - 2*om;
    double invr = 1./r;
    double invr2 = invr*invr;

    if (params.simple_visc) {
        res[0] += -visc*(m*m + 1)*invr2;
        res[1] += -visc*2*I*m*invr2;
    }
    else {
        res[0] += -visc*(m*m + gamma*eta)*invr2;
        res[1] += -visc*(3 + eta*(gamma-1))*I*m*invr2;
        res[2] += visc*I*m*drom*norm*invr;
    }

    return;
}
void viscosity_coeffs_v(int indx, double r, double complex *res, int m) {
    
    double visc = -disk.nu[indx];
    double eta = params.eta;
    double om = disk.omega[indx];
    double gamma = disk.dlnudlr[indx] + disk.dlsdlr[indx];
    double c2 = disk.c2[indx];
    double norm = (params.iso) ? 1.0 : c2;
    double drom = disk.kappa2[indx]/(2*om) - 2*om;
    double invr = 1./r;
    double invr2 = invr*invr;

    if (params.simple_visc) {
        res[0] += visc*2*I*m*invr2;
        res[1] += -visc*(1 + m*m)*invr2;
    }
    else {
        res[0] += visc*I*m*(3 - eta + gamma)*invr2;
        res[1] += -visc*(1+m*m*(2-eta)+gamma)*invr2;
        if (!params.iso) {
            res[2] += -visc*drom*params.delta*invr;
        }
    }

    return;
}
void viscosity_dcoeffs_u(int indx, double r, double complex *res, int m, double invdlr) {
    double visc = -disk.nu[indx];
    double eta = params.eta;
    double om = disk.omega[indx];
    double gamma = disk.dlnudlr[indx] + disk.dlsdlr[indx];
    double c2 = disk.c2[indx];
    double norm = (params.iso) ? 1.0 : c2;
    double drom = disk.kappa2[indx]/(2*om) - 2*om;

    double invr = 1./r;
    double invr2 = invr*invr;

    if (!params.simple_visc) {
        res[0]  += visc*(2*(gamma-1)-eta*gamma)*invr2*invdlr;
        res[1] += visc*I*m*invr2*invdlr;
    }
    return;
}
void viscosity_dcoeffs_v(int indx, double r, double complex *res, int m, double invdlr) {
    
    double visc = -disk.nu[indx];
    double eta = params.eta;
    double om = disk.omega[indx];
    double gamma = disk.dlnudlr[indx] + disk.dlsdlr[indx];
    double c2 = disk.c2[indx];
    double norm = (params.iso) ? 1.0 : c2;
    double drom = disk.kappa2[indx]/(2*om) - 2*om;

    double invr = 1./r;
    double invr2 = invr*invr;

    if (!params.simple_visc) {
        res[0] += visc*I*m*(1-eta)*invr2*invdlr;
        res[1] += visc*gamma*invr2*invdlr; 
        res[2] += visc*drom*invr*invdlr/norm;
    }

    return;
}
void viscosity_d2coeffs_u(int indx, double r, double complex *res,int m, double invdlr2) {
    
    double visc = -disk.nu[indx];
    double eta = params.eta;
    double invr2 = 1./(r*r);

    if (params.simple_visc) {
        res[0] += visc*invr2*invdlr2;
    }
    else {
        res[0] += visc*(2-eta)*invr2*invdlr2;
    }
    return;
}
void viscosity_d2coeffs_v(int indx, double r, double complex *res,int m, double invdlr2) {
    
    double visc = -disk.nu[indx];
    double eta = params.eta;
    double invr2 = 1./(r*r);

    if (params.simple_visc) {
        res[1] += visc*invr2*invdlr2;
    }
    else {
        res[1] += visc*invr2*invdlr2;
    } 
    return;
}
