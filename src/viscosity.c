#include "linearwaves.h"

void viscosity_coeffs_u(double r, double complex *res, int m) {
    
    double visc = - nu(r);
    double eta = params.eta;
    double gamma = params.nuindx + params.mu;
    double c2 = cs2(r);
    double norm = (params.iso) ? 1.0 : c2;
    double drom = kappa2(r)/(2*omega(r)) - 2*omega(r);

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
void viscosity_coeffs_v(double r, double complex *res, int m) {
    
    double visc = - nu(r);
    double eta = params.eta;
    double gamma = params.nuindx + params.mu;
    double c2 = cs2(r);
    double norm = (params.iso) ? 1.0 : c2;
    double drom = kappa2(r)/(2*omega(r)) - 2*omega(r);

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
void viscosity_dcoeffs_u(double r, double complex *res, int m, double invdlr) {
    
    double visc = -nu(r);
    double eta = params.eta;
    double gamma = params.nuindx + params.mu;
    double c2 = cs2(r);
    double norm = (params.iso) ? 1.0 : c2;
    double drom = kappa2(r)/(2*omega(r)) - 2*omega(r);

    double invr = 1./r;
    double invr2 = invr*invr;

    if (!params.simple_visc) {
        res[0]  += visc*(2*(gamma-1)-eta*gamma)*invr2*invdlr;
        res[1] += visc*I*m*invr2*invdlr;
    }
    return;
}
void viscosity_dcoeffs_v(double r, double complex *res, int m, double invdlr) {
    
    double visc = -nu(r);
    double eta = params.eta;
    double gamma = params.nuindx + params.mu;
    double c2 = cs2(r);
    double norm = (params.iso) ? 1.0 : c2;
    double drom = kappa2(r)/(2*omega(r)) - 2*omega(r);

    double invr = 1./r;
    double invr2 = invr*invr;

    if (!params.simple_visc) {
        res[0] += visc*I*m*(1-eta)*invr2*invdlr;
        res[1] += visc*gamma*invr2*invdlr; 
        res[2] += visc*drom*invr*invdlr/norm;
    }

    return;
}
void viscosity_d2coeffs_u(double r, double complex *res,int m, double invdlr2) {
    
    double visc = -nu(r);
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
void viscosity_d2coeffs_v(double r, double complex *res,int m, double invdlr2) {
    
    double visc = -nu(r);
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
