#include "linearwaves.h"

double scaleH(double x) {
    return params.h * x * pow(x,params.f);
}
double omegaK(double x) {
    return pow(x,-1.5);
}
double nu(double x) {
    return params.alpha * params.h*params.h * pow(x, params.nuindx);
}
double cs(double x) {
    return scaleH(x)*omegaK(x);
}
double cs2(double x) {
    return cs(x)*cs(x);
}

double omega2(double x) {
    if (params.iso) {
        return pow(x,-3) + cs2(x)/(x*x) * (params.delta + params.mu);
    }
    else {
        return pow(x,-3) + cs2(x)/(x*x) * params.mu;
    }
}
double kappa2(double x) {

    if (params.iso) {
        return pow(x,-3) + cs2(x)/(x*x)*( (2+params.delta)*(params.delta+params.mu));
    }
    else {
        return pow(x,-3) + cs2(x)/(x*x)*params.mu*(2+params.delta);
    }
}
double kappa(double x) {
    return pow(kappa2(x),.5);
}
double omega(double x) {
    return pow(omega2(x),.5);
}
double k2om(double x) {
    return kappa2(x)/(2*omega(x));
}
double sigma(double x) {
    return params.sig0 * pow(x,params.mu);
}
double pres(double x) {
    return cs2(x)*sigma(x);
}
double dsdr(double x) {
    return params.mu * sigma(x)/x;
}
double dpdr(double x) {
    if (params.iso) {
        return (params.delta + params.mu)*pres(x)/x;
    }
    else {
        return cs2(x)*dsdr(x);
    }
}
double dc2dr(double x) {
    return params.delta * cs2(x)/x;
}
double Dfunc(double x, double omp, int m) {
    return kappa2(x) - m*m*(omega(x) - omp)*(omega(x)-omp);
}
