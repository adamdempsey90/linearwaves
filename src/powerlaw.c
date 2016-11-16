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
double cs2(double x) {
    return params.h*params.h * pow(x,params.delta);
}
double cs(double x) {
    return pow(cs2(x),.5);
}

double omega2(double x) {
    if (!params.pcorrect) {
        return pow(x,-3);
    }
    else {
        if (params.iso) {
            return pow(x,-3) + cs2(x)/(x*x) * (params.delta + params.mu);
        }
        else {
            return pow(x,-3) + cs2(x)/(x*x) * params.mu;
        }
    }
}
double kappa2(double x) {
    if (!params.pcorrect) {
        return pow(x,-3);
    }
    else {
        if (params.iso) {
            return pow(x,-3) + cs2(x)/(x*x)*( (2+params.delta)*(params.delta+params.mu));
        }
        else {
            return pow(x,-3) + cs2(x)/(x*x)*params.mu*(2+params.delta);
        }
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
double complex sig(double x, double complex s) {
    if (params.iso) {
        return s*sigma(x);

    }
    else {
        return s*sigma(x)/cs2(x);
    }

}
