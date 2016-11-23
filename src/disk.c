#include "linearwaves.h"

void alloc_disk(void ){

    disk.c2 = (double *)malloc(sizeof(double)*params.n);
    disk.sigma = (double *)malloc(sizeof(double)*params.n);
    disk.omega = (double *)malloc(sizeof(double)*params.n);
    disk.dlsdlr = (double *)malloc(sizeof(double)*params.n);
    disk.d2lsdlr = (double *)malloc(sizeof(double)*params.n);
    disk.dlomdlr = (double *)malloc(sizeof(double)*params.n);
    disk.dlnudlr = (double *)malloc(sizeof(double)*params.n);
    disk.dlTdlr = (double *)malloc(sizeof(double)*params.n);
    disk.d2lTdlr = (double *)malloc(sizeof(double)*params.n);
    disk.kappa2 = (double *)malloc(sizeof(double)*params.n);
    disk.nu = (double *)malloc(sizeof(double)*params.n);
    disk.pres = (double *)malloc(sizeof(double)*params.n);
    disk.dpdr = (double *)malloc(sizeof(double)*params.n);
    return;

}
void free_disk(void ){

    
    free(disk.c2);
    free(disk.sigma);
    free(disk.omega);
    free(disk.dlsdlr);
    free(disk.d2lsdlr);
    free(disk.dlomdlr);
    free(disk.dlnudlr);
    free(disk.dlTdlr);
    free(disk.d2lTdlr);
    free(disk.kappa2);
    free(disk.nu);
    free(disk.pres);
    free(disk.dpdr);
    return;
}

void init_disk(char *fname, double *lr) {
    alloc_disk();
    if (params.fromfile) {
        read_sigma(fname, lr, disk.sigma,disk.dlsdlr,disk.d2lsdlr,params.n);
    }
    int i;
    double x,om2,k2; 
    for(i=0;i<params.n;i++) {
        x = exp(lr[i]);
        disk.c2[i] = cs2_func(x); 
        disk.nu[i] = nu_func(x);
        disk.dlnudlr[i] = params.nuindx;
        disk.dlTdlr[i] = params.delta;
        disk.d2lTdlr[i] = 0;
        if (params.fromfile) {
        disk.sigma[i] = exp(disk.sigma[i]);
        disk.d2lsdlr[i] = (disk.d2lsdlr[i]<0) ? fmax(disk.d2lsdlr[i],-1./(params.h*params.h)) : fmin(disk.d2lsdlr[i],1./(params.h*params.h));

        }
        else {
            disk.sigma[i] = sigma_func(x);
            disk.dlsdlr[i] = params.mu;
            disk.d2lsdlr[i] = 0;
        }
        if (!params.pcorrect) {
            disk.omega[i] = omegaK(x);
            disk.kappa2[i] = disk.omega[i]*disk.omega[i];
            disk.dlomdlr[i] = -1.5;
        }
        else {
            if (params.iso) {
                om2 = pow(x,-3);
                om2 += disk.c2[i]/(x*x) *(disk.dlTdlr[i] + disk.dlsdlr[i]) ;
                disk.omega[i] = pow(om2,.5);
                k2 = pow(x,-3);
                k2 += disk.c2[i]/(x*x) *( (2 + disk.dlTdlr[i])*(disk.dlTdlr[i] + disk.dlsdlr[i]) + disk.d2lTdlr[i] + disk.d2lsdlr[i]) ;
                disk.kappa2[i] = k2;
            }
            else {
                om2 = pow(x,-3);
                om2 += disk.c2[i]/(x*x) *(disk.dlsdlr[i]) ;
                disk.omega[i] = pow(om2,.5);
                k2 = pow(x,-3);
                k2 += disk.c2[i]/(x*x) *( (2 + disk.dlTdlr[i])*(disk.dlsdlr[i]) + disk.d2lsdlr[i]) ;
                disk.kappa2[i] = k2;
            }
            disk.dlomdlr[i] = disk.kappa2[i]/(2*disk.omega[i]*disk.omega[i]) - 2;
        }

        disk.pres[i] = disk.c2[i]*disk.sigma[i];
        if (params.iso) {
            disk.dpdr[i] = disk.pres[i]/x *(disk.dlTdlr[i] + disk.dlsdlr[i]); 
        }
        else {
            disk.dpdr[i] = disk.dlsdlr[i] * disk.pres[i]/x;
        }

    }



}



double scaleH(double x) {
    return params.h * x * pow(x,params.f);
}
double omegaK(double x) {
    return pow(x,-1.5);
}
double nu_func(double x) {
    return params.alpha * params.h*params.h * pow(x, params.nuindx);
}
double cs2_func(double x) {
    return params.h*params.h * pow(x,params.delta);
}
double cs_func(double x) {
    return pow(cs2_func(x),.5);
}

double omega2_func(double x) {
    if (!params.pcorrect) {
        return pow(x,-3);
    }
    else {
        if (params.iso) {
            return pow(x,-3) + cs2_func(x)/(x*x) * (params.delta + params.mu);
        }
        else {
            return pow(x,-3) + cs2_func(x)/(x*x) * params.mu;
        }
    }
}
double kappa2_func(double x) {
    if (!params.pcorrect) {
        return pow(x,-3);
    }
    else {
        if (params.iso) {
            return pow(x,-3) + cs2_func(x)/(x*x)*( (2+params.delta)*(params.delta+params.mu));
        }
        else {
            return pow(x,-3) + cs2_func(x)/(x*x)*params.mu*(2+params.delta);
        }
    }
}
double kappa_func(double x) {
    return pow(kappa2_func(x),.5);
}
double omega_func(double x) {
    return pow(omega2_func(x),.5);
}
double k2om_func(double x) {
    return kappa2_func(x)/(2*omega_func(x));
}
double sigma_func(double x) {
    return params.sig0 * pow(x,params.mu);
}
double pres_func(double x) {
    return cs2_func(x)*sigma_func(x);
}
double dsdr_func(double x) {
    return params.mu * sigma_func(x)/x;
}
double dpdr_func(double x) {
    if (params.iso) {
        return (params.delta + params.mu)*pres_func(x)/x;
    }
    else {
        return cs2_func(x)*dsdr_func(x);
    }
}
double dc2dr_func(double x) {
    return params.delta * cs2_func(x)/x;
}
double Dfunc(int i, double x, double omp, int m) {
    return disk.kappa2[i] - m*m*(disk.omega[i] - omp)*(disk.omega[i]-omp);
}
double complex sig(int i, double x, double complex s) {
    if (params.iso) {
        return s*disk.sigma[i];

    }
    else {
        return s*disk.sigma[i]/disk.c2[i];
    }

}
