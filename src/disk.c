#include "linearwaves.h"


double scaleH(double x, Params params, Disk *disk) {
    return params.h * x * pow(x,params.f);
}
double omegaK(double x, Params params, Disk *disk) {
    return pow(x,-1.5);
}
double nu_func(double x, Params params, Disk *disk) {
    return params.alpha * params.h*params.h * pow(x, params.nuindx);
}
double cs2_func(double x, Params params, Disk *disk) {
    return params.h*params.h * pow(x,params.delta);
}
double cs_func(double x, Params params, Disk *disk) {
    return pow(cs2_func(x, params,disk),.5);
}

double vr_func(double x, Params params, Disk *disk) {
    return -1.5*nu_func(x,params,disk)/x;
}
double dlvr_func(double x, Params params, Disk *disk) {
    return vr_func(x,params,disk) * (params.nuindx - 1);
}

double omega2_func(double x, Params params, Disk *disk) {
    if (!params.pcorrect) {
        return pow(x,-3);
    }
    else {
        if (params.iso) {
            return pow(x,-3) + cs2_func(x, params,disk)/(x*x) * (params.delta + params.mu);
        }
        else {
            return pow(x,-3) + cs2_func(x, params,disk)/(x*x) * params.mu;
        }
    }
}
double kappa2_func(double x, Params params, Disk *disk) {
    if (!params.pcorrect) {
        return pow(x,-3);
    }
    else {
        if (params.iso) {
            return pow(x,-3) + cs2_func(x,params,disk)/(x*x)*( (2+params.delta)*(params.delta+params.mu));
        }
        else {
            return pow(x,-3) + cs2_func(x,params,disk)/(x*x)*params.mu*(2+params.delta);
        }
    }
}
double kappa_func(double x, Params params, Disk *disk) {
    return pow(kappa2_func(x, params,disk),.5);
}
double omega_func(double x, Params params, Disk *disk) {
    return pow(omega2_func(x, params,disk),.5);
}
double k2om_func(double x, Params params, Disk *disk) {
    return kappa2_func(x, params,disk)/(2*omega_func(x, params,disk));
}
double sigma_func(double x, Params params, Disk *disk) {
    return params.sig0 * pow(x,params.mu);
}
double pres_func(double x, Params params, Disk *disk) {
    return cs2_func(x,params,disk)*sigma_func(x, params,disk);
}
double dsdr_func(double x, Params params, Disk *disk) {
    return params.mu * sigma_func(x, params,disk)/x;
}
double dpdr_func(double x, Params params, Disk *disk) {
    if (params.iso) {
        return (params.delta + params.mu)*pres_func(x, params,disk)/x;
    }
    else {
        return cs2_func(x,params,disk)*dsdr_func(x, params,disk);
    }
}
double dc2dr_func(double x, Params params, Disk *disk) {
    return params.delta * cs2_func(x,params,disk)/x;
}
double Dfunc(int i, double omp, int m, Params params, Disk *disk) {
    return disk->kappa2[i] - m*m*(disk->omega[i] - omp)*(disk->omega[i]-omp);
}
double complex sig(int i, double complex s, Params params, Disk *disk) {
    if (params.iso) {
        return s*disk->sigma[i];

    }
    else {
        return s*disk->sigma[i]/disk->c2[i];
    }

}

void alloc_disk(Disk *disk,int n){

    disk->c2 = (double *)malloc(sizeof(double)*n);
    disk->sigma = (double *)malloc(sizeof(double)*n);
    disk->omega = (double *)malloc(sizeof(double)*n);
    disk->dlsdlr = (double *)malloc(sizeof(double)*n);
    disk->d2lsdlr = (double *)malloc(sizeof(double)*n);
    disk->dlomdlr = (double *)malloc(sizeof(double)*n);
    disk->d2lomdlr = (double *)malloc(sizeof(double)*n);
    disk->dlnudlr = (double *)malloc(sizeof(double)*n);
    disk->dlTdlr = (double *)malloc(sizeof(double)*n);
    disk->d2lTdlr = (double *)malloc(sizeof(double)*n);
    disk->kappa2 = (double *)malloc(sizeof(double)*n);
    disk->nu = (double *)malloc(sizeof(double)*n);
    disk->pres = (double *)malloc(sizeof(double)*n);
    disk->dpdr = (double *)malloc(sizeof(double)*n);
    disk->vrbar = (double *)malloc(sizeof(double)*n);
    disk->dlvrbar = (double *)malloc(sizeof(double)*n);
    return;

}
void free_disk(Disk *disk ){

    
    SAFE_FREE(disk->c2);
    SAFE_FREE(disk->sigma);
    SAFE_FREE(disk->omega);
    SAFE_FREE(disk->dlsdlr);
    SAFE_FREE(disk->d2lsdlr);
    SAFE_FREE(disk->dlomdlr);
    SAFE_FREE(disk->d2lomdlr);
    SAFE_FREE(disk->dlnudlr);
    SAFE_FREE(disk->dlTdlr);
    SAFE_FREE(disk->d2lTdlr);
    SAFE_FREE(disk->kappa2);
    SAFE_FREE(disk->nu);
    SAFE_FREE(disk->pres);
    SAFE_FREE(disk->dpdr);
    SAFE_FREE(disk->dlvrbar);
    SAFE_FREE(disk->vrbar);
    return;
}

void init_disk(char *fname, double *lr,Disk *disk,Params params) {
    alloc_disk(disk,params.n);
    if (params.fromfile) {
        read_sigma(fname, lr, disk->sigma,disk->dlsdlr,disk->d2lsdlr,disk->omega,disk->dlomdlr, disk->d2lomdlr, params.n, params.readomega);
    }
    int i;
    double x,om2,k2; 
    for(i=0;i<params.n;i++) {
        x = exp(lr[i]);
        disk->c2[i] = cs2_func(x, params,disk); 
        disk->nu[i] = nu_func(x, params,disk);
        disk->dlnudlr[i] = params.nuindx;
        disk->dlTdlr[i] = params.delta;
        disk->d2lTdlr[i] = 0;
        disk->vrbar[i] = 0;
        disk->dlvrbar[i] = 0;
        if (params.fromfile) {
            disk->sigma[i] = exp(disk->sigma[i]);
            disk->d2lsdlr[i] = (disk->d2lsdlr[i]<0) ? fmax(disk->d2lsdlr[i],-1./(params.h*params.h)) : fmin(disk->d2lsdlr[i],1./(params.h*params.h));

        }
        else {
            disk->sigma[i] = sigma_func(x, params,disk);
            disk->dlsdlr[i] = params.mu;
            disk->d2lsdlr[i] = 0;
            disk->vrbar[i] = vr_func(x,params,disk);
            disk->vrbar[i] = dlvr_func(x,params,disk);
        }
        if ((params.readomega) && (params.fromfile)) {
                disk->omega[i] = exp(disk->omega[i]);
                disk->kappa2[i] = (2 + disk->dlomdlr[i])*2*disk->omega[i]*disk->omega[i];
        }
        else {
            if (!params.pcorrect) {
                disk->omega[i] = omegaK(x, params,disk);
                disk->kappa2[i] = disk->omega[i]*disk->omega[i];
                disk->dlomdlr[i] = -1.5;
            }
            else {
                if (params.iso) {
                    om2 = pow(x,-3);
                    om2 += disk->c2[i]/(x*x) *(disk->dlTdlr[i] + disk->dlsdlr[i]) ;
                    disk->omega[i] = pow(om2,.5);
                    k2 = pow(x,-3);
                    k2 += disk->c2[i]/(x*x) *( (2 + disk->dlTdlr[i])*(disk->dlTdlr[i] + disk->dlsdlr[i]) + disk->d2lTdlr[i] + disk->d2lsdlr[i]) ;
                    disk->kappa2[i] = k2;
                }
                else {
                    om2 = pow(x,-3);
                    om2 += disk->c2[i]/(x*x) *(disk->dlsdlr[i]) ;
                    disk->omega[i] = pow(om2,.5);
                    k2 = pow(x,-3);
                    k2 += disk->c2[i]/(x*x) *( (2 + disk->dlTdlr[i])*(disk->dlsdlr[i]) + disk->d2lsdlr[i]) ;
                    disk->kappa2[i] = k2;
                }
                disk->dlomdlr[i] = disk->kappa2[i]/(2*disk->omega[i]*disk->omega[i]) - 2;
                disk->d2lomdlr[i] = 0;
            }
        }

        disk->pres[i] = disk->c2[i]*disk->sigma[i];
        if (params.iso) {
            disk->dpdr[i] = disk->pres[i]/x *(disk->dlTdlr[i] + disk->dlsdlr[i]); 
        }
        else {
            disk->dpdr[i] = disk->dlsdlr[i] * disk->pres[i]/x;
        }

    }

    FILE *f = fopen("omega.dat","w");
    fwrite(disk->omega,sizeof(double),params.n,f);
    fwrite(disk->dlomdlr,sizeof(double),params.n,f);
    fwrite(disk->d2lomdlr,sizeof(double),params.n,f);
    fclose(f);

}

