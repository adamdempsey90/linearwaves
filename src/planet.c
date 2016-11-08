

typedef struct GSL_params {
    double x;
    int m;
} GSL_params;

double potential(double phi, double x) {
    double res = -planet.mp * pow(x*x + planet.a*planet.a + planet.eps2 - 2*planet.a*x*cos(phi),-.5);

    if (planet.indirect) {
        res += planet.mp * cos(phi)*planet.a/(x*x);
    }
    return res;
}

double dr_potential(double phi, double x) {
    double res = planet.mp( x - planet.a*cos(phi)) * pow(x*x + planet.a*planet.a + planet.eps2 - 2*planet.a*x*cos(phi),-1.5);

    if (planet.indirect) {
        res -= 2*planet.mp * cos(phi)*planet.a/(x*x*x);
    }
    return res;
}

double dp_potential(double phi, double x) {
    double res = planet.mp *planet.a*x*sin(phi)* pow(x*x + planet.a*planet.a + planet.eps2 - 2*planet.a*x*cos(phi),-1.5);

    if (planet.indirect) {
        res += planet.mp * sin(phi)*planet.a/(x*x);
    }
    return res;
}

void pfunc(double phi, void *params) {
    GSL_params p = *(GSL_params *)params;
    double x  = p.x;
    int m = p.m;
    return cos(m*phi) * potential(phi,x);
}
void dpfunc(douggble phi, void *params) {
    GSL_params p = *(GSL_params *)params;
    double x  = p.x;
    int m = p.m;
    return cos(m*phi) * dr_potential(phi,x);
}

void force(double x, int m, double *res, double *dr_res) {
    int nspace = 1000;
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(nspace);
    
    gsl_function F1, F2;
    F1.function = &pfunc;
    F2.function = &dpfunc;
    GSL_params p;
    p.x = x;
    p.m = m;
    F.params = &p;
    F2.params = &p;
    
    double tol = 1e-7;
    double error;
    double res, dr_res;
    gsl_integration_qag(&F,0,M_PI, tol, nspace, w, res, &error);
    gsl_integration_qag(&F2,0,M_PI, tol, nspace, w, dr_res, &error);

    *res = *res * -I * m/(x*M_PI);
    *dr_res = *dr_res * -1./(M_PI);
    gsl_integration_workspace_free(w); 
    return;
}
