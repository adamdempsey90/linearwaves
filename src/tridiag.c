#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void thomas_alg(double *a, double *b, double *c, double *d, int n) {
    /* 
     * Solve the tridiagonal system a_i*y_{i-1} + b_i*y_i + c_i*y_{i+1} = d_i
     * using Thomas's algorithm
     * Solution is stored in d array
     */
    int i; 
    double fac, fac2;
    c[0] /= b[0];
    d[0] /= b[0];

    for(i=1;i<n-1;i++) {
        
        fac = b[i] - a[i-1]*c[i-1];
        fac2 = d[i] - a[i-1]*d[i-1];
        c[i] = c[i] / fac;
        d[i] = fac2 / fac;

    }
    i = n-1;
    fac = b[i] - a[i-1]*c[i-1];
    fac2 = d[i] - a[i-1]*d[i-1];
    d[i] = fac2 / fac;

    d[n-1] = d[n-1];

    for(i=n-2;i>=0;i--) {
        d[i] = d[i] - c[i]*d[i+1];
    }

    return;
}


int main(void) {

    int n = 20;
    int m = 3;
    int size = m*m*;
    FILE *f = fopen("input.dat","r");

    double *md = (double *)malloc(sizeof(double)*n*size);
    double *ud = (double *)malloc(sizeof(double)*(n-1)*size);
    double *ld = (double *)malloc(sizeof(double)*(n-1)*size);
    double *fd = (double *)malloc(sizeof(double)*n*m);
    double *ans = (double *)malloc(sizeof(double)*n*m);

    fread(md, sizeof(double),n*size,f);
    fread(ud, sizeof(double),(n-1)*size,f);
    fread(ld, sizeof(double),(n-1)*size,f);
    fread(fd, sizeof(double),n*m,f);
    fread(ans, sizeof(double),n*m,f);


    thomas_alg_block(ld,md,ud,fd,n,m);


    int i;
    double err = 0;
    for(i=0;i<n;i++) err += (fd[i]-ans[i])*(fd[i]-ans[i]);

    err = pow(err/n,.5);

    printf("Error is %.5e\n",err);



    free(md);free(ud);free(ld);free(fd);free(ans);
    return 1;
}
