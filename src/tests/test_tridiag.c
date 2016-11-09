#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void matmat(double  *A, double *B, double *C,
					double alpha, double beta, int nA)
{
/* Performs \alpha * A.B + \beta * C and stores the output in C.
	A,B, and C are all matrices.
	This is essenitally a wrapper for the ZGEMM BLAS routine
*/
	int i,j;
	char TRANSA = 't';
	char TRANSB = 't';
	int m = nA;
	int n = nA;
	int k = nA;
	int LDA = nA;
	int LDB = nA;
	int LDC = nA;

    double *work = (double *)malloc(sizeof(double)*nA*nA);
	for(i=0;i<nA;i++) {
		for(j=0;j<nA;j++) work[i+nA*j] = C[j + nA*i];
	}

	for(i=0;i<nA;i++) {
		for(j=0;j<nA;j++)	C[i + nA*j] = work[i+nA*j];
	}

	dgemm_(&TRANSA, &TRANSB, &m,&n,&k,&alpha,A,&LDA,B,&LDB,&beta,C,&LDC);


	for(i=0;i<nA;i++) {
		for(j=0;j<nA;j++)  work[i+nA*j] = C[j + nA*i];
	}

	for(i=0;i<nA;i++) {
		for(j=0;j<nA;j++)	C[i + nA*j] = work[i+nA*j];
	}
    free(work);
	return;

}

void matvec(double  *A, double*B, double *C,
					double alpha, double beta, int nB)
{
/* Performs \alpha * A.B + \beta * C and stores the output in C.
	A is a matrix, B and C are vectors.
	This is essenitally a wrapper for the ZGEMV BLAS routine
*/

	char TRANS = 't';
	int m = nB;
	int n = nB;
	int LDA = nB;
	int INCX = 1;
	int INCY = 1;



	dgemv_(&TRANS, &m,&n,&alpha,A,&LDA,B,&INCX,&beta,C,&INCY);

	return;

}

void matmat3(double *A, double *B, double *C, double alpha, double beta) {

    C[0] = beta*C[0] + alpha*(A[0]*B[0] + A[1]*B[3] + A[2]*B[6]);
    C[1] = beta*C[1] + alpha*(A[0]*B[1] + A[1]*B[4] + A[2]*B[7]);
    C[2] = beta*C[2] + alpha*(A[0]*B[2] + A[1]*B[5] + A[2]*B[8]);

    C[3] = beta*C[3] + alpha*(A[3]*B[0] + A[4]*B[3] + A[5]*B[6]);
    C[4] = beta*C[4] + alpha*(A[3]*B[1] + A[4]*B[4] + A[5]*B[7]);
    C[5] = beta*C[5] + alpha*(A[3]*B[2] + A[4]*B[5] + A[5]*B[8]);

    C[6] = beta*C[6] + alpha*(A[6]*B[0] + A[7]*B[3] + A[8]*B[6]);
    C[7] = beta*C[7] + alpha*(A[6]*B[1] + A[7]*B[4] + A[8]*B[7]);
    C[8] = beta*C[8] + alpha*(A[6]*B[2] + A[7]*B[5] + A[8]*B[8]);

    return;
}

void matvec3(double *A, double *B, double *C, double alpha, double beta) {

    C[0] = beta*C[0] + alpha*(A[0]*B[0] + A[1]*B[1] + A[2]*B[2]);
    C[1] = beta*C[1] + alpha*(A[3]*B[0] + A[4]*B[1] + A[5]*B[2]);
    C[2] = beta*C[2] + alpha*(A[6]*B[0] + A[7]*B[1] + A[8]*B[2]);
    
    return;
}


void solve(double *A, double *B,int nA, int nRHS) {

    int N = nA;
    int NRHS = nRHS;
    int LDA = nA;
    int *IPIV = (int *)malloc(sizeof(int)*nA);
    int LDB = nA;
    int INFO;

    double *AT = (double *)malloc(sizeof(double)*nA*nA);
    double *BT = (double *)malloc(sizeof(double)*nA*nRHS);
    int i,j;
    for(i=0;i<nA;i++) {
        for(j=0;j<nA;j++) {
            AT[i + nA*j] = A[j + nA*i];
        }
    }
    if (nRHS > 1) {
        for(i=0;i<nA;i++) {
            for(j=0;j<nA;j++) {
                BT[i + nA*j] = B[j + nA*i];
            }
        }
    }
    else {
        for(i=0;i<nA;i++) BT[i] = B[i];
    }

    dgesv_(&nA, &NRHS, AT, &LDA, IPIV, BT, &LDB, &INFO);

    if (nRHS > 1) {
        for(i=0;i<nA;i++) {
            for(j=0;j<nA;j++) {
                B[i + nA*j] = BT[j + nA*i];
            }
        }
    }
    else {
        for(i=0;i<nA;i++) B[i] = BT[i];
    }
    free(AT);
    free(BT);
    return;
}

void thomas_alg(double *a, double *b, double *c, double *d, int n) {
    /* 
     * Solve the tridiagonal system a_i*y_{i-1} + b_i*y_i + c_i*y_{i+1} = d_i
     * using Thomas's algorithm
     * Solution is stored in d array
     */
    int i; 

    c[0] /= b[0];
    d[0] /= b[0];

    for(i=1;i<n-1;i++) {
        
        b[i] -= a[i-1] * c[i-1];

        d[i] -= a[i-1]*d[i-1];


        c[i] = c[i] / b[i];
        d[i] = d[i] / b[i];
        

    }
    i = n-1;
    b[i] -= a[i-1] * c[i-1];
    d[i] -= a[i-1]*d[i-1];
    d[i] = d[i] / b[i];


    for(i=n-2;i>=0;i--) {
        d[i] -= c[i]*d[i+1];
    }

    return;
}
void thomas_alg_block(double *a, double *b, double *c, double *d, int n, int m) {
    /* 
     * Solve the block tridiagonal system a_i*y_{i-1} + b_i*y_i + c_i*y_{i+1} = d_i
     * using Thomas's algorithm
     * Solution is stored in d array
     */
    int i; 
    int size = m*m;

    solve(&b[0],&c[0],m,m);
    solve(&b[0],&d[0],m,1);

    for(i=1;i<n-1;i++) {
        
        matmat3(&a[(i-1)*size],&c[(i-1)*size],&b[i*size],-1.0,1.0);
        matvec3(&a[(i-1)*size],&d[(i-1)*m],&d[i*m],-1.0,1.0);

        solve(&b[i*size],&c[i*size],m,m);
        solve(&b[i*size],&d[i*m],m,1);
        

    }
    i = n-1;
        
    matmat3(&a[(i-1)*size],&c[(i-1)*size],&b[i*size],-1.0,1.0);
    matvec3(&a[(i-1)*size],&d[(i-1)*m],&d[i*m],-1.0,1.0);
    solve(&b[i*size],&d[i*m],m,1);


    for(i=n-2;i>=0;i--) {
        matvec3(&c[i*size],&d[(i+1)*m],&d[i*m],-1,1);
    }

    return;
}

void test_1(int n) {

    FILE *f = fopen("input.dat","r");

    double *md = (double *)malloc(sizeof(double)*n);
    double *ud = (double *)malloc(sizeof(double)*(n-1));
    double *ld = (double *)malloc(sizeof(double)*(n-1));
    double *fd = (double *)malloc(sizeof(double)*n);
    double *ans = (double *)malloc(sizeof(double)*n);

    fread(md, sizeof(double),n,f);
    fread(ud, sizeof(double),(n-1),f);
    fread(ld, sizeof(double),(n-1),f);
    fread(fd, sizeof(double),n,f);
    fread(ans, sizeof(double),n,f);


    thomas_alg(ld,md,ud,fd,n);


    int i;
    double err = 0;
    for(i=0;i<n;i++) err += (fd[i]-ans[i])*(fd[i]-ans[i]);

    err = pow(err/n,.5);

    printf("Error is %.5e\n",err);



    free(md);free(ud);free(ld);free(fd);free(ans);

    return;

}

void test_2(void) {
    int m = 3;
    int size = 9;
    double *A = (double *)malloc(sizeof(double)*size);
    double *B = (double *)malloc(sizeof(double)*size);
    double *C = (double *)malloc(sizeof(double)*size);
    double *ansm = (double *)malloc(sizeof(double)*size);
    double *anssm = (double *)malloc(sizeof(double)*size);

    double *D = (double *)malloc(sizeof(double)*m);
    double *E = (double *)malloc(sizeof(double)*m);
    double *ansv = (double *)malloc(sizeof(double)*m);
    double *anss = (double *)malloc(sizeof(double)*m);

    double coeffs[2];

    FILE *f = fopen("input.dat","r");
    fread(A,sizeof(double),size,f);
    fread(B,sizeof(double),size,f);
    fread(C,sizeof(double),size,f);
    fread(D,sizeof(double),m,f);
    fread(E,sizeof(double),m,f);
    fread(&coeffs[0],sizeof(double),2,f);
    fread(ansm,sizeof(double),size,f);
    fread(ansv,sizeof(double),m,f);
    fread(anss,sizeof(double),m,f);
    fread(anssm,sizeof(double),size,f);


    int i;


    matmat(A,B,C,coeffs[0],coeffs[1],m);
    matvec(A,D,E,coeffs[0],coeffs[1],m);
    //matmat3(A,B,C,coeffs[0],coeffs[1]);
    //matvec3(A,D,E,coeffs[0],coeffs[1]);
    solve(A,D,m,1);
    solve(A,B,m,m);

    double err_m = 0;
    double err_v = 0;
    double err_s = 0;
    double err_sm = 0;
    for(i=0;i<size;i++) {
        err_m += (C[i] - ansm[i])*(C[i]-ansm[i]);
    }
    for(i=0;i<m;i++) {
        err_v += (E[i] - ansv[i])*(E[i]-ansv[i]);
    }
    for(i=0;i<m;i++) {
        //printf("%lg\t%lg\n",D[i],anss[i]);
        err_s += (D[i] - anss[i])*(D[i]-anss[i]);
    }
    for(i=0;i<size;i++) {
        err_sm += (B[i] - anssm[i])*(B[i]-anssm[i]);
    }
    int j;
    for(i=0;i<m;i++) {
        for(j=0;j<m;j++) {
            printf("%lg\t%lg\n",B[i+m*j],anssm[j+i*m]);
        }
    }
    err_m = pow(err_m/size,.5);
    err_v = pow(err_v/m,.5);
    err_s = pow(err_s/m,.5);
    err_sm = pow(err_sm/m,.5);

    printf("Matmat Error: %.5e\n", err_m);
    printf("Matvec Error: %.5e\n", err_v);
    printf("Solve Error: %.5e\n", err_s);
    printf("3 Solve Error: %.5e\n", err_sm);

    free(A); free(B); free(C); free(D); free(E); free(ansm); free(ansv);
    free(anss);free(anssm);
    return;
}

void test_3(int n, int m) {

    int size = m*m;
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



    int i;
    for(i=0;i<n*m;i++) printf("%lg\n",fd[i]);

    thomas_alg_block(ld,md,ud,fd,n,m);
    for(i=0;i<n*m;i++) printf("%lg\t%lg\n",fd[i],ans[i]);


    double err = 0;
    for(i=0;i<n*m;i++) err += (fd[i]-ans[i])*(fd[i]-ans[i]);

    err = pow(err/(n*m),.5);

    printf("Error is %.5e\n",err);



    free(md);free(ud);free(ld);free(fd);free(ans);
    return;
}

int main(int argc, char *argv[]) {
    int i;    
    if (argc == 1) {
        test_1(20);
    }
    else {
        switch (atoi(argv[1])) {
            case 1:
                test_1(20);
                break;
            case 2:
                test_2();
                break;
            case 3:
                test_3(20,3);
                break;
        }

    }


    return 1;
}
