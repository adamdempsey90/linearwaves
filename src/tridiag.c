#include "linearwaves.h"


extern void dgemm_(char *, char *, int *, int *,int *, double *, double *, int *, double *, int *, double *, double *, int *);
extern void dgemv_(char *, int *, int *, double *, double *, int *, double *,int *,double *,double *,int*);
extern void dgesv_(int *, int *, double *, int *, int *, double *, int *, int *);
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

void construct_total_matrix(double *ld, double *md, double *ud, double *mat, int n, int m) {

    int i,j,k;

    int size = m*m;

    i=0;
    for(j=i*m;j<(i+1)*m;j++) {
        for(k=i*m; k < (i+1)*m; k++) {
            mat[k + j*n] = md[i*size + k + j*m];
        }
        for(k=(i+1)*m; k<(i+2)*m;k++) {
            mat[k + j*n] = ud[i*size + k + j*m];
        }
    }


    for(i=1;i<n-1;i++) {

        for(j=i*m;j<(i+1)*m;j++) {
            for(k=i*m; k < (i+1)*m; k++) {
                mat[k + j*n] = md[i*size + k + j*m];
            }
            for(k=(i+1)*m; k<(i+2)*m;k++) {
                mat[k + j*n] = ud[i*size + k + j*m];
            }
            for(k=(i-1)*m;k<i*m;k++) {
                mat[k + j*n] = ld[i*size + k+j*m];
            }
        }
    }
    i = n-1;
    for(j=i*m;j<(i+1)*m;j++) {
        for(k=i*m; k < (i+1)*m; k++) {
            mat[k + j*n] = md[i*size + k + j*m];
        }
        for(k=(i-1)*m;k<i*m;k++) {
            mat[k + j*n] = ld[i*size + k+j*m];
        }
    }


}
