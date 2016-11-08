
void normalize_evectors(double complex *evecs);

void reigenvalues(double complex *A, double complex *Q, double complex *evals, double complex *evecs, int nA, int generalized)
{
/* Computes the eigenvalues and right eigenvectors of the complex matrix A which is NxN.
	A . evecs = evals Q . evecs
	This is essentially a wrapper for the ZGEEV LAPACK routine.

	INPUTS:
		The matrices M, Q, and evecs and the vector evals which are all overwritten

	OUTPUTS:
		The eigenvalues are stored in the evals array.
		The eigenvectors are stored in the ROWS of the evecs matrix
*/
	int i,j;
	char JOBVL = 'N';
	char JOBVR = 'V';
	int INFO;
	int LDA = nA;
	int LDB = nA;
	int LDVL = nA;
	int LDVR = nA;
	int LWORK = 2*nA;

	double *RWORK = (double *)malloc(sizeof(double)*8*nA);
	double complex *CWORK = (double complex *)malloc(sizeof(double complex)*2*nA);

	double complex *tA = (double complex *)malloc(sizeof(double complex)*nA*nA);
	double complex *tQ = (double complex *)malloc(sizeof(double complex)*nA*nA);

	double complex *evals_alpha = (double complex *)malloc(sizeof(double complex)*nA);
	double complex *evals_beta = (double complex *)malloc(sizeof(double complex)*nA);

	for(i=0;i<nA;i++) {
		for(j=0;j<nA;j++) {
	  		tA[i+nA*j] = A[j+nA*i];
	  		tQ[i+nA*j] = Q[j+nA*i];
		}
	}


if (generalized) {
	zggev_( &JOBVL, &JOBVR, &nA, tA, &LDA, tQ, &LDB, evals_alpha,evals_beta, NULL, &LDVL, evecs, &LDVR, CWORK, &LWORK, RWORK, &INFO );
	for(i=0;i<nA;i++) {
		if (cabs(evals_beta[i]) != 0) {
			evals[i] = evals_alpha[i]/evals_beta[i];
		}
	}
}
else {
	zgeev_( &JOBVL, &JOBVR, &nA, tA, &LDA, evals, tQ, &LDVL, evecs, &LDVR, CWORK, &LWORK, RWORK, &INFO );
}


//	normalize_evectors(evecs);


	free(tA); free(tQ);
	free(RWORK); free(CWORK);
	free(evals_alpha);
	free(evals_beta);
 	return;
}

void solve(double *A, double *B,int nA) {

    int N = nA;
    int NRHS = 1;
    int LDA = nA;
    int *IPIV = (int *)malloc(sizeof(int)*nA);
    int LDB = NA;
    int INFO;

    dgesv_(&nA, &NRHS, A, &LDA, IPIV, B, &LDB, &INFO);

    return;
}

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


	for(i=0;i<nA;i++) {
		for(j=0;j<nA;j++) work[i+N*j] = C[j + nA*i];
	}

	for(i=0;i<nA;i++) {
		for(j=0;j<nA;j++)	C[i + nA*j] = work[i+N*j];
	}

	dgemm_(&TRANSA, &TRANSB, &m,&n,&k,&alpha,A,&LDA,B,&LDB,&beta,C,&LDC);


	for(i=0;i<nA;i++) {
		for(j=0;j<nA;j++)  work[i+N*j] = C[j + nA*i];
	}

	for(i=0;i<nA;i++) {
		for(j=0;j<nA;j++)	C[i + nA*j] = work[i+N*j];
	}
	return;

}

void matmat_simple(double *A, double *B, double *C, double alpha, double beta,int n) {
 int i,j;

 for(i=0;i<n;i++) {
     for(j=0;j<n;j++) {
         res = 0;

         for(k=0;k<n;k++) {
            res += a[k+ i*n]*b[j + k*n];
         }
         c[j + i*n] = alpha * res + beta * c[j+i*n];
     }
 }

    return;
}

void cmatmat(double complex *A, double complex *B, double complex *C,
					double complex alpha, double complex beta, int nA)
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


	for(i=0;i<nA;i++) {
		for(j=0;j<nA;j++) cwork[i+N*j] = C[j + nA*i];
	}

	for(i=0;i<nA;i++) {
		for(j=0;j<nA;j++)	C[i + nA*j] = cwork[i+N*j];
	}

	zgemm_(&TRANSA, &TRANSB, &m,&n,&k,&alpha,A,&LDA,B,&LDB,&beta,C,&LDC);


	for(i=0;i<nA;i++) {
		for(j=0;j<nA;j++)  cwork[i+N*j] = C[j + nA*i];
	}

	for(i=0;i<nA;i++) {
		for(j=0;j<nA;j++)	C[i + nA*j] = cwork[i+N*j];
	}
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
void solve(double *A, double *B,int nA) {

    int N = nA;
    int NRHS = 1;
    int LDA = nA;
    int *IPIV = (int *)malloc(sizeof(int)*nA);
    int LDB = nA;
    int INFO;

    double *AT = (double *)malloc(sizeof(double)*nA*nA);
    int i,j;
    for(i=0;i<nA;i++) {
        for(j=0;j<nA;j++) {
            AT[i + nA*j] = A[j + nA*i];
        }
    }

    dgesv_(&nA, &NRHS, AT, &LDA, IPIV, B, &LDB, &INFO);
    free(AT);
    return;
}
void csolve(double complex *A, double complex *B,int nA) {

    int N = nA;
    int NRHS = 1;
    int LDA = nA;
    int *IPIV = (int *)malloc(sizeof(int)*nA);
    int LDB = nA;
    int INFO;

    double complex *AT = (double complex *)malloc(sizeof(double complex)*nA*nA);
    int i,j;
    for(i=0;i<nA;i++) {
        for(j=0;j<nA;j++) {
            AT[i + nA*j] = A[j + nA*i];
        }
    }

    zgesv_(&nA, &NRHS, AT, &LDA, IPIV, B, &LDB, &INFO);
    free(AT);
    return;
}

void cmatvec(double  complex *A, double complex *B, double complex *C,
					double complex alpha, double complex beta, int nB)
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



	zgemv_(&TRANS, &m,&n,&alpha,A,&LDA,B,&INCX,&beta,C,&INCY);

	return;

}


//
// void normalize_evectors(double complex *evecs) {
// /* Normalize the eigenvectors */
// /* Calculate the factor to normalize the disk eccentricity.
// 	 Each planet eccentricity will then be normalized by the same factor.
// */
// 	int i,j,indx;
// 	double norm;
//
//
// #ifdef OPENMP
// #pragma omp parallel private(i,j,norm,indx) shared(evecs,nrows,ncols)
// #pragma omp for schedule(static)
// #endif
// 	for(i=0;i<nrows;i++) {
//
//
// 			norm = 0;
// #ifdef NORMALIZE_INT
// 			for(j=0;j<N;j++) {
// 				indx = j + ncols*i;
//
// 				norm += weights[j]*conj(evecs[indx])*evecs[indx];
//
// 			}
// 			norm = sqrt(norm);
// #else
// #ifdef NORMALIZE_MAX
// 			for(j=0;j<N;j++) {
// 				indx = j+ncols*i;
// //				printf("%lg\t%lg\t%lg",norm,abs(evecs[indx]),fmax(norm,abs(evecs[indx])));
// 				norm = fmax(norm,abs(evecs[indx]));
// 			}
// 			if (norm == 0) norm = 1;
//
//
// #endif
// #endif
// 			for(j=0;j<ncols;j++) {
// 					indx = j + ncols*i;
// 					evecs[indx] /= norm;
// 			}
// 	}
//
// 	return;
//
//
// }

void normalize_evectors(double complex *evecs) {
/* Normalize the eigenvectors */
/* Calculate the factor to normalize the disk eccentricity.
	 Each planet eccentricity will then be normalized by the same factor.
*/
	int i,j,indx;
	double norm;


#ifdef OPENMP
#pragma omp parallel private(i,j,norm,indx) shared(evecs,nrows,ncols)
#pragma omp for schedule(static)
#endif
	for(i=0;i<nrows;i++) {


			norm = 0;
#ifdef NORMALIZE_NORM
		for(j=0;j<N;j++) {
			indx = j +ncols*i;
			norm += conj(evecs[indx])*evecs[indx];
		}
		norm = sqrt(norm);


#else
#ifdef NORMALIZE_INT
			for(j=0;j<N;j++) {
				indx = j + ncols*i;

				norm += weights[j]*conj(evecs[indx])*evecs[indx];

			}
			norm = sqrt(norm);
#else
#ifdef NORMALIZE_MAX
			for(j=0;j<N;j++) {
				indx = j+ncols*i;
//				printf("%lg\t%lg\t%lg",norm,abs(evecs[indx]),fmax(norm,abs(evecs[indx])));
				norm = fmax(norm,abs(evecs[indx]));
			}


#endif
#endif
#endif
			if (norm == 0) norm = 1;
			for(j=0;j<ncols;j++) {
					indx = j + ncols*i;
					evecs[indx] /= norm;
			}
	}

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
void matvec4(double *A, double *B, double *C, double alpha, double beta) {

    C[0] = beta*C[0] + alpha*(A[0]*B[0]  + A[1]*B[1]  + A[2]*B[2]  + A[3]*B[3]);
    C[1] = beta*C[1] + alpha*(A[4]*B[0]  + A[5]*B[1]  + A[6]*B[2]  + A[7]*B[3]);
    C[2] = beta*C[2] + alpha*(A[8]*B[0]  + A[9]*B[1]  + A[10]*B[2] + A[11]*B[3]);
    C[3] = beta*C[3] + alpha*(A[12]*B[0] + A[13]*B[1] + A[14]*B[2] + A[15]*B[3]);
    
    return;
}


void matmat4(double *A, double *B, double *C, double alpha, double beta) {

    C[0] = beta*C[0] + alpha*(A[0]*B[0] + A[1]*B[4] + A[2]*B[8] +  A[3]*B[12]);
    C[1] = beta*C[1] + alpha*(A[0]*B[1] + A[1]*B[5] + A[2]*B[9] +  A[3]*B[13]);
    C[2] = beta*C[2] + alpha*(A[0]*B[2] + A[1]*B[6] + A[2]*B[10] + A[3]*B[14]);
    C[3] = beta*C[3] + alpha*(A[0]*B[3] + A[1]*B[7] + A[2]*B[11] + A[3]*B[15]);

    C[4] = beta*C[4] + alpha*(A[4]*B[0] + A[5]*B[4] + A[6]*B[8] +  A[7]*B[12]);
    C[5] = beta*C[5] + alpha*(A[4]*B[1] + A[5]*B[5] + A[6]*B[9] +  A[7]*B[13]);
    C[6] = beta*C[6] + alpha*(A[4]*B[2] + A[5]*B[6] + A[6]*B[10] + A[7]*B[14]);
    C[7] = beta*C[7] + alpha*(A[4]*B[3] + A[5]*B[7] + A[6]*B[11] + A[7]*B[15]);

    C[8] = beta*C[8]   + alpha*(A[8]*B[0] + A[9]*B[4] + A[10]*B[8]  + A[11]*B[12]);
    C[9] = beta*C[9]   + alpha*(A[8]*B[1] + A[9]*B[5] + A[10]*B[9]  + A[11]*B[13]);
    C[10] = beta*C[10] + alpha*(A[8]*B[2] + A[9]*B[6] + A[10]*B[10] + A[11]*B[14]);
    C[11] = beta*C[11] + alpha*(A[8]*B[3] + A[9]*B[7] + A[10]*B[11] + A[11]*B[15]);

    C[12] = beta*C[12] + alpha*(A[12]*B[0] + A[13]*B[4] + A[14]*B[8]  + A[15]*B[12]);
    C[13] = beta*C[13] + alpha*(A[12]*B[1] + A[13]*B[5] + A[14]*B[9]  + A[15]*B[13]);
    C[14] = beta*C[14] + alpha*(A[12]*B[2] + A[13]*B[6] + A[14]*B[10] + A[15]*B[14]);
    C[15] = beta*C[15] + alpha*(A[12]*B[3] + A[13]*B[7] + A[14]*B[11] + A[15]*B[15]);
    return;
}

void cmatmat3(double complex *A, double complex *B, double complex *C, double complex alpha, double complex beta) {

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

void cmatvec3(double complex *A, double complex *B, double complex *C, double complex alpha, double complex beta) {

    C[0] = beta*C[0] + alpha*(A[0]*B[0] + A[1]*B[1] + A[2]*B[2]);
    C[1] = beta*C[1] + alpha*(A[3]*B[0] + A[4]*B[1] + A[5]*B[2]);
    C[2] = beta*C[2] + alpha*(A[6]*B[0] + A[7]*B[1] + A[8]*B[2]);
    
    return;
}
void cmatvec4(double complex *A, double complex *B, double complex *C, double complex alpha, double complex beta) {

    C[0] = beta*C[0] + alpha*(A[0]*B[0]  + A[1]*B[1]  + A[2]*B[2]  + A[3]*B[3]);
    C[1] = beta*C[1] + alpha*(A[4]*B[0]  + A[5]*B[1]  + A[6]*B[2]  + A[7]*B[3]);
    C[2] = beta*C[2] + alpha*(A[8]*B[0]  + A[9]*B[1]  + A[10]*B[2] + A[11]*B[3]);
    C[3] = beta*C[3] + alpha*(A[12]*B[0] + A[13]*B[1] + A[14]*B[2] + A[15]*B[3]);
    
    return;
}


void cmatmat4(double complex *A, double complex *B, double complex *C, double complex alpha, double complex beta) {

    C[0] = beta*C[0] + alpha*(A[0]*B[0] + A[1]*B[4] + A[2]*B[8] +  A[3]*B[12]);
    C[1] = beta*C[1] + alpha*(A[0]*B[1] + A[1]*B[5] + A[2]*B[9] +  A[3]*B[13]);
    C[2] = beta*C[2] + alpha*(A[0]*B[2] + A[1]*B[6] + A[2]*B[10] + A[3]*B[14]);
    C[3] = beta*C[3] + alpha*(A[0]*B[3] + A[1]*B[7] + A[2]*B[11] + A[3]*B[15]);

    C[4] = beta*C[4] + alpha*(A[4]*B[0] + A[5]*B[4] + A[6]*B[8] +  A[7]*B[12]);
    C[5] = beta*C[5] + alpha*(A[4]*B[1] + A[5]*B[5] + A[6]*B[9] +  A[7]*B[13]);
    C[6] = beta*C[6] + alpha*(A[4]*B[2] + A[5]*B[6] + A[6]*B[10] + A[7]*B[14]);
    C[7] = beta*C[7] + alpha*(A[4]*B[3] + A[5]*B[7] + A[6]*B[11] + A[7]*B[15]);

    C[8] = beta*C[8]   + alpha*(A[8]*B[0] + A[9]*B[4] + A[10]*B[8]  + A[11]*B[12]);
    C[9] = beta*C[9]   + alpha*(A[8]*B[1] + A[9]*B[5] + A[10]*B[9]  + A[11]*B[13]);
    C[10] = beta*C[10] + alpha*(A[8]*B[2] + A[9]*B[6] + A[10]*B[10] + A[11]*B[14]);
    C[11] = beta*C[11] + alpha*(A[8]*B[3] + A[9]*B[7] + A[10]*B[11] + A[11]*B[15]);

    C[12] = beta*C[12] + alpha*(A[12]*B[0] + A[13]*B[4] + A[14]*B[8]  + A[15]*B[12]);
    C[13] = beta*C[13] + alpha*(A[12]*B[1] + A[13]*B[5] + A[14]*B[9]  + A[15]*B[13]);
    C[14] = beta*C[14] + alpha*(A[12]*B[2] + A[13]*B[6] + A[14]*B[10] + A[15]*B[14]);
    C[15] = beta*C[15] + alpha*(A[12]*B[3] + A[13]*B[7] + A[14]*B[11] + A[15]*B[15]);
    return;
}
