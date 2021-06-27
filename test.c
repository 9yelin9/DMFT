// dgeev test

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include </usr/include/triqs/arrays/blas_lapack/f77/lapack.h>

// dgeev_ parameters
#define N 2 // The order of the matrix
#define LDA N // The leading dimensionof the array A
#define LDVL N // The leading dimension of the array VL(left eigenvectors)
#define LDVR N // The leading dimension of the array VR(right eigenvectors)

int main() {
	int n = N, lda = LDA, ldvl = LDVL, ldvr = LDVR;
	int info; // = 0 : successful exit
	int lwork; // The dimension of the array WORK
	double wkopt; // WORK optimal
	double *work;

	double a[LDA*N] = {
		3, 1,
		2, 4
	};
	double wr[N], wi[N]; // eigenvalues(WR : real parts, WI : imaginary parts)
	double vl[LDVL*N], vr[LDVR*N]; // eigenvectors(VL : left eigenvectors, VR : right eigenvectors)
	
	// Query and allocate the optimal workspace
	lwork = -1;
	dgeev_("V", "V", &n, a, &lda, wr, wi, vl, &ldvl, vr, &ldvr, &wkopt, &lwork, &info);
	lwork = (int)wkopt;
	work = (double*)malloc(lwork*sizeof(double));

	// Solve eigenproblem
	dgeev_("V", "V", &n, a, &lda, wr, wi, vl, &ldvl, vr, &ldvr, work, &lwork, &info);

	if(info > 0) {
		printf("dgeev_ FAIL\n");
		exit(1);
	}

	int i;
	printf("%f\t%f\n", wr[0], wr[1]);
	for(i=0; i<LDVR*N; i++) printf("%f\t", vr[i]);
	printf("\n");
	
	free(work);
}
