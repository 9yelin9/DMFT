// dgeev test

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include </usr/include/triqs/arrays/blas_lapack/f77/lapack.h>

#define pi 3.141592653689793

// dgeev_ parameters
#define N 2 // The order of the matrix
#define LDA N // The leading dimensionof the array A
#define LDVL N // The leading dimension of the array VL(left eigenvectors)
#define LDVR N // The leading dimension of the array VR(right eigenvectors)

#define TEST_PE(kx, ky, U, target_n, n1_up, n2_up) (0.5*(U*(n1_up+n2_up) + sqrt(pow(U, 2)*pow(n1_up+n2_up, 2)-4*(pow(U, 2)*n1_up*n2_up-pow(2*(cos(kx)+cos(ky)), 2)))))
#define TEST_ME(kx, ky, U, target_n, n1_up, n2_up) (0.5*(U*(n1_up+n2_up) - sqrt(pow(U, 2)*pow(n1_up+n2_up, 2)-4*(pow(U, 2)*n1_up*n2_up-pow(2*(cos(kx)+cos(ky)), 2)))))

int itv = 128; // interval

void EigenCal(double w[N], double kx, double ky, double U, double target_n, double n1_up, double n2_up) { // Eigenvalue, eigenvector calculator
	/*
	   a[0] = U<n1_up>	a[1] = tr
	   a[2] = tr*		a[3] = U<n2_up>

	   eigenvalue e에 대하여
	   (U<n1_up>-e)(U<n2_up>-e) - (t^2*|r|^2) = 0

	   따라서 r는 절대값을 사용해도 됨
	*/
	int n = N, lda = LDA, ldvl = LDVL, ldvr = LDVR;
	int info; // = 0 : successful exit
	int lwork; // The dimension of the array WORK
	double wkopt; // WORK optimal
	double *work;

	double a[LDA*N];
	double wr[N], wi[N]; // eigenvalues(WR : real parts, WI : imaginary parts)
	double vl[LDVL*N], vr[LDVR*N]; // eigenvectors(VL : left eigenvectors, VR : right eigenvectors)
	
	a[0] = U*n1_up;
	a[1] = 2.0*(cos(kx)+cos(ky));
	a[2] = 2.0*(cos(kx)+cos(ky));
	a[3] = U*n2_up;

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
	for(i=0; i<N; i++) w[i] = wr[i];
	
	free(work);
}

int main() {
	FILE *fp;

	int x, y;
	double target_n, n1_up, n2_up, U, kx, ky;

	double w[N];

	fp = fopen("data/dgeev_test.txt", "w");
    fprintf(fp, "wr[0]\twr[1]\n");
    printf("wr[0]\twr[1]\n");

	target_n = 1;
	n1_up = n2_up = target_n/4;
	U = 1;

	for(x=0; x<itv; x++) {
		kx = -pi + (2*pi*x/(double)itv);
		for(y=0; y<itv; y++) {
			ky = -pi + (2*pi*y/(double)itv);

			EigenCal(w, kx, ky, U, target_n, n1_up, n2_up);
			printf("%f\t%f\n", w[0], w[1]);
			fprintf(fp, "%f\t%f\n", w[0], w[1]);
			
			//printf("%f\t%f\n\n", TEST_PE(kx, ky, U, target_n, n1_up, n2_up), TEST_ME(kx, ky, U, target_n, n1_up, n2_up));
			
			//printf("error : %f\n", wr[0]-TEST_PE(kx, ky, U, target_n, n1_up, n2_up));
		}
	}

	fclose(fp);

	return 0;
}
