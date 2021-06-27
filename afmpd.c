// Magnetic Phase Diagram - AFM/PM Transition

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

int itv = 128; // interval

void EigenCal(double w[N], v[LDVL*N], double kx, double ky, double U, double target_n, double n1_up, double n2_up) { // Eigenvalue, eigenvector calculator
	/*
	   a[0]=U<n1_up>	a[2]=tr
	   a[1]=tr*			a[3]=U<n2_up>

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
	for(i=0; i<LDVR*N; i++) v[i] = vr[i];

	free(work);
}

double DoubleMuCal(double target_n, double U) { // Double cell mu calculator
	int i, j, x, y;
	double n, n1_up, n2_up, mu, kx, ky, sum;

	double w[N], v[LDVR*N];
	
	n = 1;
	n1_up = n2_up = target_n/4;
	mu = U;

	while(target_n < n) {
		sum = 0;

		for(x=0; x<itv; x++) {
			kx = -pi + (2*pi*x/(double)itv);
			for(y=0; y<itv; y++) {
				 ky = -pi + (2*pi*y/(double)itv);

				 EigenCal(w, v, kx, ky, U, target_n, n1_up, n2_up);

				 for(i=0; i<N; i++) {
					 if(w[i] < mu) {
						 for(j=i*N; j<LDVR+i*N; j++) sum += pow(v[j], 2);
					 }
				 }
			}
		}
		n = (double)sum/(itv*itv);
		mu -= 0.01;
	}

	return mu;
}

double DoubleMCal(double *n, double U, double mu) { // Double cell m calculator
	int i, x, y, up_cnt, down_cnt;
	double m, m1, m2, kx, ky;

	m1 = 0.1;
	m2 = 0.1;

	for(i=0; i<128; i++) {
		up_cnt = 0;
		down_cnt = 0;

		for(x=0; x<itv; x++) {
			kx = -pi + (2*pi*x/(double)itv);
			for(y=0; y<itv; y++) {
				ky = -pi + (2*pi*y/(double)itv);

				//if(DUE(kx, ky, U, *n, m1, m2) < mu) up_cnt++;
				//if(DDE(kx, ky, U, *n, m1, m2) < mu) down_cnt++;
			}
		}
		m = ((double)up_cnt/(itv*itv) - (double)down_cnt/(itv*itv))/2;
		m2 = m - m1;
	}

	return m;
}

int main() {
	FILE *fp;

	double target_n, U, mu, m, time;

	fp = fopen("data/afmpd.txt", "w");
	fprintf(fp, "target_n\tt/U\n");

	//printf("target_n\tt/U\telapsed time(s)\n");
	
	for(target_n=1.0; target_n>0.1; target_n-=0.1) {
		time = clock();

		for(U=5; U<10; U+=1) {
			mu = DoubleMuCal(target_n, U);
			//m = DoubleMCal(target_n, U, mu);

			//if(fabs(m) > 1e-1) break;
		}
		//printf("%.3f\t%f\t%.3f\n", target_n, 1/U, (clock()-time)*0.000001);
		//fprintf(fp, "%f\t%f\n", target_n, 1/U);
	}

	fclose(fp);

	return 0;
}
