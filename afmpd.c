// Magnetic Phase Diagram - AFM/PM Transition

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <lapack.h>

#define pi 3.141592653689793
#define EU(U, target_n, n1_up, n1_down, n2_up, n2_down) (U*target_n*(n1_up*n1_down + n2_up*n2_down))

// dgeev_ parameters
#define N 2 // The order of the matrix
#define LDA N // The leading dimensionof the array A
#define LDVL N // The leading dimension of the array VL(left eigenvectors)
#define LDVR N // The leading dimension of the array VR(right eigenvectors)

int itv = 128; // interval

// 변수를 구조체로 만들어보기 (종류별로)
void EigenCal(double w[N], double v[LDVL*N], double kx, double ky, double U, double target_n, double n1_down, double n2_down) { // Eigenvalue, eigenvector calculator
	/*
	   a[0]=U<n1_down>	a[2]=tr
	   a[1]=tr*			a[3]=U<n2_down>
	*/
	int n = N, lda = LDA, ldvl = LDVL, ldvr = LDVR;
	int info; // = 0 : successful exit
	int lwork; // The dimension of the array WORK
	double wkopt; // WORK optimal
	double *work;

	double a[LDA*N];
	double wr[N], wi[N]; // eigenvalues(WR : real parts, WI : imaginary parts)
	double vl[LDVL*N], vr[LDVR*N]; // eigenvectors(VL : left eigenvectors, VR : right eigenvectors)
	
	a[0] = U*n1_down;
	a[1] = 2.0*(cos((kx+ky)/2)+cos((ky-kx)/2));
	a[2] = 2.0*(cos((kx+ky)/2)+cos((ky-kx)/2));
	a[3] = U*n2_down;

	// Query and allocate the optimal workspace
	// 독립 함수로 만들기
	lwork = -1;
	dgeev_("V", "V", &n, a, &lda, wr, wi, vl, &ldvl, vr, &ldvr, &wkopt, &lwork, &info);
	lwork = (int)wkopt; // 이걸 얻는 함수
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
	// 반복을 줄이면 좋다

	free(work);
}

double DoubleMuCal(double target_n, double U) { // Double cell mu calculator
	int i, j, x, y;
	double n1_up, n2_up, n1_down, n2_down, mu, kx, ky, sum1, sum2;

	double w[N], v[LDVR*N];
	
	FILE *fp;
	fp = fopen("data/afmpd_energy.txt", "w");
	fprintf(fp, "kx\tky\teigenvalue1\teigenvalue2\n");

	n1_up = n2_up = 1;
	n1_down = n2_down = target_n/4;
	mu = U/2;

	//while(target_n < n1_up*4) {
		sum1 = 0;
		sum2 = 0;

		for(x=0; x<itv; x++) {
			kx = -pi + (2*pi*x/(double)itv);

			for(y=0; y<itv; y++) {
				ky = -pi + (2*pi*y/(double)itv);
				fprintf(fp, "%f\t%f\t", kx, ky);

				EigenCal(w, v, kx, ky, U, target_n, n1_down, n2_down);
				fprintf(fp, "%f\t%f\n", w[0], w[1]);
				printf("%f\t%f\n", w[0], w[1]);

				/*
				for(i=0; i<N; i++) {
					if(w[i] - EU(U, target_n, n1_up, n1_down, n2_up, n2_down) < mu) {
						for(j=0; j<LDVR*N; j+=LDVR) sum1 += pow(v[j], 2);
						for(j=1; j<LDVR*N; j+=LDVR) sum2 += pow(v[j], 2);
					}
				}
				*/
			}
		}
		n1_up = sum1/(itv*itv);
		n2_up = sum2/(itv*itv);
		mu -= 0.01;
	//}

	fclose(fp);

	return mu;
}

double DoubleMCal(double target_n, double U, double mu) { // Double cell m calculator
	int i, j, k, x, y;
	double n1_up, n2_up, n1_down, n2_down, kx, ky, sum1, sum2, m;

	double w[N], v[LDVR*N];

	n1_up = n2_up = 1;	
	n1_down = 0.1;
	n2_down = 0.3;

	for(k=0; k<128; k++) {
		sum1 = 0;
		sum2 = 0;

		for(x=0; x<itv; x++) {
			kx = -pi + (2*pi*x/(double)itv);
			for(y=0; y<itv; y++) {
				ky = -pi + (2*pi*y/(double)itv);

				EigenCal(w, v, kx, ky, U, target_n, n1_down, n2_down);

				for(i=0; i<N; i++) {
					if(w[i] - EU(U, target_n, n1_up, n1_down, n2_up, n2_down) < mu) {
						for(j=0; j<LDVR*N; j+=LDVR) sum1 += pow(v[j], 2);
						for(j=1; j<LDVR*N; j+=LDVR) sum2 += pow(v[j], 2);
					}
				}
			}
		}
		n1_up = sum1/(itv*itv);
		n1_down = target_n/2 - n1_up;

		n2_up = sum2/(itv*itv);
		n2_down = target_n/2 - n2_up;

		m = ((n1_up+n2_up) - (n1_down+n2_down))/4;
	}

	return m;
}

int main() {
	FILE *fp;

	double target_n, U, mu, m;
	//double time;

	fp = fopen("data/afmpd.txt", "w");
	fprintf(fp, "target_n\tt/U\n");

	//printf("target_n\tt/U\telapsed time(s)\n");
	
	for(target_n=1.0; target_n>0.9; target_n-=0.1) {
		//time = clock();

		for(U=0; U<1; U+=1) {
			mu = DoubleMuCal(target_n, U);
			//m = DoubleMCal(target_n, U, mu);
			//printf("%f\t%f\t%f\n", target_n, U, m);

			//if(fabs(m) > 1e-1) break;
		}
		//printf("%.3f\t%f\t%.3f\n", target_n, 1/U, (clock()-time)*0.000001);
		//fprintf(fp, "%f\t%f\n", target_n, 1/U);
	}

	fclose(fp);

	return 0;
}
