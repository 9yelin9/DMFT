// Magnetic Phase Diagram - AFM/PM Transition

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <lapack.h>

#define _USE_MATH_DEFINES

#define LN 2
#define LDA LN
#define LDVL LN
#define LDVR LN

#define EU(u, n, u1, u2, d1, d2) (u*n*(u1*d1 + u2*d2))

const int itv = 128; // interval
lapack_int lwork = 1;

typedef struct _HubbardModel{
	double n;
	double mu;
	double u1;
	double u2;
	double d1;
	double d2;
	double m;
	double m1;
	double m2;
	lapack_complex_double w[LN];
	lapack_complex_double v[LDVR*LN];
} _HM;

void EigenCalOpt() { // EigenCal optimizer
	lapack_int ln = LN, lda = LDA, ldvl = LDVL, ldvr = LDVR, info;
	lapack_complex_double wkopt;
	double rwork[2*LN];

	lapack_complex_double a[LDA*LN] = {
		1+1*I, 1+1*I,
		1+1*I, 1+1*I
	};

	lwork = -1;
	LAPACK_zgeev("N", "V", &ln, a, &lda, 0, 0, &ldvl, 0, &ldvr, &wkopt, &lwork, rwork, &info);
	lwork = (int)creal(wkopt);
}

void UEigenCal(_HM *hm, double u, double kx, double ky) { // Spin-up eigenproblem calculator
	lapack_int ln = LN, lda = LDA, ldvl = LDVL, ldvr = LDVR, info;
	lapack_complex_double *vl, *work;
	double rwork[2*LN];
	double a0real, a1real, a2real, a1imag, a2imag, a3real;

	a0real = u*(hm->d1 - hm->m1);
	a1real = a2real = -1-cos(-kx-ky)-cos(-kx)-cos(-ky);
	a1imag = +sin(-kx-ky)+sin(-kx)+sin(-ky);
	a2imag = -sin(-kx-ky)-sin(-kx)-sin(-ky);
	a3real = u*(hm->d2 - hm->m2);

	lapack_complex_double a[LDA*LN] = {
		a0real+0*I, a1real+a1imag*I,
		a2real+a2imag*I, a3real+0*I
	};

	vl = (lapack_complex_double*)malloc((LDVL*LN)*sizeof(lapack_complex_double));
	work = (lapack_complex_double*)malloc(lwork*sizeof(lapack_complex_double));

	LAPACK_zgeev("N", "V", &ln, a, &lda, hm->w, vl, &ldvl, hm->v, &ldvr, work, &lwork, rwork, &info);

	if(info > 0) {
		printf("LAPACK_zgeev FAIL\n");
		exit(1);
	}

	free(vl);
	free(work);
}

void DEigenCal(_HM *hm, double u, double kx, double ky) { // Spin-down eigenproblem calculator
	lapack_int ln = LN, lda = LDA, ldvl = LDVL, ldvr = LDVR, info;
	lapack_complex_double *vl, *work;
	double rwork[2*LN];
	double a0real, a1real, a2real, a1imag, a2imag, a3real;

	a0real = u*(hm->u1 + hm->m1);
	a1real = a2real = -1-cos(-kx-ky)-cos(-kx)-cos(-ky);
	a1imag = +sin(-kx-ky)+sin(-kx)+sin(-ky);
	a2imag = -sin(-kx-ky)-sin(-kx)-sin(-ky);
	a3real = u*(hm->u2 + hm->m2);

	lapack_complex_double a[LDA*LN] = {
		a0real+0*I, a1real+a1imag*I,
		a2real+a2imag*I, a3real+0*I
	};

	vl = (lapack_complex_double*)malloc((LDVL*LN)*sizeof(lapack_complex_double));
	work = (lapack_complex_double*)malloc(lwork*sizeof(lapack_complex_double));

	LAPACK_zgeev("N", "V", &ln, a, &lda, hm->w, vl, &ldvl, hm->v, &ldvr, work, &lwork, rwork, &info);

	if(info > 0) {
		printf("LAPACK_zgeev FAIL\n");
		exit(1);
	}

	free(vl);
	free(work);
}

void DoubleMuCal(_HM *hm, double u, double n_target) { // Double cell mu calculator
	int i, x, y;
	double kx, ky, u1_sum, u2_sum, d1_sum, d2_sum;

	hm->n = 1.1;
	hm->mu = u/2;
	hm->m1 = hm->m2 = 0;
	hm->u1 = hm->u2 = hm->d1 = hm->d2 = 1;

	while(hm->n > n_target) {
		u1_sum = 0;
		u2_sum = 0;
		d1_sum = 0;
		d2_sum = 0;

		for(x=0; x<itv; x++) {
			kx = -M_PI+ (2*M_PI*x/(double)itv);
			for(y=0; y<itv; y++) {
				ky = -M_PI+ (2*M_PI*y/(double)itv);

				UEigenCal(hm, u, kx, ky);
				for(i=0; i<LN; i++) {
					if(creal(hm->w[i]) - EU(u, hm->n, hm->u1, hm->u2, hm->d1, hm->d2) < hm->mu) {
						u1_sum += pow(creal(hm->v[2*i]), 2) + pow(cimag(hm->v[2*i]), 2);
						u2_sum += pow(creal(hm->v[2*i+1]), 2) + pow(cimag(hm->v[2*i+1]), 2);
					}
				}
				DEigenCal(hm, u, kx, ky);
				for(i=0; i<LN; i++) {
					if(creal(hm->w[i]) - EU(u, hm->n, hm->u1, hm->u2, hm->d1, hm->d2) < hm->mu) {
						d1_sum += pow(creal(hm->v[2*i]), 2) + pow(cimag(hm->v[2*i]), 2);
						d2_sum += pow(creal(hm->v[2*i+1]), 2) + pow(cimag(hm->v[2*i+1]), 2);
					}
				}
			}
		}
		hm->n = (u1_sum + u2_sum + d1_sum + d2_sum)/(itv*itv);	
		hm->u1 = u1_sum/(hm->n*itv*itv);
		hm->u2 = u2_sum/(hm->n*itv*itv);
		hm->d1 = d1_sum/(hm->n*itv*itv);
		hm->d2 = d2_sum/(hm->n*itv*itv);
		
		//printf("%.1f\t%.1f\t%f\t%f\t%f\t%f\n", n_target, u, hm->u1, hm->u2, hm->d1, hm->d2);
		//printf("%.1f\t%.1f\t%f\n", n_target, u, hm->n);

		hm->mu -= 0.01;
	}
	//printf("\n");
}

void DoubleMCal(_HM *hm, double u) { // Double cell m calculator
	int itr, i, x, y;
	double kx, ky, u1_sum, d1_sum, u2_sum, d2_sum; 

	hm->m1 = 0.1;
	hm->m2 = 0.1;

	for(itr=0; itr<32; itr++) {
		u1_sum = 0;
		u2_sum = 0;
		d1_sum = 0;
		d2_sum = 0;

		for(x=0; x<itv; x++) {
			kx = -M_PI+ (2*M_PI*x/(double)itv);
			for(y=0; y<itv; y++) {
				ky = -M_PI+ (2*M_PI*y/(double)itv);

				UEigenCal(hm, u, kx, ky);
				for(i=0; i<LN; i++) {
					if(creal(hm->w[i]) - EU(u, hm->n, hm->u1, hm->u2, hm->d1, hm->d2) < hm->mu) {
						u1_sum += pow(creal(hm->v[2*i]), 2) + pow(cimag(hm->v[2*i]), 2);
						u2_sum += pow(creal(hm->v[2*i+1]), 2) + pow(cimag(hm->v[2*i+1]), 2);
					}
				}
				DEigenCal(hm, u, kx, ky);
				for(i=0; i<LN; i++) {
					if(creal(hm->w[i]) - EU(u, hm->n, hm->u1, hm->u2, hm->d1, hm->d2) < hm->mu) {
						d1_sum += pow(creal(hm->v[2*i]), 2) + pow(cimag(hm->v[2*i]), 2);
						d2_sum += pow(creal(hm->v[2*i+1]), 2) + pow(cimag(hm->v[2*i+1]), 2);
					}
				}
			}
		}
		hm->u1 = u1_sum/(hm->n*itv*itv);
		hm->u2 = u2_sum/(hm->n*itv*itv);
		hm->d1 = d1_sum/(hm->n*itv*itv);
		hm->d2 = d2_sum/(hm->n*itv*itv);

		hm->m1 = (hm->u1 - hm->d1)/2;
		hm->m2 = (hm->u2 - hm->d2)/2;
		hm->m = hm->m1 + hm->m2;
		
		//printf("%.1f\t%f\t%f\n", u, hm->n, hm->u1+hm->u2+hm->d1+hm->d2); 
		printf("%.3f\t%.1f\t%f\t%f\t%f\t%f\n", hm->n, u, hm->u1, hm->u2, hm->d1, hm->d2); 
		//printf("%.3f\t%.1f\t%f\t%f\t%f\n", hm->n, u, hm->m1, hm->m2, hm->m);
	}
	printf("\n");
}

int main() {
	double n_target, u, u_old;
	double time;

	//FILE *fp;
	//fp = fopen("data/afmpd.txt", "w");
	//fprintf(fp, "n\tt/u\n");
	//printf("n\tt/u\telapsed time(s)\n");

	_HM *hm = malloc(sizeof(_HM));
	EigenCalOpt();
	u_old = 1;

	for(n_target=1.0; n_target>0.1; n_target-=0.1) {
		time = clock();

		for(u=u_old; u<100; u+=1) {
			DoubleMuCal(hm, u, n_target);
			DoubleMCal(hm, u);

			if(fabs(hm->m) > 1e-1) break;
		}
		u_old = u;
		//printf("%.3f\t%f\t%.3f\n", hm->n, 1/u, (clock()-time)*0.000001);
		//fprintf(fp, "%f\t%f\n", hm->n, 1/u);
	}

	//fclose(fp);

	return 0;
}
