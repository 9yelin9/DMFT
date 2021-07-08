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

#define EU(u, n, m) (u*n*2*(pow(n/4, 2) - pow(m, 2)))

const int itv = 128; // interval
lapack_int lwork = 1;

typedef struct _HubbardModel{
	double n;
	double m;
	double mu;
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
	lwork = (lapack_int)creal(wkopt);
}

void EigenCal(_HM *hm, double u, double k1, double k2) { // Eigenproblem calculator
	lapack_int ln = LN, lda = LDA, ldvl = LDVL, ldvr = LDVR, info;
	lapack_complex_double *vl, *work;
	double rwork[2*LN];
	double a0real, a1real, a2real, a1imag, a2imag, a3real;

	a0real = u*(hm->n/4 - hm->m); // u*d1
	a1real = a2real = -(1+cos(k1+k2)+cos(k1)+cos(k2));
	a1imag = -(sin(k1+k2)+sin(k1)+sin(k2));
	a2imag = sin(k1+k2)+sin(k1)+sin(k2);
	a3real = u*(hm->n/4 + hm->m); // u*d2

	lapack_complex_double a[LDA*LN] = {
		a0real+0*I, a1real+a1imag*I,
		a2real+a2imag*I, a3real+0*I
	};

	vl = (lapack_complex_double*)malloc((LDVL*LN)*sizeof(lapack_complex_double));
	work = (lapack_complex_double*)malloc(lwork*sizeof(lapack_complex_double));

	LAPACK_zgeev("V", "V", &ln, a, &lda, hm->w, vl, &ldvl, hm->v, &ldvr, work, &lwork, rwork, &info);

	if(info > 0) {
		printf("LAPACK_zgeev FAIL\n");
		exit(1);
	}

	free(vl);
	free(work);
}

void DoubleMuCal(_HM *hm, double u, double mu_old, double n_target) { // Double cell mu calculator
	int i, j, ln;
	double k1, k2, u1_sum, u2_sum, n;

	n = 0;
	hm->n = n_target;
	hm->m = 0;
	hm->mu = mu_old;

	while(n < n_target) {
		u1_sum = 0;
		u2_sum = 0;

		for(i=0; i<itv; i++) {
			k1 = -M_PI+ (2*M_PI*i/(double)itv);
			for(j=0; j<itv; j++) {
				k2 = -M_PI+ (2*M_PI*j/(double)itv);

				EigenCal(hm, u, k1, k2);
				for(ln=0; ln<LN; ln++) {
					if(creal(hm->w[ln]) - EU(u, hm->n, hm->m) < hm->mu) {
						u1_sum += pow(creal(hm->v[2*ln]), 2) + pow(cimag(hm->v[2*ln]), 2);
						u2_sum += pow(creal(hm->v[2*ln+1]), 2) + pow(cimag(hm->v[2*ln+1]), 2);
					}
				}
			}
		}
		n = 2*(u1_sum + u2_sum)/(itv*itv);
		hm->mu += 0.01;
		//printf("%f\t%.1f\t%f\t%f\n", n, u, hm->mu, hm->m);
	}
	hm->n = n;
}

void DoubleMCal(_HM *hm, double u) { // Double cell m calculator
	int i, j, itr, ln;
	double k1, k2, u1_sum, u2_sum;

	hm->m = 0.1;

	for(itr=0; itr<32; itr++) {
		u1_sum = 0;
		u2_sum = 0;

		for(i=0; i<itv; i++) {
			k1 = -M_PI+ (2*M_PI*i/(double)itv);
			for(j=0; j<itv; j++) {
				k2 = -M_PI+ (2*M_PI*j/(double)itv);

				EigenCal(hm, u, k1, k2);
				for(ln=0; ln<LN; ln++) {
					if(creal(hm->w[ln]) - EU(u, hm->n, hm->m) < hm->mu) {
						u1_sum += pow(creal(hm->v[2*ln]), 2) + pow(cimag(hm->v[2*ln]), 2);
						u2_sum += pow(creal(hm->v[2*ln+1]), 2) + pow(cimag(hm->v[2*ln+1]), 2);
					}
				}
			}
		}
		hm->m = (u1_sum - u2_sum)/(2*itv*itv);
		//printf("%f\t%.1f\t%f\n", hm->n, u, hm->m);
	}
}

void EnergyCal(double n_target) { // Energy calculator
	int i, p;
	double k1, k2, u, mu_old;

	FILE *fp;
	char buf[128];

	_HM *hm = malloc(sizeof(_HM));
	EigenCalOpt();
	mu_old = -2;

	for(u=1; u<10; u+=0.1) {
		sprintf(buf, "data/afm_n%.1fu%.1f.txt", n_target, u);
		fp = fopen(buf, "w");
		fprintf(fp, "path\tenergy1\tenergy2\tmu");

		DoubleMuCal(hm, u, mu_old, n_target);
		DoubleMCal(hm, u);

		p = 0;

		// rb(0, 0) ~ Mb(pi, -pi)
		for(i=0; i<itv; i++) {
			k1 = M_PI*i/(double)itv;
			k2 = -M_PI*i/(double)itv;

			EigenCal(hm, u, k1, k2);
			fprintf(fp, "%d\t%f\t%f\t%f\n", p, creal(hm->w[0]), creal(hm->w[1]), hm->mu);
			p++;
		}

		// Mb(pi, -pi) ~ rb(2*pi, 0)
		for(i=0; i<itv; i++) {
			k1 = M_PI + M_PI*i/(double)itv;
			k2 = -M_PI + M_PI*i/(double)itv;

			EigenCal(hm, u, k1, k2);
			fprintf(fp, "%d\t%f\t%f\t%f\n", p, creal(hm->w[0]), creal(hm->w[1]), hm->mu);
			p++;
		}

		// rb(2*pi, 0) ~ Xb(pi, 0) ~ r(0, 0)
		for(i=0; i<itv; i++) {
			k1 = 2*M_PI - 2*M_PI*i/(double)itv;

			EigenCal(hm, u, k1, k2);
			fprintf(fp, "%d\t%f\t%f\t%f\n", p, creal(hm->w[0]), creal(hm->w[1]), hm->mu);
			p++;
		}

		fclose(fp);
		
		mu_old = hm->mu - 0.5;
		if(fabs(hm->m) > 1.2e-1) break;
	}

	printf("Complete EnergyCal\n");
}


int main() {
	double n_target, u, u_old, mu_old;
	double time;

	_HM *hm = malloc(sizeof(_HM));
	EigenCalOpt();
	u_old = 1;

	// Band structure
	//EnergyCal(2);

	// 

	// AFM/PM Transition
	//FILE *fp;
	//fp = fopen("data/afmpd.txt", "w");
	//fprintf(fp, "n/2\tt/u\n");
	//printf("n/2\tt/u\telapsed time(s)\n");

	for(n_target=2; n_target>1.0; n_target-=0.2) {
		time = clock();
		mu_old = -(u_old+1);

		for(u=u_old; u<10; u+=0.1) {
			DoubleMuCal(hm, u, mu_old, n_target);
			DoubleMCal(hm, u);
			mu_old = hm->mu - 0.5;

			if(fabs(hm->m) > 1.2e-1) break;
		}
		u_old = u;
		printf("%.3f\t%f\t%.3f\n", hm->n/2, 1/u, (clock()-time)*0.000001);
		//fprintf(fp, "%f\t%f\n", hm->n/2, 1/u);
	}

	//fclose(fp);

	return 0;
}
