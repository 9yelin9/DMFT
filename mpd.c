// Magnetic Phase Diagram

#define _USE_MATH_DEFINES

#define LN 2
#define LDA 2 
#define LDVL 2 
#define LDVR 2 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <lapack.h>
#include <omp.h>

const int k = 100; // k interval
lapack_int lwork = -1;

typedef struct HubbardModel {
	double n;
	double m;
	double mu;
	double e;
	double eu;
} HM;

void EigenCalOpt() { // EigenCal optimizer
	lapack_int ln = LN, lda = LDA, ldvl = LDVL, ldvr = LDVR, info;
	lapack_complex_double w[LN], vl[LDVL*LN], vr[LDVR*LN], wkopt;
	double rwork[2*LN];

	lapack_complex_double a[LDA*LN] = {
		1 + 1*I, 1 + 1*I,
		1 + 1*I, 1 + 1*I
	};

	LAPACK_zgeev("Vectors", "Vectors", &ln, a, &lda, w, vl, &ldvl, vr, &ldvr, &wkopt, &lwork, rwork, &info);
	if(info != 0) {
		printf("EigenCalOpt FAIL\n");
		exit(1);
	}
	lwork = (lapack_int)creal(wkopt);
}

void EigenCal(char aorf_c, char uord_c, HM *hm, double u, double k1, double k2, lapack_complex_double *w, lapack_complex_double *v) { // FM, AFM eigenproblem calculator
	lapack_int ln = LN, lda = LDA, ldvl = LDVL, ldvr = LDVR, info;
	lapack_complex_double vl[LDVL*LN], *work;
	int aorf, uord;
	double rwork[2*LN];
	double a0re, a1re, a2re, a3re, a1im, a2im;

	aorf = (2*aorf_c - 135)/5;
	uord = (2*uord_c - 153)/17;

	a0re = u*(hm->n/4 - uord*hm->m);
	a3re = u*(hm->n/4 - aorf*uord*hm->m);
	a1re = -cos(k1+k2)-cos(k1)-cos(k2)-1;
	a1im =  sin(k1+k2)+sin(k1)+sin(k2);
	a2re =  a1re;
	a2im = -a1im;

	lapack_complex_double a[LDA*LN] = {
		a0re + 0*I,    a2re + a2im*I,
		a1re + a1im*I, a3re + 0*I
	};

	work = (lapack_complex_double*)malloc(lwork*sizeof(lapack_complex_double));
	LAPACK_zgeev("Vectors", "Vectors", &ln, a, &lda, w, vl, &ldvl, v, &ldvr, work, &lwork, rwork, &info);
	if(info != 0) {
		printf("EigenCal FAIL\n");
		exit(1);
	}
	free(work);
}

void OccCnt(char aorf, HM *hm, double u, double *n, double *m) { // Occupation counter
	lapack_complex_double w[LN], v[LDVR*LN];
	double k1, k2;
	double up1_sum = 0, up2_sum = 0, dn1_sum = 0, dn2_sum = 0;

	for(int i=0; i<k*k; i++) {
		k1 = -M_PI + 2*M_PI*(i/k)/k;
		k2 = -M_PI + 2*M_PI*(i%k)/k;

		EigenCal(aorf, 'U', hm, u, k1, k2, w, v);
		if(creal(w[0]) < hm->mu) {
			up1_sum += pow(creal(v[0]), 2) + pow(cimag(v[0]), 2);
			up2_sum += pow(creal(v[1]), 2) + pow(cimag(v[1]), 2);
		}	    
		if(creal(w[1]) < hm->mu) {
			up1_sum += pow(creal(v[2]), 2) + pow(cimag(v[2]), 2);
			up2_sum += pow(creal(v[3]), 2) + pow(cimag(v[3]), 2);
		}

		EigenCal(aorf, 'D', hm, u, k1, k2, w, v);
		if(creal(w[0]) < hm->mu) {
			dn1_sum += pow(creal(v[0]), 2) + pow(cimag(v[0]), 2);
			dn2_sum += pow(creal(v[1]), 2) + pow(cimag(v[1]), 2);
		}	
		if(creal(w[1]) < hm->mu) {
			dn1_sum += pow(creal(v[2]), 2) + pow(cimag(v[2]), 2);
			dn2_sum += pow(creal(v[3]), 2) + pow(cimag(v[3]), 2);
		}
	}
	*n = (up1_sum + up2_sum + dn1_sum + dn2_sum)/(k*k);
	*m = (up1_sum - dn1_sum)/(2*k*k);
}

void MCal(char aorf, HM *hm, double u, double n_target) { // FM, AFM m calculator
	double n, m, m_cvg[3], itv;

	hm->n = n_target;
	hm->m = n_target/8;
	hm->mu = -3.5;

	printf("#%c\tn\tm\n", aorf);
	for(int i=0; i<100; i++) {
		n = 0;
		itv = 0.1;

		while(itv > 1e-5) {
			OccCnt(aorf, hm, u, &n, &m);
			if(fabs(n - n_target) < 1e-4) break;

			if(n > n_target - itv) {
				hm->mu -= itv;
				itv *= 0.1;
			}
			hm->mu += itv;
		}
		printf("%2d\t%f\t%f\n", i, n, m);
		m_cvg[i%3] = m;
		hm->m = m;
		hm->mu = floor(hm->mu);

		//if(fabs(m) < n_target/20) break;
		if(fabs((m_cvg[0] + m_cvg[1] + m_cvg[2])/3 - m) < 1e-6) break;
	}
	hm->n = n;
}

void ECal(char aorf, HM *hm, double u, double n_target) { // FM, AFM energy calculator
	lapack_complex_double w[LN], v[LDVR*LN];
	double k1, k2;
	double up1_sum = 0, up2_sum = 0, dn1_sum = 0, dn2_sum = 0;
	double eup1_sum = 0, eup2_sum = 0, edn1_sum = 0, edn2_sum = 0;

	MCal(aorf, hm, u, n_target);

	for(int i=0; i<k*k; i++) {
		k1 = -M_PI + 2*M_PI*(i/k)/k;
		k2 = -M_PI + 2*M_PI*(i%k)/k;

		EigenCal(aorf, 'U', hm, u, k1, k2, w, v);
		if(creal(w[0]) < hm->mu) {
			up1_sum += pow(creal(v[0]), 2) + pow(cimag(v[0]), 2);
			up2_sum += pow(creal(v[1]), 2) + pow(cimag(v[1]), 2);

			eup1_sum += creal(w[0]) * (pow(creal(v[0]), 2) + pow(cimag(v[0]), 2));
			eup2_sum += creal(w[0]) * (pow(creal(v[1]), 2) + pow(cimag(v[1]), 2));
		}	    
		if(creal(w[1]) < hm->mu) {
			up1_sum += pow(creal(v[2]), 2) + pow(cimag(v[2]), 2);
			up2_sum += pow(creal(v[3]), 2) + pow(cimag(v[3]), 2);

			eup1_sum += creal(w[1]) * (pow(creal(v[2]), 2) + pow(cimag(v[2]), 2));
			eup2_sum += creal(w[1]) * (pow(creal(v[3]), 2) + pow(cimag(v[3]), 2));
		}

		EigenCal(aorf, 'D', hm, u, k1, k2, w, v);
		if(creal(w[0]) < hm->mu) {
			dn1_sum += pow(creal(v[0]), 2) + pow(cimag(v[0]), 2);
			dn2_sum += pow(creal(v[1]), 2) + pow(cimag(v[1]), 2);

			edn1_sum += creal(w[0]) * (pow(creal(v[0]), 2) + pow(cimag(v[0]), 2));
			edn2_sum += creal(w[0]) * (pow(creal(v[1]), 2) + pow(cimag(v[1]), 2));
		}	
		if(creal(w[1]) < hm->mu) {
			dn1_sum += pow(creal(v[2]), 2) + pow(cimag(v[2]), 2);
			dn2_sum += pow(creal(v[3]), 2) + pow(cimag(v[3]), 2);

			edn1_sum += creal(w[1]) * (pow(creal(v[2]), 2) + pow(cimag(v[2]), 2));
			edn2_sum += creal(w[1]) * (pow(creal(v[3]), 2) + pow(cimag(v[3]), 2));
		}
	}
	hm->e = (eup1_sum + eup2_sum + edn1_sum + edn2_sum)/(k*k);
	hm->eu = u*2*(up1_sum*dn1_sum + up2_sum*dn2_sum)/(k*k*k*k); 
}

void EPrt(char aorf, char uord, double n_target, double u) { // FM, AFM energy printor
	HM hm;
	lapack_complex_double w[LN], v[LDVR*LN];
	FILE *fp;
	char buf[128];
	int p = 0;
	double k1, k2;

	sprintf(buf, "data/%c%c_n%.1fu%.1f.txt", aorf, uord, n_target, u);
	fp = fopen(buf, "w");
	fprintf(fp, "path\tenergy1\tenergy2\tmu\n");

	MCal(aorf, &hm, u, n_target);

	for(int i=0; i<k; i++) { // rb(0, 0) ~ Mb(pi, -pi)
		k1 =  M_PI*i/(double)k;
		k2 = -M_PI*i/(double)k;

		EigenCal(aorf, uord, &hm, u, k1, k2, w, v);
		fprintf(fp, "%d\t%f\t%f\t%f\n", p, creal(w[0]), creal(w[1]), hm.mu);
		p++;
	}

	for(int i=0; i<k; i++) { // Mb(pi, -pi) ~ rb(2*pi, 0)
		k1 =  M_PI + M_PI*i/(double)k;
		k2 = -M_PI + M_PI*i/(double)k;

		EigenCal(aorf, uord, &hm, u, k1, k2, w, v);
		fprintf(fp, "%d\t%f\t%f\t%f\n", p, creal(w[0]), creal(w[1]), hm.mu);
		p++;
	}

	for(int i=0; i<k; i++) { // rb(2*pi, 0) ~ Xb(pi, 0) ~ r(0, 0)
		k1 = 2*M_PI - 2*M_PI*i/(double)k;

		EigenCal(aorf, uord, &hm, u, k1, k2, w, v);
		fprintf(fp, "%d\t%f\t%f\t%f\n", p, creal(w[0]), creal(w[1]), hm.mu);
		p++;
	}

	printf("%s is done!\n", buf);
	fclose(fp);
}

void PGraphTest(char aorf, double n_target) { // FM/PM, AFM/PM transition graph test
	HM hm;
	double u_start = 1;

	clock_t t0 = clock();
	for(double u=u_start; u<15; u+=0.1) {
		printf("#1/u = %f\n", 1/u);
		MCal(aorf, &hm, u, n_target);
		printf("\n");

		if(fabs(hm.m) > n_target/20) break;
	}
	clock_t t1 = clock();

	printf("#elapsed time(s) : %f\n\n", (double)(t1-t0)/CLOCKS_PER_SEC);
}

void FAGraphTest(double n_target) { // FM/AFM transition graph test
	HM hm_f, hm_a;

	clock_t t0 = clock();
	for(double u=5; u<10; u+=0.5) {
		printf("#1/u = %f\n", 1/u);
		ECal('F', &hm_f, u, n_target);
		printf("\n");
		ECal('A', &hm_a, u, n_target);
		printf("\n#1/u\tf e\ta e\tf eu\ta eu\tf e-eu\ta e-eu\n%f\t%f\t%f\t%f\t%f\t%f\t%f\t\n\n", 1/u, hm_f.e, hm_a.e, hm_f.eu, hm_a.eu, hm_f.e-hm_f.eu, hm_a.e-hm_a.eu);

		//if(hm_f.e < hm_a.e) break;
	}
	clock_t t1 = clock();

	printf("#elapsed time(s) : %f\n", (double)(t1-t0)/CLOCKS_PER_SEC);

}

void PGraph(char aorf) { // FM/PM, AFM/PM transition graph
	HM hm;
	FILE *fp;
	char buf[128];
	double n_target, n_start, n_stop, u, u_start;

	//sprintf(buf, "data/%cP.txt", aorf);
	fp = fopen(buf, "w");
	fprintf(fp, "#n/2\t1/u\tm\telapsed time(s)\n");

	if(aorf > 67) { // FM
		n_start = 1.4;
		n_stop  = 0.1;
		u_start = 5;
	}
	else { // AFM
		n_start = 2.0;
		n_stop  = 0.9;
		u_start = 1;
	}

	clock_t tt0 = clock();
	for(n_target=n_start; n_target>n_stop; n_target-=0.2) {

		clock_t t0 = clock();
		printf("#n_target = %.1f\n", n_target);
		for(u=u_start; u<15; u+=0.1) {
			printf("#1/u = %f\n", 1/u);
			MCal(aorf, &hm, u, n_target);
			printf("\n");
			if(fabs(hm.m) > n_target/20) break;
		}
		clock_t t1 = clock();

		fprintf(fp, "%f\t%f\t%f\t%f\n", hm.n/2, 1/u, hm.m, (double)(t1-t0)/CLOCKS_PER_SEC);
	}
	clock_t tt1 = clock();

	printf("#total elapsed time(s) : %f\n", (double)(tt1-tt0)/CLOCKS_PER_SEC);
	fclose(fp);
}

void FAGraph() { // FM/AFM transition graph
	//HM hm_f, hm_a;
	FILE *fp;
	char buf[128];

	//sprintf(buf, "data/FA.txt");
	fp = fopen(buf, "w");
	fprintf(fp, "#n/2\t1/u\tm\telapsed time(s)\n");
}

int main(int argc, char *argv[]) {
	if(argc != 2) {
		printf("Usage : %s <n_target>\n", argv[0]);
		exit(1);
	}
	double n_target = atof(argv[1]);
	printf("#./mpd %.1f\n", n_target);

	EigenCalOpt();

	// Test
	//EPrt('A', 'U', n_target, 5);
	//PGraphTest('F', n_target);
	//PGraphTest('A', n_target);
	FAGraphTest(n_target);

	// Data print
	//PGraph('F');
	//PGraph('A');
	//FAGraph();

	return 0;
}
