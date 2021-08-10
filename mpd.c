// Magnetic Phase Diagram

#define _USE_MATH_DEFINES

#define LN 2
#define LDA 2 
#define LDVL 2 
#define LDVR 2 

#define EU(u, n, m) (u*2*2*(pow(n/4, 2)-pow(m, 2)))

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

void EigenCal(char aorf, char uord, HM *hm, double u, double k1, double k2, lapack_complex_double *w, lapack_complex_double *v) { // FM, AFM eigenproblem calculator
	lapack_int ln = LN, lda = LDA, ldvl = LDVL, ldvr = LDVR, info;
	lapack_complex_double vl[LDVL*LN], *work;
	int aorf_int, uord_int;
	double rwork[2*LN];
	double a0re, a1re, a2re, a3re, a1im, a2im;

	aorf_int = (2*aorf - 135)/5;  // 'A' = -1, 'F' =  1
	uord_int = (2*uord - 153)/17; // 'U' =  1, 'D' = -1

	a0re = u*(hm->n/4 - uord_int*hm->m);
	a3re = u*(hm->n/4 - aorf_int*uord_int*hm->m);
	a1re = -cos(k1+k2)-cos(k1)-cos(k2)-1;
	a1im =  sin(k1+k2)+sin(k1)+sin(k2);
	a2re =  a1re;
	a2im = -a1im;

	lapack_complex_double a[LDA*LN] = {
		a0re + 0*I,    a2re + a2im*I,
		a1re + a1im*I, a3re + 0*I
	};
	//printf("n : %f\tm : %f\n%cA : %f\t%f\t%f\t%f\n", hm->n, hm->m, uord, a0re, a1re, a2re, a3re);

	work = (lapack_complex_double*)malloc(lwork*sizeof(lapack_complex_double));
	LAPACK_zgeev("Vectors", "Vectors", &ln, a, &lda, w, vl, &ldvl, v, &ldvr, work, &lwork, rwork, &info);
	if(info != 0) {
		printf("EigenCal FAIL\n");
		exit(1);
	}
	//printf("W : %f\t%f\t\nV : %f\t%f / %f\t%f\n", creal(w0[0]), creal(w0[1]), creal(v0[0]), creal(v0[1]), creal(v0[2]), creal(v0[3]));
	free(work);
}

void OccCnt(char aorf, HM *hm, double u, double *m) { // Occupation counter
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
	hm->n = (up1_sum + up2_sum + dn1_sum + dn2_sum)/(2*k*k);
	*m = (up1_sum - dn1_sum)/(2*k*k);
	
	printf("%f\t%f\t%f\t%f\n", hm->mu, hm->n, *m, hm->m);
}

void MCal(char aorf, HM *hm, double u, double n_target) { // FM, AFM m calculator
	double m, itv;

	hm->n = 0.1;
	hm->m = 0.01;
	hm->mu = 0;

	itv = 0.1;
	while(1) {
		OccCnt(aorf, hm, u, &m);
		if(hm->n > n_target - itv) break;
		hm->mu += itv;
	}
	
	itv = 0.01;
	for(int itr=0; itr<30; itr++) {
		while(1) {
			OccCnt(aorf, hm, u, &m);
			if(hm->n > n_target) break;
			hm->mu += itv;
		}
		if(fabs(m) < 1e-4) break;
		printf("\n");
		//printf("%f\t%f\t%f\n", hm->mu, hm->n, hm->m);
		hm->m = m;
		hm->mu -= 0.1;
	}
}

void ECal(char aorf, HM *hm, double u, double n_target) { // FM, AFM energy calculator
	lapack_complex_double w[LN], v[LDVR*LN];
	double k1, k2;
	double up1_sum = 0, up2_sum = 0, dn1_sum = 0, dn2_sum = 0;
	
	MCal(aorf, hm, u, n_target);

	for(int i=0; i<k*k; i++) {
		k1 = -M_PI + 2*M_PI*(i/k)/k;
		k2 = -M_PI + 2*M_PI*(i%k)/k;

		EigenCal(aorf, 'U', hm, u, k1, k2, w, v);
		if(creal(w[0]) < hm->mu) {
			up1_sum += (creal(w[0])) * (pow(creal(v[0]), 2) + pow(cimag(v[0]), 2));
			up2_sum += (creal(w[0])) * (pow(creal(v[1]), 2) + pow(cimag(v[1]), 2));
		}	    
		if(creal(w[1]) < hm->mu) {
			up1_sum +=  (creal(w[1])) * (pow(creal(v[2]), 2) + pow(cimag(v[2]), 2));
			up2_sum +=  (creal(w[1])) * (pow(creal(v[3]), 2) + pow(cimag(v[3]), 2));
		}

		EigenCal(aorf, 'D', hm, u, k1, k2, w, v);
		if(creal(w[0]) < hm->mu) {
			dn1_sum += (creal(w[0])) * (pow(creal(v[0]), 2) + pow(cimag(v[0]), 2));
			dn2_sum += (creal(w[0])) * (pow(creal(v[1]), 2) + pow(cimag(v[1]), 2));
		}	
		if(creal(w[1]) < hm->mu) {
			dn1_sum += (creal(w[1])) * (pow(creal(v[2]), 2) + pow(cimag(v[2]), 2));
			dn2_sum += (creal(w[1])) * (pow(creal(v[3]), 2) + pow(cimag(v[3]), 2));
		}
	}
	hm->e = (up1_sum + up2_sum + dn1_sum + dn2_sum)/(4*k*k);
}

void EPrt(char aorf, char uord, double n_target, double u_start, double u_stop) { // FM, AFM energy printor
	HM hm;
	lapack_complex_double w[LN], v[LDVR*LN];
	FILE *fp;
	char buf[128];
	int p = 0;
	double k1, k2, u;

	for(u=u_start; u<u_stop; u+=2) {
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
	printf("EPrt is done!\n");
}

void PMGraph(char aorf) { // FM/PM, AFM/PM transition graph
	HM hm;
	double n_target, u;

	printf("n/2\t1/u\tm\n");
	for(int n_target_int=20; n_target_int>12; n_target_int-=2) {
		clock_t t0 = clock();
		n_target = (double)n_target_int*0.1;

		for(u=1; u<15; u+=1) {
			MCal(aorf, &hm, u, n_target);
			printf("%f\t%f\t%f\n", hm.n/2, 1/u, hm.m);

			if(fabs(hm.m) > hm.n/2) break;
		}
		clock_t t1 = clock();
		printf("elapsed time : %f\n", (double)(t1-t0)/CLOCKS_PER_SEC);
	}
}

void FMAFMGraph() { // FM/AFM transition graph
	HM hm_fm, hm_afm;
	double n_target, u;

	printf("1/u\tfm n/2\tafm n/2\tfm energy\tafm energy\tfm m\tafm m\tfm mu\tafm mu\tfm EU\tafm EU\n");
	for(int n_target_int=16; n_target_int>15; n_target_int-=2) {
		clock_t t0 = clock();
		n_target = (double)n_target_int*0.1;

		for(u=5; u<15; u+=2) {
			ECal('F', &hm_fm, u, n_target);
			ECal('A', &hm_afm, u, n_target);
			printf("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", 1/u, hm_fm.n/2, hm_afm.n/2, hm_fm.e, hm_afm.e, hm_fm.m, hm_afm.m, hm_fm.mu, hm_afm.mu, EU(u, hm_fm.n, hm_fm.m), EU(u, hm_afm.n, hm_afm.m));

			//if(hm_fm.e > hm_afm.e) break;
		}
		clock_t t1 = clock();
		printf("elapsed time : %f\n", (double)(t1-t0)/CLOCKS_PER_SEC);
	}
}

int main() {
	//omp_set_num_threads(16);
	EigenCalOpt();
	EPrt('F', 'U', 1.0, 16, 20);
	//PMGraph('F');
	//PMGraph('A');
	//FMAFMGraph();

	return 0;
}
