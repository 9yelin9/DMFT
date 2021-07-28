// Magnetic Phase Diagram

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <lapack.h>
#include <omp.h>

#define _USE_MATH_DEFINES

#define LN 2
#define LDA LN
#define LDVL LN
#define LDVR LN

#define EU(u, n, m) (u*n*2*(pow(n/4, 2) - pow(m, 2)))

const int itv = 128; // interval
double tol = 1e-1; // tolerance
lapack_int lwork = 1;

typedef struct _HubbardModel {
	double n;
	double m;
	double mu;
	double e;

	// 구조체에서 빼기(함수 내에서 선언)
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
	LAPACK_zgeev("V", "V", &ln, a, &lda, 0, 0, &ldvl, 0, &ldvr, &wkopt, &lwork, rwork, &info);
	lwork = (lapack_int)creal(wkopt);
}

void EigenCal(char aorf, _HM *hm, double u, double k1, double k2) { // FM, AFM eigenproblem calculator
	lapack_int ln = LN, lda = LDA, ldvl = LDVL, ldvr = LDVR, info;
	lapack_complex_double *vl, *work;
	int aorf_int;
	double rwork[2*LN];
	double a0real, a1real, a2real, a1imag, a2imag, a3real;

	aorf_int = (2*aorf - 135)/5; // 'A' = -1, 'F' = 1

	a0real = u*(hm->n/4 - hm->m); // u*d1
	a1real = a2real = -(1+cos(k1+k2)+cos(k1)+cos(k2));
	a1imag = -(sin(k1+k2)+sin(k1)+sin(k2));
	a2imag = sin(k1+k2)+sin(k1)+sin(k2);
	a3real = u*(hm->n/4 - aorf_int*hm->m); // u*d2

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

void MuCal(char aorf, _HM *hm, double u, double n_target) { // FM, AFM mu calculator
	int i, j;
	double k1, up1, n, tmp;
// w, v 함수 내에서 정의하고 각각의 cpu에 공간 할당하기
	n = 0;
	hm->n = n_target;
	hm->m = 0;
	hm->mu = -u;

	while(n < n_target) {
		up1 = 0;

		for(i=0; i<itv; i++) {
			k1 = -M_PI + 2*M_PI*i/(double)itv;
			tmp = 0;
//#pragma omp parallel for firstprivate(hm) reduction(+ : tmp)
			for(j=0; j<itv; j++) {
				double k2 = -M_PI + 2*M_PI*j/(double)itv;

				EigenCal(aorf, hm, u, k1, k2);
				if(creal(hm->w[0]) - EU(u, hm->n, hm->m) < hm->mu) {
					tmp += pow(creal(hm->v[0]), 2) + pow(cimag(hm->v[0]), 2);
				}	
				if(creal(hm->w[1]) - EU(u, hm->n, hm->m) < hm->mu) {
					tmp += pow(creal(hm->v[2]), 2) + pow(cimag(hm->v[2]), 2); 
				}
			}
			up1 += tmp;
		}
		n = 4*up1/(itv*itv);
		hm->mu += 0.1;
	}
	hm->n = n;
}

void MCal(char aorf, _HM *hm, double u) { // FM, AFM m calculator
	int i, j, itr;
	double k1, k2, up1;

	hm->m = 0.1;

	for(itr=0; itr<128; itr++) {
		up1 = 0;

		for(i=0; i<itv; i++) {
			k1 = -M_PI + 2*M_PI*i/(double)itv;
			for(j=0; j<itv; j++) {
				k2 = -M_PI + 2*M_PI*j/(double)itv;

				EigenCal(aorf, hm, u, k1, k2);
				if(creal(hm->w[0]) - EU(u, hm->n, hm->m) < hm->mu) {
					up1 += pow(creal(hm->v[0]), 2) + pow(cimag(hm->v[0]), 2);
				}	
				if(creal(hm->w[1]) - EU(u, hm->n, hm->m) < hm->mu) {
					up1 += pow(creal(hm->v[2]), 2) + pow(cimag(hm->v[2]), 2); 
				}
			}
		}
		hm->m = up1/(itv*itv) - hm->n/4;
	}
}

void ECal(char aorf, _HM *hm, double u) { // FM, AFM energy calculator
	int i, j;
	double k1, k2, e_up1, e_up2, e_dn1, e_dn2;

	e_up1 = 0;
	e_up2 = 0;
	e_dn1 = 0;
	e_dn2 = 0;

	for(i=0; i<itv; i++) {
		k1 = -M_PI + 2*M_PI*i/(double)itv;
		for(j=0; j<itv; j++) {
			k2 = -M_PI + 2*M_PI*j/(double)itv;

			EigenCal(aorf, hm, u, k1, k2);
			if(creal(hm->w[0]) - EU(u, hm->n, hm->m) < hm->mu) {
				e_up1 += (creal(hm->w[0]) + hm->mu) * (pow(creal(hm->v[0]), 2) + pow(cimag(hm->v[0]), 2));
				e_up2 += (creal(hm->w[0]) + hm->mu) * (pow(creal(hm->v[1]), 2) + pow(cimag(hm->v[1]), 2));
			}	
			if(creal(hm->w[1]) - EU(u, hm->n, hm->m) < hm->mu) {
				e_up1 += (creal(hm->w[1]) + hm->mu) * (pow(creal(hm->v[2]), 2) + pow(cimag(hm->v[2]), 2));
				e_up2 += (creal(hm->w[1]) + hm->mu) * (pow(creal(hm->v[3]), 2) + pow(cimag(hm->v[3]), 2));
			}
		}
	}

	hm->m = -hm->m;
	for(i=0; i<itv; i++) {
		k1 = -M_PI + 2*M_PI*i/(double)itv;
		for(j=0; j<itv; j++) {
			k2 = -M_PI + 2*M_PI*i/(double)itv;

			EigenCal(aorf, hm, u, k1, k2);
			if(creal(hm->w[0]) - EU(u, hm->n, hm->m) < hm->mu) {
				e_dn1 += (creal(hm->w[0]) + hm->mu) * (pow(creal(hm->v[0]), 2) + pow(cimag(hm->v[0]), 2));
				e_dn2 += (creal(hm->w[0]) + hm->mu) * (pow(creal(hm->v[1]), 2) + pow(cimag(hm->v[1]), 2));
			}	
			if(creal(hm->w[1]) - EU(u, hm->n, hm->m) < hm->mu) {
				e_dn1 += (creal(hm->w[1]) + hm->mu) * (pow(creal(hm->v[2]), 2) + pow(cimag(hm->v[2]), 2));
				e_dn2 += (creal(hm->w[1]) + hm->mu) * (pow(creal(hm->v[3]), 2) + pow(cimag(hm->v[3]), 2));
			}
		}
	}
	hm->e = (e_up1 + e_up2 + e_dn1 + e_dn2)/(itv*itv);
	hm->m = -hm->m;
}

void EnergyPrt(double n_target) { // FM, AFM energy printor
	FILE *fp;
	char buf[128];
	int i, p;
	double k1, k2, u;

	_HM *hm_fm = malloc(sizeof(_HM));
	_HM *hm_fm_ = malloc(sizeof(_HM));
	_HM *hm_afm = malloc(sizeof(_HM));

	for(u=1; u<15; u+=2) {
		sprintf(buf, "data/fmafm_n%.1fu%.1f.txt", n_target, u);
		fp = fopen(buf, "w");
		fprintf(fp, "path\tfm energy1\tfm energy2\tfm_ energy1\tfm_ energy2\tafm energy1\tafm energy2\n");

		MuCal('F', hm_fm, u, n_target);
		MCal('F', hm_fm, u);
		ECal('F', hm_fm, u);

		hm_fm_->n = hm_fm->n;
		hm_fm_->m = hm_fm->m;

		MuCal('A', hm_afm, u, n_target);
		MCal('A', hm_afm, u);
		ECal('A', hm_afm, u);

		p = 0;

		for(i=0; i<itv; i++) { // rb(0,0) ~ Mb(pi,-pi)
			k1 = M_PI*i/(double)itv;
			k2 = -M_PI*i/(double)itv;

			EigenCal('F', hm_fm, u, k1, k2);
			EigenCal('F', hm_fm_, u, k1, k2);
			EigenCal('A', hm_afm, u, k1, k2);
			fprintf(fp, "%d\t%f\t%f\t%f\t%f\t%f\t%f\n", p, creal(hm_fm->w[0]), creal(hm_fm->w[1]), creal(hm_fm_->w[0]), creal(hm_fm_->w[1]), creal(hm_afm->w[0]), creal(hm_afm->w[1]));
			p++;
		}

		for(i=0; i<itv; i++) { // Mb(pi,-pi) ~ rb(2*pi,0)
			k1 = M_PI + M_PI*i/(double)itv;
			k2 = -M_PI + M_PI*i/(double)itv;

			EigenCal('F', hm_fm, u, k1, k2);
			EigenCal('F', hm_fm_, u, k1, k2);
			EigenCal('A', hm_afm, u, k1, k2);
			fprintf(fp, "%d\t%f\t%f\t%f\t%f\t%f\t%f\n", p, creal(hm_fm->w[0]), creal(hm_fm->w[1]), creal(hm_fm_->w[0]), creal(hm_fm_->w[1]), creal(hm_afm->w[0]), creal(hm_afm->w[1]));
			p++;
		}

		for(i=0; i<itv; i++) { // rb(2*pi,0) ~ Xb(pi,0) ~ r(0,0)
			k1 = 2*M_PI - 2*M_PI*i/(double)itv;

			EigenCal('F', hm_fm, u, k1, k2);
			EigenCal('F', hm_fm_, u, k1, k2);
			EigenCal('A', hm_afm, u, k1, k2);
			fprintf(fp, "%d\t%f\t%f\t%f\t%f\t%f\t%f\n", p, creal(hm_fm->w[0]), creal(hm_fm->w[1]), creal(hm_fm_->w[0]), creal(hm_fm_->w[1]), creal(hm_afm->w[0]), creal(hm_afm->w[1]));
			p++;
		}

		fclose(fp);
		free(hm_fm);
		free(hm_afm);
		printf("%s is done!\n", buf);
	}
	printf("EnergyPrt is done!\n");
}

void FMPMGraph() { // FM/PM transition graph
	//FILE *fp;
	//double time;
	double n_target, u;

	//fp = fopen("data/fmpm.txt", "w");
	//fprintf(fp, "FM/PM Transition\nn/2\tt/u\telapsed time(s)\n");

	for(int n_target_int=20; n_target_int>0; n_target_int-=2) {
		//time = clock();
		n_target = (double)n_target_int*0.1;
		_HM *hm = malloc(sizeof(_HM));

		for(u=5; u<15; u+=1) {
			MuCal('F', hm, u, n_target);
			MCal('F', hm, u);
			printf("%f\t%f\t%f\n", hm->n/2, 1/u, hm->m);

			//if(fabs(hm->m) > tol) break;
		}
		//fprintf(fp, "%f\t%f\t%f\n", hm->n/2, 1/u, (clock()-time)*0.000001);
		free(hm);
	}

	printf("FMPMGraph is done!\n");
	//fclose(fp);
}

void AFMPMGraph() { // AFM/PM transition graph
	//FILE *fp;
	//double time;
	double n_target, u;

	//fp = fopen("data/afmpm.txt", "w");
	//fprintf(fp, "AFM/PM Transition\nn/2\tt/u\telapsed time(s)\n");

	for(int n_target_int=20; n_target_int>10; n_target_int-=2) {
		//time = clock();
		n_target = (double)n_target_int*0.1;
		_HM *hm = malloc(sizeof(_HM));

		for(u=1; u<15; u+=0.1) {
			MuCal('A', hm, u, n_target);
			MCal('A', hm, u);

			if(fabs(hm->m) > tol) break;
		}
		//fprintf(fp, "%f\t%f\t%f\n", hm->n/2, 1/u, (clock()-time)*0.000001);
		free(hm);
	}

	printf("AFMPMGraph is done!\n");
	//fclose(fp);
}

void FMAFMGraph() { // FM/AFM transition graph
	//FILE *fp;
	//double time;
	double n_target, u;

	//fp = fopen("data/fmafm.txt", "w");
	//fprintf(fp, "FM/AFM Transition\nn/2\tt/u\telapsed time(s)\n");

	for(int n_target_int=16; n_target_int>15; n_target_int-=2) {
		//time = clock();
		n_target = (double)n_target_int*0.1;
		_HM *hm_fm = malloc(sizeof(_HM));
		_HM *hm_afm = malloc(sizeof(_HM));
		printf("n/2(fm/afm)\t1/u\tfm energy\tafm energy\tfm m\tafm m\n");

		for(u=5; u<15; u+=2) {
			MuCal('F', hm_fm, u, n_target);
			MCal('F', hm_fm, u);
			ECal('F', hm_fm, u);

			MuCal('A', hm_afm, u, n_target);
			MCal('A', hm_afm, u);
			ECal('A', hm_afm, u);

			printf("%f/%f\t%f\t%f\t%f\t%f\t%f\n", hm_fm->n/2, hm_afm->n/2, 1/u, hm_fm->e, hm_afm->e, hm_fm->m, hm_afm->m);
			//if(hm_fm->energy > hm_afm->energy) break;
		}
		//fprintf(fp, "%f\t%f\t%f\n", hm_fm->n/2, 1/u, (clock()-time)*0.000001);
		free(hm_fm);
		free(hm_afm);
	}

	printf("FMAFMGraph is done!\n");
	//fclose(fp);
}



int main() {
	//omp_set_num_threads(4);
	EigenCalOpt();
	//EnergyPrt(16);
	//FMPMGraph();
	//AFMPMGraph();
	FMAFMGraph();

	return 0;
}
