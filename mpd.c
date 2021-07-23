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

typedef struct _HubbardModel{
	double n;
	double m;
	double mu;
	double energy;
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

void FMEigenCal(_HM *hm, double u, double k1, double k2) { // FM eigenproblem calculator
	lapack_int ln = LN, lda = LDA, ldvl = LDVL, ldvr = LDVR, info;
	lapack_complex_double *vl, *work;
	double rwork[2*LN];
	double a0real, a1real, a2real, a1imag, a2imag, a3real;

	a0real = u*(hm->n/4 - hm->m); // u*d1
	a1real = a2real = -(1+cos(k1+k2)+cos(k1)+cos(k2));
	a1imag = -(sin(k1+k2)+sin(k1)+sin(k2));
	a2imag = sin(k1+k2)+sin(k1)+sin(k2);
	a3real = u*(hm->n/4 - hm->m); // u*d2

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

void FMMuCal(_HM *hm, double u, double mu_old, double n_target) { // FM mu calculator
	int i, j, ln;
	double k1, k2, u1_sum, n;

	n = 0;
	hm->n = n_target;
	hm->m = 0;
	hm->mu = mu_old;

	while(n < n_target) {
		u1_sum = 0;

		for(i=0; i<itv; i++) {
			k1 = -M_PI + 2*M_PI*i/(double)itv;
			for(j=0; j<itv; j++) {
				k2 = -M_PI + 2*M_PI*j/(double)itv;

				FMEigenCal(hm, u, k1, k2);
				for(ln=0; ln<LN; ln++) {
					if(creal(hm->w[ln]) - EU(u, hm->n, hm->m) < hm->mu) {
						u1_sum += pow(creal(hm->v[ln*2]), 2) + pow(cimag(hm->v[ln*2]), 2); 
					}
				}

			}
		}
		n = 4*u1_sum/(itv*itv);
		hm->mu += 0.01;
	}
	hm->n = n;
}

void FMMCal(_HM *hm, double u) { // FM m calculator
	int i, j, itr, ln;
	double k1, k2, u1_sum;

	hm->m = 0.1;

	for(itr=0; itr<128; itr++) {
		u1_sum = 0;

		for(i=0; i<itv; i++) {
			k1 = -M_PI + 2*M_PI*i/(double)itv;
			for(j=0; j<itv; j++) {
				k2 = -M_PI + 2*M_PI*j/(double)itv;

				FMEigenCal(hm, u, k1, k2);
				for(ln=0; ln<LN; ln++) {
					if(creal(hm->w[ln]) - EU(u, hm->n, hm->m) < hm->mu) {
						u1_sum += pow(creal(hm->v[ln*2]), 2) + pow(cimag(hm->v[ln*2]), 2); 
					}
				}
			}
		}
		hm->m = u1_sum/(itv*itv) - hm->n/4;
	}
}

void FMEnergyCal(_HM *hm, double u) { // FM energy calculator
	int i, j, ln;
	double k1, k2, energy_u1, energy_u2, energy_d1, energy_d2;

	energy_u1 = 0;
	energy_u2 = 0;
	energy_d1 = 0;
	energy_d2 = 0;

	for(i=0; i<itv; i++) {
		k1 = -M_PI + 2*M_PI*i/(double)itv;
		for(j=0; j<itv; j++) {
			k2 = -M_PI + 2*M_PI*j/(double)itv;

			FMEigenCal(hm, u, k1, k2);
			for(ln=0; ln<LN; ln++) {
				if(creal(hm->w[ln]) - EU(u, hm->n, hm->m) < hm->mu) {
					energy_u1 += (creal(hm->w[ln]) - EU(u, hm->n, hm->m)) * (pow(creal(hm->v[ln*2]), 2) + pow(cimag(hm->v[ln*2]), 2)); 
					energy_u2 += (creal(hm->w[ln]) - EU(u, hm->n, hm->m)) * (pow(creal(hm->v[ln*2+1]), 2) + pow(cimag(hm->v[ln*2+1]), 2)); 
				}
			}

			hm->m = -hm->m;
			FMEigenCal(hm, u, k1, k2);
			for(ln=0; ln<LN; ln++) {
				if(creal(hm->w[ln]) - EU(u, hm->n, hm->m) < hm->mu) {
					energy_d1 += (creal(hm->w[ln]) - EU(u, hm->n, hm->m)) * (pow(creal(hm->v[ln*2]), 2) + pow(cimag(hm->v[ln*2]), 2)); 
					energy_d2 += (creal(hm->w[ln]) - EU(u, hm->n, hm->m)) * (pow(creal(hm->v[ln*2+1]), 2) + pow(cimag(hm->v[ln*2+1]), 2)); 
				}
			}

			hm->m = -hm->m;
		}
	}
	hm->energy = (energy_u1 + energy_u2 + energy_d1 + energy_d2)/(itv*itv);
}

void AFMEigenCal(_HM *hm, double u, double k1, double k2) { // AFM eigenproblem calculator
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

void AFMMuCal(_HM *hm, double u, double mu_old, double n_target) { // AFM mu calculator
	int i, j, ln;
	double k1, k2, u1_sum, n;

	n = 0;
	hm->n = n_target;
	hm->m = 0;
	hm->mu = mu_old;

	while(n < n_target) {
		u1_sum = 0;

		for(i=0; i<itv; i++) {
			k1 = -M_PI + 2*M_PI*i/(double)itv;
			for(j=0; j<itv; j++) {
				k2 = -M_PI + 2*M_PI*j/(double)itv;

				AFMEigenCal(hm, u, k1, k2);
				for(ln=0; ln<LN; ln++) {
					if(creal(hm->w[ln]) - EU(u, hm->n, hm->m) < hm->mu) {
						u1_sum += pow(creal(hm->v[ln*2]), 2) + pow(cimag(hm->v[ln*2]), 2); 
					}
				}
			}
		}
		n = 4*u1_sum/(itv*itv);
		hm->mu += 0.01;
	}
	hm->n = n;
}

void AFMMCal(_HM *hm, double u) { // AFM m calculator
	int i, j, itr, ln;
	double k1, k2, u1_sum;

	hm->m = 0.1;

	for(itr=0; itr<128; itr++) {
		u1_sum = 0;

		for(i=0; i<itv; i++) {
			k1 = -M_PI + 2*M_PI*i/(double)itv;
			for(j=0; j<itv; j++) {
				k2 = -M_PI + 2*M_PI*j/(double)itv;

				AFMEigenCal(hm, u, k1, k2);
				for(ln=0; ln<LN; ln++) {
					if(creal(hm->w[ln]) - EU(u, hm->n, hm->m) < hm->mu) {
						u1_sum += pow(creal(hm->v[ln*2]), 2) + pow(cimag(hm->v[ln*2]), 2); 
					}
				}
			}
		}
		hm->m = u1_sum/(itv*itv) - hm->n/4;
	}
}

void AFMEnergyCal(_HM *hm, double u) { // AFM energy calculator
	int i, j, ln;
	double k1, k2, energy_u1, energy_u2, energy_d1, energy_d2;

	energy_u1 = 0;
	energy_u2 = 0;
	energy_d1 = 0;
	energy_d2 = 0;

	for(i=0; i<itv; i++) {
		k1 = -M_PI + 2*M_PI*i/(double)itv;
		for(j=0; j<itv; j++) {
			k2 = -M_PI + 2*M_PI*j/(double)itv;

			AFMEigenCal(hm, u, k1, k2);
			for(ln=0; ln<LN; ln++) {
				if(creal(hm->w[ln]) - EU(u, hm->n, hm->m) < hm->mu) {
					energy_u1 += (creal(hm->w[ln]) - EU(u, hm->n, hm->m)) * (pow(creal(hm->v[ln*2]), 2) + pow(cimag(hm->v[ln*2]), 2)); 
					energy_u2 += (creal(hm->w[ln]) - EU(u, hm->n, hm->m)) * (pow(creal(hm->v[ln*2+1]), 2) + pow(cimag(hm->v[ln*2+1]), 2)); 
				}
			}

			hm->m = -hm->m;
			AFMEigenCal(hm, u, k1, k2);
			for(ln=0; ln<LN; ln++) {
				if(creal(hm->w[ln]) - EU(u, hm->n, hm->m) < hm->mu) {
					energy_d1 += (creal(hm->w[ln]) - EU(u, hm->n, hm->m)) * (pow(creal(hm->v[ln*2]), 2) + pow(cimag(hm->v[ln*2]), 2)); 
					energy_d2 += (creal(hm->w[ln]) - EU(u, hm->n, hm->m)) * (pow(creal(hm->v[ln*2+1]), 2) + pow(cimag(hm->v[ln*2+1]), 2)); 
				}
			}

			hm->m = -hm->m;
		}
	}
	hm->energy = (energy_u1 + energy_u2 + energy_d1 + energy_d2)/(itv*itv);
}

void FMPMGraph() { // FM/PM transition graph
	//FILE *fp;
	//double time;
	double n_target, u;
	double mu_old, u_old;
	u_old = 7;

	//fp = fopen("data/fmpm.txt", "w");
	//fprintf(fp, "FM/PM Transition\nn/2\tt/u\telapsed time(s)\n");

	//#pragma omp parallel
	//#pragma omp for private(fp, n_target, u, time)
	for(int n_target_int=14; n_target_int>0; n_target_int-=2) {
		//time = clock();
		n_target = (double)n_target_int*0.1;
		_HM *hm = malloc(sizeof(_HM));
		mu_old = -u_old/2;

		for(u=u_old; u<10; u+=0.1) {
			FMMuCal(hm, u, mu_old, n_target);
			FMMCal(hm, u);
			mu_old = hm->mu - 0.5;
			//printf("thread : %d\tn : %f\tm : %f\n", omp_get_thread_num(), hm->n, hm->m);

			if(fabs(hm->m) > tol) break;
		}
		//fprintf(fp, "%f\t%f\t%f\n", hm->n/2, 1/u, (clock()-time)*0.000001);
		u_old = u;
		free(hm);
	}

	printf("FMPMGraph is done!\n");
	//fclose(fp);
}

void AFMPMGraph() { // AFM/PM transition graph
	//FILE *fp;
	//double time;
	double n_target, u;
	double mu_old, u_old;
	u_old = 1;

	//fp = fopen("data/afmpm.txt", "w");
	//fprintf(fp, "AFM/PM Transition\nn/2\tt/u\telapsed time(s)\n");

	//#pragma omp parallel
	//#pragma omp for private(n_target, u, time)
	for(int n_target_int=20; n_target_int>10; n_target_int-=2) {
		//time = clock();
		n_target = (double)n_target_int*0.1;
		_HM *hm = malloc(sizeof(_HM));
		mu_old = -u_old;

		for(u=u_old; u<10; u+=0.1) {
			AFMMuCal(hm, u, mu_old, n_target);
			AFMMCal(hm, u);
			mu_old = hm->mu - 0.5;
			//printf("thread : %d\tn : %f\tm : %f\n", omp_get_thread_num(), hm->n, hm->m);

			if(fabs(hm->m) > tol) break;
		}
		//fprintf(fp, "%f\t%f\t%f\n", hm->n/2, 1/u, (clock()-time)*0.000001);
		u_old = u;
		free(hm);
	}

	printf("AFMPMGraph is done!\n");
	//fclose(fp);
}

void FMAFMGraph() { // FM/AFM transition graph
	//FILE *fp;
	//double time;
	double n_target, u;
	double u_old, fmmu_old, afmmu_old;
	u_old = 12;

	//fp = fopen("data/fmafm.txt", "w");
	//fprintf(fp, "FM/AFM Transition\nn/2\tt/u\telapsed time(s)\n");

	//#pragma omp parallel
	//#pragma omp for private(n_target, u, time)
	for(int n_target_int=16; n_target_int>10; n_target_int-=2) {
		//time = clock();
		n_target = (double)n_target_int*0.1;
		_HM *hm_fm = malloc(sizeof(_HM));
		_HM *hm_afm = malloc(sizeof(_HM));
		fmmu_old = -u_old/2;
		afmmu_old = -u_old;
		//printf("fm n/afm n\tt/u\tfm energy\tafm energy\tfm m\tafm m\n");

		for(u=u_old; u>5; u-=1) {
			FMMuCal(hm_fm, u, fmmu_old, n_target);
			FMMCal(hm_fm, u);
			FMEnergyCal(hm_fm, u);

			AFMMuCal(hm_afm, u, afmmu_old, n_target);
			AFMMCal(hm_afm, u);
			AFMEnergyCal(hm_afm, u);
			
			fmmu_old = hm_fm->mu - 0.5;
			afmmu_old = hm_afm->mu - 0.5;

			//printf("%f/%f\t%f\t%f\t%f\t%f\t%f\n", hm_fm->n, hm_afm->n, 1/u, hm_fm->energy, hm_afm->energy, hm_fm->m, hm_afm->m);
			//if(hm_fm->energy > hm_afm->energy) break;
			printf("%f\t%f\n", hm_fm->energy+hm_fm->mu, hm_afm->energy+hm_afm->mu);
		}
		u_old = u;
		//fprintf(fp, "%f\t%f\t%f\n", hm_fm->n/2, 1/u, (clock()-time)*0.000001);
		free(hm_fm);
		free(hm_afm);
	}

	printf("FMAFMGraph is done!\n");
	//fclose(fp);
}

void EnergyPrt() { // FM, AFM energy printor
	FILE *fp;
	char buf[128];

	int i, p;
	double k1, k2, u, n_target;
	double u_old, fmmu_old, afmmu_old;
	u_old = 13;

	for(int n_target_int=18; n_target_int>17; n_target_int-=2) {
		n_target = (double)n_target_int*0.1;
		_HM *hm_fm = malloc(sizeof(_HM));
		_HM *hm_fm_ = malloc(sizeof(_HM));
		_HM *hm_afm = malloc(sizeof(_HM));
		fmmu_old = -u_old/2;
		afmmu_old = -u_old;

		for(u=u_old; u<17; u+=2) {
			sprintf(buf, "data/fmafm_n%.1fu%.1f.txt", n_target, u);
			fp = fopen(buf, "w");
			fprintf(fp, "path\tfm energy1\tfm energy2\tfm_ energy1\tfm_ energy2\tafm energy1\tafm energy2\n");

			FMMuCal(hm_fm, u, fmmu_old, n_target);
			FMMCal(hm_fm, u);
			hm_fm_->n = hm_fm->n - 0.5;
			hm_fm_->m = -hm_fm->m - 0.5;

			AFMMuCal(hm_afm, u, afmmu_old, n_target);
			AFMMCal(hm_afm, u);

			fmmu_old = hm_fm->mu;
			afmmu_old = hm_afm->mu;

			p = 0;

			for(i=0; i<itv; i++) { // rb(0,0) ~ Mb(pi,-pi)
				k1 = M_PI*i/(double)itv;
				k2 = -M_PI*i/(double)itv;

				FMEigenCal(hm_fm, u, k1, k2);
				FMEigenCal(hm_fm_, u, k1, k2);
				AFMEigenCal(hm_afm, u, k1, k2);
				fprintf(fp, "%d\t%f\t%f\t%f\t%f\t%f\t%f\n", p, creal(hm_fm->w[0]), creal(hm_fm->w[1]), creal(hm_fm_->w[0]), creal(hm_fm_->w[1]), creal(hm_afm->w[0]), creal(hm_afm->w[1]));
				p++;
			}

			for(i=0; i<itv; i++) { // Mb(pi,-pi) ~ rb(2*pi,0)
				k1 = M_PI + M_PI*i/(double)itv;
				k2 = -M_PI + M_PI*i/(double)itv;

				FMEigenCal(hm_fm, u, k1, k2);
				FMEigenCal(hm_fm_, u, k1, k2);
				AFMEigenCal(hm_afm, u, k1, k2);
				fprintf(fp, "%d\t%f\t%f\t%f\t%f\t%f\t%f\n", p, creal(hm_fm->w[0]), creal(hm_fm->w[1]), creal(hm_fm_->w[0]), creal(hm_fm_->w[1]), creal(hm_afm->w[0]), creal(hm_afm->w[1]));
				p++;
			}

			for(i=0; i<itv; i++) { // rb(2*pi,0) ~ Xb(pi,0) ~ r(0,0)
				k1 = 2*M_PI - 2*M_PI*i/(double)itv;

				FMEigenCal(hm_fm, u, k1, k2);
				FMEigenCal(hm_fm_, u, k1, k2);
				AFMEigenCal(hm_afm, u, k1, k2);
				fprintf(fp, "%d\t%f\t%f\t%f\t%f\t%f\t%f\n", p, creal(hm_fm->w[0]), creal(hm_fm->w[1]), creal(hm_fm_->w[0]), creal(hm_fm_->w[1]), creal(hm_afm->w[0]), creal(hm_afm->w[1]));
				p++;
			}

			fclose(fp);
			free(hm_fm);
			free(hm_afm);
			printf("%s is done!\n", buf);
		}
		u_old = u;
	}
	printf("EnergyPrt is done!\n");
}

int main() {
	EigenCalOpt();
	//FMPMGraph();
	//AFMPMGraph();
	FMAFMGraph();
	//EnergyPrt();

	return 0;
}
