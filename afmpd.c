// The Magnetic Phase Diagram

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define pi 3.141592653689793
#define ue(U, u1, d1, u2, d2, kx, ky) (U*0.5*(u1+u2) - U*(u1*d1 + u2*d2) + 0.5*sqrt(pow(U, 2)*pow(u1-u2, 2)+4*pow(2.0*(cos(kx)+cos(ky)), 2)))
#define de(U, u1, d1, u2, d2, kx, ky) (U*0.5*(d1+d2) - U*(u1*d1 + u2*d2) - 0.5*sqrt(pow(U, 2)*pow(d1-d2, 2)+4*pow(2.0*(cos(kx)+cos(ky)), 2)))

int itv = 128; // interval

double mu_iter(double target_n, double U, double *n) {
	int x, y, ucnt, dcnt;
	double u1, mu, kx, ky;

	*n = 1;
	mu = U/2;

	while(target_n < *n) {
		ucnt = 0;
		dcnt = 0;
		u1 = *n/4;

		for(x=0; x<itv; x++) {
			kx = -pi + (2*pi*x/(double)itv);
			for(y=0; y<itv; y++) {
				 ky = -pi + (2*pi*y/(double)itv);
				
				 if(ue(U, u1, u1, u1, u1, kx, ky) < mu) ucnt++;
				 if(de(U, u1, u1, u1, u1, kx, ky) < mu) dcnt++;
			}
		}
		*n = (double)ucnt/(itv*itv) + (double)dcnt/(itv*itv);
		mu -= 0.001;
	}

	return mu;
}

double m_iter(double *n, double U, double mu) {
	int x, y, ucnt, dcnt;
	double u1, u2, m, kx, ky;

	m = U/2;
	u1 = 1;
	// u1, u2 초기값에 따라 afm, fm 이 갈릴 것
	// 에너지 교차로 인해 나타나는 1차 상전이 때문

	while(fabs(m) > 1e-16) {
		ucnt = 0;
		dcnt = 0;

		for(x=0; x<itv; x++) {
			kx = -pi + (2*pi*x/(double)itv);
			for(y=0; y<itv; y++) {
				ky = -pi + (2*pi*y/(double)itv);
				
				if(ue(U, u1, (*n/2-u1), (*n/2-u1), u1, kx, ky) < mu) ucnt++;
				if(de(U, u1, (*n/2-u1), (*n/2-u1), u1, kx, ky) < mu) dcnt++;
			}
		}
		m = 0.5*((double)ucnt/(itv*itv) - (double)dcnt/(itv*itv));
		u1 -= 0.001;
	}

	return m;
}

int main() {
	FILE *fp;

	double target_n, n, U, mu, m, time;

	fp = fopen("data/mpd_test.txt", "w");
	fprintf(fp, "n\tt/U\n");

	printf("n\tt/U\telapsed time(s)\n");
	
	for(target_n=1; target_n>0.1; target_n-=0.1) {
		time = clock();

		for(U=0; U<10; U+=1) {
			mu = mu_iter(target_n, U, &n);
			m = m_iter(&n, U, mu);

			if(fabs(m) > 1e-8) break;
		}
		printf("%.3f\t%f\t%.3f\n", n, 1/U, (clock()-time)*0.000001);
		//fprintf(fp, "%f\t%f\n", n, 1/U);
	}

	fclose(fp);

	return 0;
}
