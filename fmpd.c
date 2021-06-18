// The Magnetic Phase Diagram

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define pi 3.141592653689793
#define uenergy(kx, ky, U, n, m) ((-2.0*(cos(kx)+cos(ky))) + (U*(n/2.0-m))) // spin-up energy
#define denergy(kx, ky, U, n, m) ((-2.0*(cos(kx)+cos(ky))) + (U*(n/2.0+m))) // spin-down energy

int itv = 128; // interval

double mu_iter(double target_n, double U, double *n) {
	int x, y, ucnt, dcnt;
	double m, mu, kx, ky;

	*n = 1;
	m = 0;
	mu = U/2;

	while(target_n < *n) {
		ucnt = 0;
		dcnt = 0;

		for(x=0; x<itv; x++) {
			kx = -pi + (2*pi*x/(double)itv);
			for(y=0; y<itv; y++) {
				 ky = -pi + (2*pi*y/(double)itv);
				
				 if(uenergy(kx, ky, U, target_n, m) < mu) ucnt++;
				 if(denergy(kx, ky, U, target_n, m) < mu) dcnt++;
			}
		}
		*n = (double)ucnt/(itv*itv) + (double)dcnt/(itv*itv);
		mu -= 0.001;
	}

	return mu;
}

double m_iter(double *n, double U, double mu) {
	int i, x, y, ucnt, dcnt;
	double m, kx, ky;

	i = 0;
	m = 0.1;

	while(i != 128) {
		ucnt = 0;
		dcnt = 0;

		for(x=0; x<itv; x++) {
			kx = -pi + (2*pi*x/(double)itv);
			for(y=0; y<itv; y++) {
				ky = -pi + (2*pi*y/(double)itv);

				if(uenergy(kx, ky, U, *n, m) < mu) ucnt++;
				if(denergy(kx, ky, U, *n, m) < mu) dcnt++;
			}
		}
		m = 0.5 * ((double)ucnt/(itv*itv) - (double)dcnt/(itv*itv));
		i++;
	}

	return m;
}

int main() {
	FILE *fp;

	double target_n, n, U, mu, m, time;

	fp = fopen("data/mpd_test.txt", "w");
	fprintf(fp, "n\tt/U\n");

	printf("n\tt/U\telapsed time(s)\n");
	
	for(target_n=0.3; target_n>0.1; target_n-=0.1) {
		time = clock();

		for(U=0; U<10; U+=0.1) {
			mu = mu_iter(target_n, U, &n);
			m = m_iter(&n, U, mu);
			printf("%f\t%f\t%f\n", n, U, m);

			if(fabs(m) > 1e-1) break;
		}
		printf("%.3f\t%f\t%.3f\n", n, 1/U, (clock()-time)*0.000001);
		fprintf(fp, "%f\t%f\n", n, 1/U);
	}

	fclose(fp);

	return 0;
}
