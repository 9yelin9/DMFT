// Magnetic Phase Diagram - FM/PM Transition

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define pi 3.141592653689793
#define TBE(kx, ky) (-2.0*(cos(kx)+cos(ky))) // Tight-binding energy
#define SUE(kx, ky, U, n, m) (TBE(kx, ky) + (U*(n/2.0-m))) // Single cell spin-up energy eigenvalue
#define SDE(kx, ky, U, n, m) (TBE(kx, ky) + (U*(n/2.0+m))) // Single cell spin-down energy eigenvalue

int itv = 128; // interval

double SingleMuCal(double target_n, double U) { // Single cell mu calculator
	int x, y, up_cnt, down_cnt;
	double n, m, mu, kx, ky;

	n = 1;
	m = 0;
	mu = U/2;

	while(target_n < n) {
		up_cnt = 0;
		down_cnt = 0;

		for(x=0; x<itv; x++) {
			kx = -pi + (2*pi*x/(double)itv);
			for(y=0; y<itv; y++) {
				 ky = -pi + (2*pi*y/(double)itv);
				
				 if(SUE(kx, ky, U, target_n, m) < mu) up_cnt++;
				 if(SDE(kx, ky, U, target_n, m) < mu) down_cnt++; 
			}
		}
		n = (double)up_cnt/(itv*itv) + (double)down_cnt/(itv*itv);
		mu -= 0.001;
	}

	return mu;
}

double SingleMCal(double target_n, double U, double mu) { // Single cell m calculator
	int i, x, y, up_cnt, down_cnt;
	double m, kx, ky;

	m = 0.1;

	for(i=0; i<128; i++) {
		up_cnt = 0;
		down_cnt = 0;

		for(x=0; x<itv; x++) {
			kx = -pi + (2*pi*x/(double)itv);
			for(y=0; y<itv; y++) {
				ky = -pi + (2*pi*y/(double)itv);

				if(SUE(kx, ky, U, target_n, m) < mu) up_cnt++;
				if(SDE(kx, ky, U, target_n, m) < mu) down_cnt++;
			}
		}
		m = ((double)up_cnt/(itv*itv) - (double)down_cnt/(itv*itv))/2;
	}

	return m;
}

int main() {
	FILE *fp;

	double target_n, U, mu, m, time;

	fp = fopen("data/mpd.txt", "w");
	fprintf(fp, "target_n\tt/U\n");

	printf("target_n\tt/U\telapsed time(s)\n");
	
	for(target_n=1.0; target_n>0.1; target_n-=0.1) {
		time = clock();

		for(U=0; U<10; U+=0.1) {
			mu = SingleMuCal(target_n, U);
			m = SingleMCal(target_n, U, mu);

			if(fabs(m) > 1e-1) break;
		}
		printf("%.1f\t%f\t%.3f\n", target_n, 1/U, (clock()-time)*0.000001);
		fprintf(fp, "%f\t%f\n", target_n, 1/U);
	}

	fclose(fp);

	return 0;
}
