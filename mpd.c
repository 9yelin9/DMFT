// Magnetic Phase Diagram - FM/PM Transition

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define pi 3.141592653689793
#define TBE(kx, ky) (-2.0*(cos(kx)+cos(ky))) // Tight-binding energy
#define SUE(kx, ky, u, n, m) (TBE(kx, ky) + (u*(n/2.0-m))) // Single cell spin-up energy eigenvalue
#define SDE(kx, ky, u, n, m) (TBE(kx, ky) + (u*(n/2.0+m))) // Single cell spin-down energy eigenvalue

int itv = 128; // interval

double SingleMuCal(double *n, double u, double n_target) { // Single cell mu calculator
	int x, y, up_cnt, down_cnt;
	double m, mu, kx, ky;

	*n = 1;
	m = 0;
	mu = u/2;

	while(*n < n_target) {
		up_cnt = 0;
		down_cnt = 0;

		for(x=0; x<itv; x++) {
			kx = -pi + (2*pi*x/(double)itv);
			for(y=0; y<itv; y++) {
				 ky = -pi + (2*pi*y/(double)itv);
				
				 if(SUE(kx, ky, u, n_target, m) < mu) up_cnt++;
				 if(SDE(kx, ky, u, n_target, m) < mu) down_cnt++; 
			}
		}
		*n = (double)up_cnt/(itv*itv) + (double)down_cnt/(itv*itv);
		mu -= 0.001;
	}

	return mu;
}

double SingleMCal(double n, double u, double mu) { // Single cell m calculator
	int i, x, y, up_cnt, down_cnt;
	double m, kx, ky;

	m = 0.1;

	for(i=0; i<32; i++) {
		up_cnt = 0;
		down_cnt = 0;

		for(x=0; x<itv; x++) {
			kx = -pi + (2*pi*x/(double)itv);
			for(y=0; y<itv; y++) {
				ky = -pi + (2*pi*y/(double)itv);

				if(SUE(kx, ky, u, n, m) < mu) up_cnt++;
				if(SDE(kx, ky, u, n, m) < mu) down_cnt++;
			}
		}
		m = ((double)up_cnt/(itv*itv) - (double)down_cnt/(itv*itv))/2;
	}

	return m;
}

int main() {
	double n_target, n, u, u_old, mu, m, time;

	//FILE *fp;
	//fp = fopen("data/mpd.txt", "w");
	//fprintf(fp, "n\tt/u\n");
	printf("n\tt/u\telapsed time(s)\n");

	u_old = 1;
	
	for(n_target=1.0; n_target>0.1; n_target-=0.1) {
		time = clock();

		for(u=u_old; u<100; u+=0.1) {
			mu = SingleMuCal(&n, u, n_target);
			m = SingleMCal(n, u, mu);

			if(fabs(m) > 1e-1) break;
		}
		u_old = u;
		printf("%.1f\t%f\t%.3f\n", n, 1/u, (clock()-time)*0.000001);
		//fprintf(fp, "%f\t%f\n", n, 1/u);
	}

	//fclose(fp);

	return 0;
}
