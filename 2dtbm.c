// 2D Tight-Binding Model

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define pi 3.141592653689793
#define energy_up(kx, ky, U, m) ((-2*(cos(kx)+cos(ky))) + (U*(0.5-m))) // spin-up energy
#define energy_down(kx, ky, U, m) ((-2*(cos(kx)+cos(ky))) + (U*(0.5+m))) // spin-down energy

int itv = 128; // interval

void save_dispersion(double U, double m, int i) {
	FILE *fp;
	char fname[1024];
	int x, y;
	double kx, ky;

	sprintf(fname, "data/U%.0f_%d.txt", U, i);
	fp = fopen(fname, "w");
	fprintf(fp, "kx\tky\tenergy\n");

	for(x=0; x<itv; x++) {
		kx = -pi + (2*pi*x/(double)itv);
		for(y=0; y<itv; y++) {
			ky = -pi + (2*pi*y/(double)itv);
			fprintf(fp, "%f\t%f\t%f\t%f\n", kx, ky, energy_up(kx, ky, U, m), energy_down(kx, ky, U, m));
		}
	}

	fclose(fp);
}

int main() {
	int i, x, y, cnt;
	double n_up, U, mu, m, kx, ky;

	n_up = 0.5;

	for(U=0; U<10; U+=1) {
		mu = U/2;
		save_dispersion(U, m, 0);
		
		for(i=0; i<128; i++) { // iteration
			m = ((2*n_up)-1) / 2;
			cnt = 0;
			
			for(x=0; x<itv; x++) {
				kx = -pi + (2*pi*x/(double)itv);
				for(y=0; y<itv; y++) {
					ky = -pi + (2*pi*y/(double)itv);

					if(energy_up(kx, ky, U, m) < mu) cnt++;
				}
			}
			n_up = (double)cnt/(itv*itv);
		}
		save_dispersion(U, m, i);
	}
	
	return 0;
}
