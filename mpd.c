// The Magnetic Phase Diagram

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define pi 3.141592653689793
#define energy_up(kx, ky, U, n, m) ((-2.0*(cos(kx)+cos(ky))) + (U*(n/2.0-m))) // spin-up energy
#define energy_down(kx, ky, U, n, m) ((-2.0*(cos(kx)+cos(ky))) + (U*(n/2.0+m))) // spin-down energy

int itv = 128; // interval

int main() {
	FILE *fp;

	int i, x, y, cnt;
	double target_n, n, n_up, U, mu, m, dm, kx, ky, time;

	fp = fopen("data/mpd.txt", "w");
	fprintf(fp, "n\tt/U\n");

	printf("%s\t%s\t%s\t%s\t%s\n", "target_n", "n", "error", "transition", "elapsed time(s)");

	for(target_n=1; target_n>0; target_n-=0.1) { // 목표 n값 설정
		time = clock();
		
		for(U=0; U<10; U+=0.1) {
			mu = U/2;
			n = 1;
			m = 0;
			
			// mu iteration
			while(target_n < n) {
				cnt = 0;

				for(x=0; x<itv; x++) {
					kx = -pi + (2*pi*x/(double)itv);
					for(y=0; y<itv; y++) {
						 ky = -pi + (2*pi*y/(double)itv);
						
						 if(energy_up(kx, ky, U, target_n, m) < mu) cnt++;
					}
				}
				n = 2.0*(double)cnt/(itv*itv);
				mu -= 0.001;
			}

			// m iteration
			dm = 1;

			while(dm != 0) {
				cnt = 0;

				for(x=0; x<itv; x++) {
					kx = -pi + (2*pi*x/(double)itv);
					for(y=0; y<itv; y++) {
						ky = -pi + (2*pi*y/(double)itv);

						if(energy_up(kx, ky, U, n, m) < mu) cnt++;
					}
				}
				dm = m - (2.0*(double)cnt/(itv*itv)-1)/2.0;
				m = (2.0*(double)cnt/(itv*itv)-1)/2.0;
			}

			if(m == -0.5) break;
		}
		printf("%.1f\t%f\t%f\t%.1f\t%.3f\n", target_n, n, target_n-n, U, (clock()-time)*0.000001);
		fprintf(fp, "%f\t%f\n", n, 1.0/U);
	}	

	fclose(fp);

	return 0;
}
