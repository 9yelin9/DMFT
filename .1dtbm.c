// 1D Tight-Binding Model

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define pi 3.141592653689793
#define energy_up(k, U, m) ((-2 * cos(k)) + (U * (0.5 - m))) // spin-up energy

int save_dispersion(double U, double m, int iter){
	FILE *fp;

	char fname[1024];
	sprintf(fname, "data/init_u%.2lf_%dth.txt", U, iter);
	fp = fopen(fname, "w");

	int i;
	double k;
	for(i=0; i<=128; i++){
		k = -pi + 2*pi*i/128.;
		fprintf(fp, "%lf\t%lf\n", k, energy_up(k, U, m));
	}

	fclose(fp);

	return 0;
}

int main() {
	int i;
	double n_up, U, mu, m, k;

	printf("U\tm\n");

	n_up = 0.5;
	for(U = 7; U > 3; U-= 0.1) {
		

		//n_up = 0.5; // n_down = 1 - n_up
		// 초기값
		// 1. m = 0 일 경우
		// 2. 처음부터 벌어져 있을 경우
		// 3. 애매하게 시작했으면 수렴해야되는데
		
		mu = U / 2;
		m = ((2 * n_up) - 1) / 2;
		//save_dispersion(U, m, 0);
		for(i = 0; i < 100; i++) {
			//if( i%10 == 0 ) save_dispersion(U, m, i);

			m = ((2 * n_up) - 1) / 2;
	
			for(k = 0; k < pi; k += 0.0001) {
				if(energy_up(k, U, m) > mu) {
					n_up = (2 * k) / (2 * pi);
					break;
				}
				n_up = 1.;
			}
		}
		printf("%.1lf\t%f\n", U, m);
	}

	printf("n_up = %f\n", n_up);

	return 0;
}
