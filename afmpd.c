// Magnetic Phase Diagram - AFM/PM Transition

/*
in double cell...

u1 : spin-up at site 1
u2 : spin-up at site 2
d1 : spin-down at site 1
d2 : spin-down at site 2

u1 + d1 = n/2
u2 + d2 = n/2

n = u1+d1+u2+d2

m1 = (u1-d1)/2 = u1 - n/4
m2 = (u2-d2)/2 = u2 - n/4

m = m1 + m2

u1 = d2 // AFM
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define pi 3.141592653689793
#define TBE(kx, ky) (-2.0*(cos(kx)+cos(ky))) // Tight-binding energy
#define DUE(kx, ky, U, n, m1, m2) (U*0.5*( m1+m2+n/2) + U*(2*pow(n/4, 2)-pow(m1, 2)-pow(m2, 2)) + sqrt(pow(U, 2)*pow(m1-m2, 2)+4*pow(TBE(kx, ky), 2))/2) // Double cell spin-up energy eigenvalue
#define DDE(kx, ky, U, n, m1, m2) (U*0.5*(-m1-m2+n/2) + U*(2*pow(n/4, 2)-pow(m1, 2)-pow(m2, 2)) + sqrt(pow(U, 2)*pow(m2-m1, 2)+4*pow(TBE(kx, ky), 2))/2) // Double cell spin-down energy eigenvalue

int itv = 128; // interval

double DoubleMuIter(double target_n, double U, double *n) {
	int x, y, up_cnt, down_cnt;
	double m1, m2, mu, kx, ky;

	*n = 1;
	m1 = 0;
	m2 = 0;
	mu = U+2;

	while(target_n < *n) {
		up_cnt = 0;
		down_cnt = 0;

		for(x=0; x<itv; x++) {
			kx = -pi + (2*pi*x/(double)itv);
			for(y=0; y<itv; y++) {
				 ky = -pi + (2*pi*y/(double)itv);
				
				 if(DUE(kx, ky, U, target_n, m1, m2) < mu) up_cnt++;
				 if(DDE(kx, ky, U, target_n, m1, m2) < mu) down_cnt++; 
			}
		}
		*n = (double)up_cnt/(itv*itv) + (double)down_cnt/(itv*itv);
		mu -= 0.01;
	}

	return mu;
}

double DoubleMIter(double *n, double U, double mu) {
	int i, x, y, up_cnt, down_cnt;
	double m, m1, m2, kx, ky;

	m1 = 0.1;
	m2 = 0.1;

	for(i=0; i<128; i++) {
		up_cnt = 0;
		down_cnt = 0;

		for(x=0; x<itv; x++) {
			kx = -pi + (2*pi*x/(double)itv);
			for(y=0; y<itv; y++) {
				ky = -pi + (2*pi*y/(double)itv);

				if(DUE(kx, ky, U, *n, m1, m2) < mu) up_cnt++;
				if(DDE(kx, ky, U, *n, m1, m2) < mu) down_cnt++;
			}
		}
		m = ((double)up_cnt/(itv*itv) - (double)down_cnt/(itv*itv))/2;
		m2 = m - m1;
	}

	return m;
}

int main() {
	FILE *fp;

	double target_n, n, U, mu, m, time;

	fp = fopen("data/afmpd.txt", "w");
	fprintf(fp, "n\tt/U\n");

	printf("n\tt/U\telapsed time(s)\n");
	
	for(target_n=1.0; target_n>0.1; target_n-=0.1) {
		time = clock();

		for(U=0; U<10; U+=1) {
			mu = DoubleMuIter(target_n, U, &n);
			m = DoubleMIter(&n, U, mu);

			if(fabs(m) > 1e-1) break;
		}
		printf("%.3f\t%f\t%.3f\n", n, 1/U, (clock()-time)*0.000001);
		fprintf(fp, "%f\t%f\n", n, 1/U);
	}

	fclose(fp);

	return 0;
}
