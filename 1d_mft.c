// 정적 평균장 이론

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define pi 3.141592

double energy(double k, double mu);
double compute_energy(double mu);

int main() {
	double mu;

	compute_energy(0);
	compute_energy(2);
	compute_energy(1);
	compute_energy(-1);

	return 0;
}

double energy(double k, double mu) {
	double t, a, res;

	t = 1;
	a = 1;

	res = -2 * t * cos(k * a) - mu;

	return round(res * 1000000) / 1000000; // 소수점 여섯째자리까지 반올림
}

double compute_energy(double mu) {
	int n;
	double i, k, k0, k1, h, sum, res, t0;

	t0 = clock();

	printf("mu = %.1f, ", mu);

	i = 0.0000001; // k interval

	// k = mu인 k값 찾기
	for(k = 0; k < pi; k += i) {
		if(energy(k, mu) == 0) {
			k0 = -k;
			k1 = k;

			printf("k = %f\n", k);

			break;
		}
	}

	// 심프슨 공식으로 -k ~ k 범위 적분하기
	sum = energy(k0, mu) + energy(k1, mu);
	n = 0;

	for(k = k0 + i; k < k1; k += i) {
		if(n % 2 == 0) {
			sum += 4 * energy(k, mu);
		}
		else {
			sum += 2 * energy(k, mu);
		}

		n++;
	}

	h = (k1 - k0) / n;
	res = fabs((h / 3) * sum);

	printf("-k ~ k energy integral = %f, elapsed time = %.3fs\n\n", res, (clock() - t0) * 0.000001);

	return res;
}
