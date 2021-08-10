#include <stdio.h>
#include <math.h>

#define _USE_MATH_DEFINES

int itv = 10;

void plus(double *up, double *dn) {
	*up += 1;
	*dn += 1;
}

void test(double *up, double *dn, double *arr) {
	plus(up, dn);

	arr[0] = 1;
	arr[1] = 1;
}

int main() {
	int i = 0 , j = 0;
	double k1 = -M_PI + 2*M_PI*i/(double)itv;
	double k2 = -M_PI + 2*M_PI*j/(double)itv;

	for(i=0; i<itv*itv; i++) {
		printf("(%d,%d)\t", i/itv, i%itv);
		if(i%itv == itv-1) printf("\n");
	}

	double up = 0, dn = 0;
	double arr[2] = {0, };

	printf("b : %f, %f\n", up, dn);
	test(&up, &dn, arr);
	printf("a : %f, %f\n", up, dn);

	char a = 'A', b = 'B';
	printf("result : %c, %c\n", a, b);

	return 0;
}
