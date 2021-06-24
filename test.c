// LAPACK Test

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <lapacke.h>

int main() {
	int i;
	double A[2][2] = {{1, 1}, {1, 1}};

	dgeev('N', 'N', 2, A);

	return 0;
}
