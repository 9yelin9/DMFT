// LAPACK dgeev test

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include </usr/include/triqs/arrays/blas_lapack/f77/lapack.h>

#define N 2
#define LDA N
#define LDVL N
#define LDVR N

void print_eigenvalues( char* desc, int n, double* wr, double* wi ) {
        int j;
        printf( "\n %s\n", desc );
   for( j = 0; j < n; j++ ) {
      if( wi[j] == (double)0.0 ) {
         printf( " %6.2f", wr[j] );
      } else {
         printf( " (%6.2f,%6.2f)", wr[j], wi[j] );
      }
   }
   printf( "\n" );
}

void print_eigenvectors( char* desc, int n, double* wi, double* v, int ldv ) {
        int i, j;
        printf( "\n %s\n", desc );
   for( i = 0; i < n; i++ ) {
      j = 0;
      while( j < n ) {
         if( wi[j] == (double)0.0 ) {
            printf( " %6.2f", v[i+j*ldv] );
            j++;
         } else {
            printf( " (%6.2f,%6.2f)", v[i+j*ldv], v[i+(j+1)*ldv] );
            printf( " (%6.2f,%6.2f)", v[i+j*ldv], -v[i+(j+1)*ldv] );
            j += 2;
         }
      }
      printf( "\n" );
   }
}

int main() {
	int n = N, lda = LDA, ldvl = LDVL, ldvr = LDVR, info, lwork;
	double wkopt;
	double *work;

	double wr[N], wi[N], vl[LDVL*N], vr[LDVR*N];
        double a[LDA*N] = {
			1, 1,
			1, 1
        };
	
	lwork = -1;
	dgeev_("Vectors", "Vectors", &n, a, &lda, wr, wi, vl, &ldvl, vr, &ldvr, &wkopt, &lwork, &info );
	
	lwork = (int)wkopt;
	work = (double*)malloc(lwork*sizeof(double));
	dgeev_("Vectors", "Vectors", &n, a, &lda, wr, wi, vl, &ldvl, vr, &ldvr, work, &lwork, &info );
	
	if(info > 0) {
		printf("FAIL");
		exit(1);
	}

        /* Print eigenvalues */
        print_eigenvalues( "Eigenvalues", n, wr, wi );
        /* Print left eigenvectors */
        print_eigenvectors( "Left eigenvectors", n, wi, vl, ldvl );
        /* Print right eigenvectors */
        print_eigenvectors( "Right eigenvectors", n, wi, vr, ldvr );

	return 0;
}
