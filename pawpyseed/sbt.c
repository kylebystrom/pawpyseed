#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <omp.h>
#include <time.h>
#include <mkl.h>
#include <mkl_types.h>
#include "utils.h"

double complex* spherical_bessel_transform(double* r, double drho, double* f, double encut, int N) {
	MKL_Complex16* x = mkl_malloc(N*sizeof(MKL_Complex16), 64);
	double t[N];
	double dt = 2 * PI / N / drho;
	for (int i = 0; i < N; i++) {
		t[i] = dt * i;
	}
	for (int m = 0; m < N; m++) {
		x[m].real = pow(r[m], 1.5) * f[m];
		x[m].imag = 0;
	}
	fft(x, N);
	for (int n = 0; n < N/2; n++) {
		rp = x[n].real * creal(M(l, t[n])) - x[n].imag * cimag(M(l, t[n]));
		ip = x[n].imag * creal(M(l, t[n])) + x[n].real * cimag(M(l, t[n]));
		phase = I * (ketamin + rhomin) * n * dt
		M(l, t[n]);
		x.real = rp * cos(phase) - ip * sin(phase);
		x.imag = ip * cos(phase) + rp * sin(phase);
		if (n > N/2) {
			x.real = 0;
			x.imag = 0;
		}
	}
	fft(x, N);
	for (int p = 0; p < N; p++) {
		vals[p] = x[p].real + I * x[p].imag;
		vals[p] *= 2 / N / pow(k[p], 1.5);
	}
	return;
}