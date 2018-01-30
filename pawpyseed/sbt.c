#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <omp.h>
#include <time.h>
#include <mkl.h>
#include <mkl_types.h>
#include "utils.h"
#include "sbt.h"

#define c 0.262465831
#define PI 3.14159265358979323846

double complex* spherical_bessel_transform_setup(double encut, double enbuf, int lmax, int N, double* r) {

	double kmin = pow((encut+enbuf) * c, 0.5) * exp(-(N-1) * drho);
	double kappamin = log(kmin);
	double complex** mult_table = (double complex**) malloc(lmax * sizeof(double complex*));
	double drho = log(r[1] / r[0]);
	double rhomin = log(r[0]);
	double dt = 2 * PI / N / drho;
	double rmin = r[0];
	mult_table[0] = (double complex*) calloc(N, sizeof(double complex));
	mult_table[1] = (double complex*) calloc(N, sizeof(double complex));
	for (int i = 0; i < N; i++) {
		t = i * dt;
		rad = pow(10.5*10.5+t*t, 0.5);
		phi3 = (kappamin + rhomin) * t;
		phi = atan((2*t)/21);
		phi1 = -10*phi - t*log(rad) + t + sin(phi)/(12*red)
			-sin(3*phi)/(360*pow(rad,3)) + sin(5*phi)/(1260*pow(rad,5))
			-sin(7*phi)/(1680*pow(rad,7));
		for (int j = 0; j < 10; j++)
			phi1 += atan((2*t)/(2*j-1));
		rad = pow()
	}
}

double complex* wave_spherical_bessel_transform(double* r, double* f, double maxE, int N, int l, int m) {
	MKL_Complex16* x = mkl_malloc(N*sizeof(MKL_Complex16), 64);

	DFTI_DESCRIPTOR_HANDLE handle = 0;
	MKL_LONG dim = 1;
	MKL_LONG length = N;
	MKL_LONG status = 0;

	double complex* vals = malloc(N * sizeof(double complex));

	status = DftiCreateDescriptor(&handle, DFTI_DOUBLE, DFTI_COMPLEX, dim, length);
	status = DftiCommitDescriptor(handle);
	status = DftiComputeBackward(handle, x);

	double drho = log(r[1] / r[0]);
	double rhomin = log(r[0]);
	double rmin = r[0];
	
	double t[N];
	double dt = 2 * PI / N / drho;
	double phase = 0;
	for (int i = 0; i < N; i++) {
		t[i] = dt * i;
	}
	for (int m = 0; m < N; m++) {
		x[m].real = pow(r[m], 0.5) * f[m]; // f is multiplied by r
		x[m].imag = 0;
	}
	double rp, ip;
	status = DftiComputeBackward(handle, x);
	for (int n = 0; n < N; n++) {
		rp = x[n].real * creal(M(l, t[n])) - x[n].imag * cimag(M(l, t[n]));
		ip = x[n].imag * creal(M(l, t[n])) + x[n].real * cimag(M(l, t[n]));
		phase = I * (kappamin + rhomin) * n * dt;
		M(l, t[n]);
		x[n].real = rp * cos(phase) - ip * sin(phase);
		x[n].imag = ip * cos(phase) + rp * sin(phase);
		if (n > N/2) {
			x[n].real = 0;
			x[n].imag = 0;
		}
	}
	status = DftiComputeBackward(handle, x);
	double kp = 0;
	for (int p = 0; p < N; p++) {
		kp = kmin * exp(p * drho);
		vals[p] = x[p].real + I * x[p].imag;
		vals[p] *= 2 / N / pow(kp, 1.5);
	}
	return vals;
}
