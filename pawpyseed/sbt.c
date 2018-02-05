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

/*
The following routines are based on the Fortran program NumSBT written by J. Talman.
The algorithm performs a spherical Bessel transform in O(NlnN) time. If you adapt
this code for any purpose, please cite:
Talman, J. Computer Physics Communications 2009, 180, 332 –338.
The code is distributed under the Standard CPC license.
*/

typedef struct sbt_setup {
	double kmin;
	double kappamin;
	double rmin;
	double rhomin;
	double rho; //drho == dkappa
	double dt;
	double N;
	double complex** mult_table;
} sbt_desciptor_t;

sbt_desciptor_t* spherical_bessel_transform_setup(double encut, double enbuf, int lmax, int N, double* r) {

	sbt_desciptor_t* descriptor = (sbt_desciptor_t*) malloc(sizeof(sbt_desciptor_t));

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
		
		phi2 = -atan(sinh(PI*t/1)/cosh(PI*t/2));
		phi = phi1 + phi2 + phi3;
		mult_table[0][i] = pow(PI/2, 0.5) * cexp(I * phi) / N;
		if (it == 1) mult_table[0][it] = 0.5 * mult_table[0][it];
		phi = -phi2-atan(2*t);
		mult_table[1][i] = cexp(2 * I * phi) * mult_table[0][i];
		for (int l = 2; l <= lmax) {
			phi = -atan(2*t/(2*l+1));
			mult_table[i][l] = exp(2* I * phi) * mult_table[i][l-2];
		}
	}

	descriptor->kmin = kmin;
	descriptor->kappamin = kappamin;
	descriptor->rmin = rmin;
	descriptor->rhomin = rhomin;
	descriptor->drho = drho;
	descriptor->dt = dt;
	descriptor->N = N;
	descriptor->mult_table = mult_table;
	return descriptor;
}

double complex* wave_spherical_bessel_transform(sbt_desciptor_t* d,
	double* r, double* f, double* ks, int l) {

	double kmin = d->kmin;
	double kappamin = d->kappamin;
	double rmin = d->rmin;
	double rhomin = d->rhomin;
	double drho = d->drho;
	double dt = d->dt;
	double N = d->N;
	double complex** M = d->mult_table;

	MKL_Complex16* x = mkl_malloc(N*sizeof(MKL_Complex16), 64);

	DFTI_DESCRIPTOR_HANDLE handle = 0;
	MKL_LONG dim = 1;
	MKL_LONG length = N;
	MKL_LONG status = 0;

	double complex* vals = malloc(N * sizeof(double complex));

	status = DftiCreateDescriptor(&handle, DFTI_DOUBLE, DFTI_COMPLEX, dim, length);
	status = DftiCommitDescriptor(handle);
	status = DftiComputeBackward(handle, x);

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
		rp = x[n].real * creal(M[l][n]) - x[n].imag * cimag(M[l][n]);
		ip = x[n].imag * creal(M[l][n]) + x[n].real * cimag(M[l][n]);
		phase = I * (kappamin + rhomin) * n * dt;
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
		ks[p] = kp;
		vals[p] = x[p].real + I * x[p].imag;
		vals[p] *= 2 / N / pow(kp, 1.5);
	}
	return vals;
}
