#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <float.h>
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
Talman, J. Computer Physics Communications 2009, 180, 332-338.
The code is distributed under the Standard CPC license.
*/

sbt_descriptor_t* spherical_bessel_transform_setup(double encut, double enbuf, int lmax,
	int N, double* r, double* inpks) {

	N = 2 * N;
	setbuf(stdout, NULL);
	sbt_descriptor_t* descriptor = (sbt_descriptor_t*) malloc(sizeof(sbt_descriptor_t));

	if (lmax == 0) lmax = 1;
	double* ks = (double*) malloc(N * sizeof(double));
	double* rs = (double*) malloc(N * sizeof(double));
	double complex** mult_table = (double complex**) malloc((lmax+1) * sizeof(double complex*));
	CHECK_ALLOCATION(ks);
	CHECK_ALLOCATION(rs);
	CHECK_ALLOCATION(mult_table);
	double drho = log(r[1] / r[0]);
	double rhomin = log(r[0]);
	double dt = 2 * PI / N / drho;
	double rmin = r[0];
	double kmin = pow((encut+enbuf) * c, 0.5) * exp(-(N/2-1) * drho);
	double kappamin = log(kmin);
	for (int i = 0; i < N; i++) {
		ks[i] = kmin * exp(i*drho);
		rs[i] = rmin * exp((i-N/2)*drho);
	}
	for (int i = 0; i < N/2; i++) {
		inpks[i] = ks[i];
	}
	rhomin = log(rs[0]);
	mult_table[0] = (double complex*) calloc(N, sizeof(double complex));
	mult_table[1] = (double complex*) calloc(N, sizeof(double complex));
	for (int i = 2; i <= lmax; i++)
		mult_table[i] = (double complex*) calloc(N, sizeof(double complex));
	for (int i = 0; i <= lmax; i++)
		CHECK_ALLOCATION(mult_table[i]);
	double t=0.0, rad=0.0, phi=0.0, phi1=0.0, phi2=0.0, phi3=0.0;
	for (int i = 0; i < N; i++) {
		t = i * dt;
		rad = pow(10.5*10.5+t*t, 0.5);
		phi3 = (kappamin + rhomin) * t;
		phi = atan((2*t)/21);
		phi1 = -10*phi - t*log(rad) + t + sin(phi)/(12*rad)
			-sin(3*phi)/(360*pow(rad,3)) + sin(5*phi)/(1260*pow(rad,5))
			-sin(7*phi)/(1680*pow(rad,7));
		for (int j = 1; j <= 10; j++)
			phi1 += atan((2*t)/(2*j-1));
		
		phi2 = -atan(tanh(PI*t/2));
		phi = phi1 + phi2 + phi3;
		mult_table[0][i] = pow(PI/2, 0.5) * cexp(I * phi) / N;
		if (i == 0) mult_table[0][i] = 0.5 * mult_table[0][i];
		phi = -phi2-atan(2*t);
		mult_table[1][i] = cexp(2 * I * phi) * mult_table[0][i];
		for (int l = 1; l < lmax; l++) {
			phi = -atan(2*t/(2*l+1));
			mult_table[l+1][i] = cexp(2 * I * phi) * mult_table[l-1][i];
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
	descriptor->ks = ks;
	descriptor->rs = rs;
	descriptor->lmax = lmax;
	return descriptor;
}

double* wave_spherical_bessel_transform(sbt_descriptor_t* d, double* f, int l) {

	int N = d->N;
	double complex** M = d->mult_table;
	double* ks = d->ks;
	double* r = d->rs;
	double* fs = (double*) malloc(N*sizeof(double));
	CHECK_ALLOCATION(fs);
	double C = f[0] / pow(r[N/2], l+1);
	for (int i = 0; i < N/2; i++) {
		fs[i] = C * pow(r[i], l+1);
	}
	for (int i = N/2; i < N; i++) {
		fs[i] = f[i-N/2];
	}

	double complex* x = mkl_calloc(N, sizeof(double complex), 64);

	DFTI_DESCRIPTOR_HANDLE handle = 0;
	MKL_LONG dim = 1;
	MKL_LONG length = N;
	MKL_LONG status = 0;

	double* vals = malloc(N / 2 * sizeof(double));

	status = DftiCreateDescriptor(&handle, DFTI_DOUBLE, DFTI_COMPLEX, dim, length);
	status = DftiCommitDescriptor(handle);

	for (int m = 0; m < N; m++) {
		x[m] = pow(r[m], 0.5) * fs[m];
		// f is the radial part of the function times r,
		// so only multiply by r^0.5 instead of r^1.5
	}
	status = DftiComputeBackward(handle, x);
	CHECK_STATUS(status);
	for (int n = 0; n < N; n++) {
		x[n] *= M[l][n];
		if (n >= N/2) {
			x[n] = 0;
		}
	}
	status = DftiComputeBackward(handle, x);
	CHECK_STATUS(status);
	for (int p = 0; p < N / 2; p++) {
		vals[p] = creal(x[p]);
		vals[p] *= 2 / pow(ks[p], 1.5);
	}

	mkl_free(x);
	free(fs);
	return vals;
}

double* inverse_wave_spherical_bessel_transform(sbt_descriptor_t* d, double* f, int l) {

	//double kmin = d->kmin;
	//double kappamin = d->kappamin;
	//double rmin = d->rmin;
	//double rhomin = d->rhomin;
	//double drho = d->drho;
	int N = d->N;
	double complex** M = d->mult_table;
	//double* ks = d->ks;
	//double* r = d->rs;
	double* ks = d->rs;
	double* r = d->ks;
	double* fs = (double*) malloc(N*sizeof(double));
	CHECK_ALLOCATION(fs);
	for (int i = 0; i < N/2; i++) {
		fs[i] = f[i];
	}
	for (int i = N/2; i < N; i++) {
		fs[i] = 0;
	}

	double complex* x = mkl_malloc(N * sizeof(double complex), 64);

	DFTI_DESCRIPTOR_HANDLE handle = 0;
	MKL_LONG dim = 1;
	MKL_LONG length = N;
	MKL_LONG status = 0;

	double* vals = malloc(N / 2 * sizeof(double));

	status = DftiCreateDescriptor(&handle, DFTI_DOUBLE, DFTI_COMPLEX, dim, length);
	status = DftiCommitDescriptor(handle);

	for (int m = 0; m < N; m++) {
		x[m] = pow(r[m], 1.5) * fs[m];
	}
	
	status = DftiComputeBackward(handle, x);
	CHECK_STATUS(status);
	for (int n = 0; n < N; n++) {
		x[n] *= M[l][n];
		if (n >= N/2) {
			x[n] = 0;
		}
	}
	status = DftiComputeBackward(handle, x);
	CHECK_STATUS(status);
	for (int p = 0; p < N / 2; p++) {
		vals[p] = creal(x[p + N / 2]) / PI * 2;
		//vals[p] *= 2 / pow(ks[p], 1.5);
		vals[p] *= 2 / pow(ks[p + N / 2], 1.5);
	}

	mkl_free(x);
	free(fs);
	return vals;
}

void free_sbt_descriptor(sbt_descriptor_t* d) {
	for (int l = 0; l <= d->lmax; l++) {
		free(d->mult_table[l]);
	}
	free(d->ks);
	free(d->rs);
	free(d->mult_table);
	free(d);
}
