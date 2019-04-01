#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <omp.h>
#include <time.h>
#include <mkl.h>
#include <mkl_types.h>
#include "utils.h"
#include "linalg.h"

#define PI 3.14159265359

void fft3d(double complex* x, int* G_bounds, double* lattice,
	double* kpt, int* Gs, float complex* Cs, int num_waves, int* fftg) {

	MKL_LONG status = 0;
	DFTI_DESCRIPTOR_HANDLE handle = 0;
	MKL_LONG dim = 3;
	MKL_LONG length[3] = {fftg[0], fftg[1], fftg[2]};

	int gridsize = fftg[0] * fftg[1] * fftg[2];
	for (int w = 0; w < gridsize; w++) {
		x[w] = 0;
	}
	int g1, g2, g3;
	for (int w = 0; w < num_waves; w++) {
		g1 = (Gs[3*w+0]+fftg[0]) % fftg[0];
		g2 = (Gs[3*w+1]+fftg[1]) % fftg[1];
		g3 = (Gs[3*w+2]+fftg[2]) % fftg[2];
		x[g1*fftg[1]*fftg[2] + g2*fftg[2] + g3] = Cs[w];
	}
	double inv_sqrt_vol = pow(determinant(lattice), -0.5);

	status = DftiCreateDescriptor(&handle, DFTI_DOUBLE, DFTI_COMPLEX, dim, length);
	CHECK_STATUS(status);
	status = DftiSetValue(handle, DFTI_BACKWARD_SCALE, inv_sqrt_vol);
	CHECK_STATUS(status);
	status = DftiCommitDescriptor(handle);
	CHECK_STATUS(status);
	status = DftiComputeBackward(handle, x);
	CHECK_STATUS(status);
	DftiFreeDescriptor(&handle);
}

void fwd_fft3d(double complex* x, int* G_bounds, double* lattice,
	double* kpt, int* Gs, float complex* Cs, int num_waves, int* fftg) {

	MKL_LONG status = 0;
	DFTI_DESCRIPTOR_HANDLE handle = 0;
	MKL_LONG dim = 3;
	MKL_LONG length[3] = {fftg[0], fftg[1], fftg[2]};

	int gridsize = fftg[0] * fftg[1] * fftg[2];
	int g1, g2, g3;

	double sqrt_vol = pow(determinant(lattice), 0.5);

	status = DftiCreateDescriptor(&handle, DFTI_DOUBLE, DFTI_COMPLEX, dim, length);
	CHECK_STATUS(status);
	status = DftiSetValue(handle, DFTI_FORWARD_SCALE, sqrt_vol/fftg[0]/fftg[1]/fftg[2]);
	CHECK_STATUS(status);
	status = DftiCommitDescriptor(handle);
	CHECK_STATUS(status);
	status = DftiComputeForward(handle, x);
	CHECK_STATUS(status);

	for (int w = 0; w < num_waves; w++) {
		g1 = (Gs[3*w+0]+fftg[0]) % fftg[0];
		g2 = (Gs[3*w+1]+fftg[1]) % fftg[1];
		g3 = (Gs[3*w+2]+fftg[2]) % fftg[2];
		Cs[w] = x[g1*fftg[1]*fftg[2] + g2*fftg[2] + g3];
	}
	
	DftiFreeDescriptor(&handle);
}
