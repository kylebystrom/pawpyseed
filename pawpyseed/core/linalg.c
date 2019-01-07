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

	//double test_total = 0;
	for (int w = 0; w < num_waves; w++) {
		int g1 = Gs[3*w]-G_bounds[0], g2 = Gs[3*w+1]-G_bounds[2], g3 = Gs[3*w+2]-G_bounds[4];
		x[g1*fftg[1]*fftg[2] + g2*fftg[2] + g3] = Cs[w];
		//test_total += cabs(Cs[w]) * cabs(Cs[w]);
	}

	MKL_LONG status1 = DftiCreateDescriptor(&handle, DFTI_DOUBLE, DFTI_COMPLEX, dim, length);
	MKL_LONG status2 = DftiCommitDescriptor(handle);
	MKL_LONG status3 = DftiComputeBackward(handle, x);
	//printf("%s\n%s\n%s\n", DftiErrorMessage(status1), DftiErrorMessage(status2), DftiErrorMessage(status3));

	//double kmins[3] = {G_bounds[0] + kpt[0], G_bounds[2] + kpt[1], G_bounds[4] + kpt[2]};
	double kmins[3] = {G_bounds[0], G_bounds[2], G_bounds[4]};
	double dv = determinant(lattice) / fftg[0] / fftg[1] / fftg[2];
	double inv_sqrt_vol = pow(determinant(lattice), -0.5);

	double frac[3] = {0,0,0};
	double kdotr = 0;

	double total = 0;
	double rp, ip;
	for (int i = 0; i < fftg[0]; i++) {
		for (int j = 0; j < fftg[1]; j++) {
			for (int k = 0; k  < fftg[2]; k++) {
				frac[0] = ((double)i)/fftg[0];
				frac[1] = ((double)j)/fftg[1];
				frac[2] = ((double)k)/fftg[2];
				kdotr = 2 * PI * dot(kmins, frac);
				x[i*fftg[1]*fftg[2] + j*fftg[2] + k] *= inv_sqrt_vol * cexp(I * kdotr);
				total += (pow(creal(x[i*fftg[1]*fftg[2] + j*fftg[2] + k]), 2)
					+ pow(cimag(x[i*fftg[1]*fftg[2] + j*fftg[2] + k]), 2)) * dv;
			}
		}
	}
	
	DftiFreeDescriptor(&handle);
}

double complex* fft_calloc(int num_items, int item_size) {
        return (double complex*) mkl_calloc(num_items, item_size, 64);
}


