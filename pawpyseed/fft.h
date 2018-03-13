#ifndef FFT_H
#define FFT_H
#include <mkl.h>
#include <mkl_types.h>

void trilinear_interpolate_values(MKL_Complex16* x, double* frac, int* fftg, double complex* values);

MKL_Complex16* fft_calloc(int num_items, int item_size) {
	return (MKL_Complex16*) mkl_calloc(num_items, item_size, 64);
}

double complex* fft_mult(int i, MKL_Complex16* x, double complex y) {
	return (x[i].real + I * x[i].imag) * y;
}

void fft3d(MKL_Complex16* x, int* G_bounds, double* lattice,
	double* kpt, int* Gs, float complex* Cs, int num_waves, int* fftg);

#endif
