#ifndef FFT_H
#define FFT_H
#include <mkl.h>
#include <mkl_types.h>

#define fft_complex MKL_Complex16

void trilinear_interpolate_values(MKL_Complex16* x, double* frac, int* fftg, double complex* values);

MKL_Complex16* fft_calloc(int num_items, int item_size);

double complex fft_mult(int i, MKL_Complex16* x, double complex y);

void fft3d(MKL_Complex16* x, int* G_bounds, double* lattice,
	double* kpt, int* Gs, float complex* Cs, int num_waves, int* fftg);

#endif
