#ifndef GSL_FFT_H
#define GSL_FFT_H
#include <gsl/gsl_fft_complex.h>

#define fft_complex double

double* fft_calloc(int num_items, int item_size);

double complex fft_mult(int i, double* x, double complex y);

void fft3d(double* x, int* G_bounds, double* lattice,
	double* kpt, int* Gs, float complex* Cs, int num_waves, int* fftg);

#endif