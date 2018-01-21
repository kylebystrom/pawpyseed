#ifndef FFT_H
#define FFT_H
#include <mkl.h>
#include <mkl_types.h>

void trilinear_interpolate_values(MKL_Complex16* x, double* frac, int* fftg, double complex* values);

void fft3d(MKL_Complex16* x, int* G_bounds, double* lattice,
	double* kpt, int* Gs, float complex* Cs, int num_waves, int* fftg);

#endif
