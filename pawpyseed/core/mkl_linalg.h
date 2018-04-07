/** \file
Linear algebra routines performed by interacing with the Intel Math Kernel Library
*/

#ifndef FFT_H
#define FFT_H
#include <mkl.h>
#include <mkl_types.h>

#define fft_complex MKL_Complex16

/**
DEPRECATED, DO NOT USE
*/
void trilinear_interpolate_values(MKL_Complex16* x, double* frac, int* fftg, double complex* values);

/**
Utility allocation function to allocate memory for an MKL complex double.
*/
MKL_Complex16* fft_calloc(int num_items, int item_size);

/**
Utility function to multiply a complex MKL double and a complex C double.
*/
double complex fft_mult(int i, MKL_Complex16* x, double complex y);

/**
Uses the 3D fast fourier transform to calculate the wavefunction
defined by plane-wave coefficients Cs in real space. These 
values get stored in x. The fast index is z (i.e. the third
lattice direction) for storage and computation.
*/
void fft3d(MKL_Complex16* x, int* G_bounds, double* lattice,
	double* kpt, int* Gs, float complex* Cs, int num_waves, int* fftg);

#endif
