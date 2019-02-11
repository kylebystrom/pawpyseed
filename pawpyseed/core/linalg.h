/** \file
Linear algebra routines performed by interacing with the Intel Math Kernel Library
*/

#ifndef FFT_H
#define FFT_H
#include <mkl.h>
#include <mkl_types.h>

/**
Uses the 3D fast fourier transform to calculate the wavefunction
defined by plane-wave coefficients Cs in real space. These 
values get stored in x. The fast index is z (i.e. the third
lattice direction) for storage and computation.
*/
void fft3d(double complex* x, int* G_bounds, double* lattice,
	double* kpt, int* Gs, float complex* Cs, int num_waves, int* fftg);

void fwd_fft3d(double complex* x, int* G_bounds, double* lattice,
	double* kpt, int* Gs, float complex* Cs, int num_waves, int* fftg);

#endif
