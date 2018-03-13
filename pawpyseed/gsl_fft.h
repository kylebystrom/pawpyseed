#ifndef FFT_H
#define FFT_H

double* fft_calloc(int num_items, int item_size) {
	return (double*) mkl_calloc(2*num_items, item_size, 64);
}

double complex* fft_mult(int i, MKL_Complex16* x, double complex y) {
	return (x[2*i] + I * x[2*i+1]) * y;
}

void fft3d(MKL_Complex16* x, int* G_bounds, double* lattice,
	double* kpt, int* Gs, float complex* Cs, int num_waves, int* fftg);

#endif