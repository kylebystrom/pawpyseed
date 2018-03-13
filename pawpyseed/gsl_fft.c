#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <omp.h>
#include <time.h>
#include "utils.h"
#include "gsl_fft.h"

void fft3d(double* x, int* G_bounds, double* lattice,
	double* kpt, int* Gs, float complex* Cs, int num_waves, int* fftg) {

	int status = 0;
	int dim = 3;
	int length[3] = {fftg[0], fftg[1], fftg[2]};

	//double test_total = 0;
	for (int w = 0; w < num_waves; w++) {
		int g1 = Gs[3*w]-G_bounds[0], g2 = Gs[3*w+1]-G_bounds[2], g3 = Gs[3*w+2]-G_bounds[4];
		x[2 * (g1*fftg[1]*fftg[2] + g2*fftg[2] + g3) + 1].real = creal(Cs[w]);
		x[2 * (g1*fftg[1]*fftg[2] + g2*fftg[2] + g3) + 1].imag = cimag(Cs[w]);
		//test_total += cabs(Cs[w]) * cabs(Cs[w]);
	}

	//printf("%s\n%s\n%s\n", DftiErrorMessage(status1), DftiErrorMessage(status2), DftiErrorMessage(status3));

	//double kmins[3] = {G_bounds[0] + kpt[0], G_bounds[2] + kpt[1], G_bounds[4] + kpt[2]};
	double kmins[3] = {G_bounds[0], G_bounds[2], G_bounds[4]};
	double dv = determinant(lattice) / fftg[0] / fftg[1] / fftg[2];
	double inv_sqrt_vol = pow(determinant(lattice), -0.5);

	gsl_fft_complex_wavetable* tab =
		(gsl_fft_complex_wavetable*) gsl_fft_complex_wavetable_alloc(fftg[0]);
	gsl_fft_complex_wavetable* tab =
		(gsl_fft_complex_wavetable*) gsl_fft_complex_wavetable_alloc(fftg[1]);
	gsl_fft_complex_wavetable* tab =
		(gsl_fft_complex_wavetable*) gsl_fft_complex_wavetable_alloc(fftg[2]);

	gsl_fft_complex_workspace* wrk =
		(gsl_fft_complex_workspace*) gsl_fft_complex_workspace_alloc(fftg[0]);
	gsl_fft_complex_workspace* wrk =
		(gsl_fft_complex_workspace*) gsl_fft_complex_workspace_alloc(fftg[1]);
	gsl_fft_complex_workspace* wrk =
		(gsl_fft_complex_workspace*) gsl_fft_complex_workspace_alloc(fftg[2]);

	int status = 0;

	status = gsl_fft_complex_backward(x, fftg[1]*fftg[2], fftg[0], tab, wrk);
	for (int p = 0; p < fftg[0]; p++) {
		status = gsl_fft_complex_backward(x + p * fftg[1] * fftg[2], fftg[2], fftg[1], tabp, wrkp);
		for (int k = 0; k < fftg[1]; k++) {
			status = gsl_fft_complex_backward(x + p * fftg[1] * fftg[2] + q * fftg[2], 1, fftg[2], tabq, wrkq);
		}
	}
	gsl_fft_complex_backward(x, 1, fftg[0]*fftg[1]*fftg[2], tab, wrk);

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
				int index = i*fftg[1]*fftg[2] + j*fftg[2] + k;
				kdotr = 2 * PI * dot(kmins, frac);
				rp = x[2*index] * inv_sqrt_vol;
				ip = x[2*index+1] * inv_sqrt_vol;
				x[2*index] = rp * cos(kdotr) - ip * sin(kdotr);
				x[2*index+1] = ip * cos(kdotr) + rp * sin(kdotr);
				total += (pow(x[2*index], 2) + pow(x[2*index+1].imag, 2)) * dv;
			}
		}
	}
	
	gsl_fft_complex_wavetable_free(tab);
	gsl_fft_complex_workspace_free(wrk);
	gsl_fft_complex_wavetable_free(tabp);
	gsl_fft_complex_workspace_free(wrkp);
	gsl_fft_complex_wavetable_free(tabq);
	gsl_fft_complex_workspace_free(wrkq);
	//printf("total %lf\n",total);
}

double* fft_calloc(int num_items, int item_size) {
	return (double*) callc(2*num_items, item_size);
}

double complex fft_mult(int i, double* x, double complex y) {
	return (x[2*i] + I * x[2*i+1]) * y;
}
