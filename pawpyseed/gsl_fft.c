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
		x[g1*fftg[1]*fftg[2] + g2*fftg[2] + g3].real = creal(Cs[w]);
		x[g1*fftg[1]*fftg[2] + g2*fftg[2] + g3].imag = cimag(Cs[w]);
		//test_total += cabs(Cs[w]) * cabs(Cs[w]);
	}

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
				rp = x[i*fftg[1]*fftg[2] + j*fftg[2] + k].real * inv_sqrt_vol;
				ip = x[i*fftg[1]*fftg[2] + j*fftg[2] + k].imag * inv_sqrt_vol;
				x[i*fftg[1]*fftg[2] + j*fftg[2] + k].real = rp * cos(kdotr) - ip * sin(kdotr);
				x[i*fftg[1]*fftg[2] + j*fftg[2] + k].imag = ip * cos(kdotr) + rp * sin(kdotr);
				total += (pow(x[i*fftg[1]*fftg[2] + j*fftg[2] + k].real, 2)
					+ pow(x[i*fftg[1]*fftg[2] + j*fftg[2] + k].imag, 2)) * dv;
			}
		}
	}
	
	DftiFreeDescriptor(&handle);
	//printf("total %lf\n",total);
}