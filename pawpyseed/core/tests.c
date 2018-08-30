#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include "tests.h"
#include "reader.h"
#include "utils.h"
#include "linalg.h"
#include "projector.h"
#include "sbt.h"
#include <mkl.h>
#include <mkl_types.h>
#include <assert.h>

#define PI 3.14159265359

double Ylmr(int l, int m, double theta, double phi) { return creal(Ylm(l, m, theta, phi)); }

double Ylmi(int l, int m, double theta, double phi) { return cimag(Ylm(l, m, theta, phi)); }

double* get_sbtd_ks(sbt_descriptor_t* d) {

	return d->ks;
}

int fft_check(char* wavecar, double* kpt_weights, int* fftg) {
	
	setbuf(stdout, NULL);

	pswf_t* wf = read_wavefunctions("WAVECAR", kpt_weights);
	double complex* x = (double complex*) mkl_calloc(fftg[0]*fftg[1]*fftg[2],
		sizeof(double complex), 64);
	fft3d(x, wf->G_bounds, wf->lattice, wf->kpts[0]->k, wf->kpts[0]->Gs,
		wf->kpts[0]->bands[0]->Cs, wf->kpts[0]->bands[0]->num_waves, fftg);
	int* Gs = wf->kpts[0]->Gs;
	float complex* Cs = wf->kpts[0]->bands[0]->Cs;
	double inv_sqrt_vol = pow(determinant(wf->lattice), -0.5);
	double dv = determinant(wf->lattice) / fftg[0] / fftg[1] / fftg[2];
	double* kpt = wf->kpts[0]->k;
	double total1 = 0;
	double total2 = 0;
	double total3 = 0;
	for (int i =0; i < fftg[0]; i++) {
		for (int j = 0; j < fftg[1]; j++) {
			for (int k = 0; k < fftg[2]; k++) {
				double f1 = (double)i / fftg[0];
				double f2 = (double)j / fftg[1];
				double f3 = (double)k / fftg[2];
				double complex temp = 0;
				for (int w = 0; w < wf->kpts[0]->bands[0]->num_waves; w++) {
					temp += Cs[w] * cexp((f1 * (Gs[3*w]) +
							f2 * (Gs[3*w+1]) +
							f3 * (Gs[3*w+2])) * 2 * PI * I);
					if (i == 0 && j == 0 && k == 0) total3 += pow(cabs(Cs[w]), 2);
				}
				temp *= inv_sqrt_vol;
				int ind = i*fftg[1]*fftg[2]+j*fftg[2]+k;
				total1 += pow(cabs(x[ind]), 2);
				total2 += pow(cabs(temp), 2);
				assert (cabs(x[ind] - temp) < 0.00001);
			}
		}
	}
	printf("made it %lf %lf %lf\n", total1 * dv, total2 * dv, total3);
	mkl_free(x);
	return 0;
}
