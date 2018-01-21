#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include "tests.h"
#include "reader.h"
#include "utils.h"
#include "fft.h"
#include "projector.h"
#include <mkl.h>
#include <mkl_types.h>
#include <assert.h>

#define PI 3.14159265359

double Ylmr(int l, int m, double theta, double phi) { return creal(Ylm(l, m, theta, phi)); }

double Ylmi(int l, int m, double theta, double phi) { return cimag(Ylm(l, m, theta, phi)); }

int fft_check(char* wavecar, double* kpt_weights, int* fftg) {
	
	freopen("mytest.out", "w", stdout);
	setbuf(stdout, NULL);

	printf("made it 0");
	pswf_t* wf = read_wavefunctions("WAVECAR", kpt_weights);
	MKL_Complex16* x = (MKL_Complex16*) mkl_calloc(fftg[0]*fftg[1]*fftg[2], sizeof(MKL_Complex16), 64);
	printf("%d\n", fftg[0]);
	printf("%d\n", fftg[1]);
	printf("%d\n", fftg[2]);
	printf("made it 1");
	fft3d(x, wf->G_bounds, wf->lattice, wf->kpts[0]->k, wf->kpts[0]->Gs, wf->kpts[0]->bands[0]->Cs, wf->kpts[0]->bands[0]->num_waves, fftg);
	int* Gs = wf->kpts[0]->Gs;
	float complex* Cs = wf->kpts[0]->bands[0]->Cs;
	double inv_sqrt_vol = pow(determinant(wf->lattice), -0.5);
	double dv = determinant(wf->lattice) / fftg[0] / fftg[1] / fftg[2];
	double* kpt = wf->kpts[0]->k;
	printf("made it 2\n");
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
					temp += Cs[w] * cexp((f1 * (Gs[3*w] + kpt[0]) +
							f2 * (Gs[3*w+1] + kpt[1]) +
							f3 * (Gs[3*w+2] + kpt[2])) * 2 * PI * I);
					if (i == 0 && j == 0 && k == 0) total3 += pow(cabs(Cs[w]), 2);
				}
				temp *= inv_sqrt_vol;
				int ind = i*fftg[1]*fftg[2]+j*fftg[2]+k;
				total1 += x[ind].real * x[ind].real + x[ind].imag * x[ind].imag;
				total2 += cabs(temp) * cabs(temp);
				printf("check %lf %lf %lf %lf\n", creal(temp), cimag(temp), x[ind].real, x[ind].imag);
				assert (cabs(x[ind].real + I * x[ind].imag - temp) < 0.00001);
			}
		}
	}
	printf("made it %lf %lf %lf\n", total1 * dv, total2 * dv, total3);
	mkl_free(x);
	return 0;
}
