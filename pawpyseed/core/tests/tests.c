#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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

	pswf_t* wf = read_wavefunctions(wavecar, kpt_weights);
	double complex* x = (double complex*) mkl_calloc(fftg[0]*fftg[1]*fftg[2],
		sizeof(double complex), 64);
	fft3d(x, wf->G_bounds, wf->lattice, wf->kpts[0]->k, wf->kpts[0]->Gs,
		wf->kpts[0]->bands[0]->Cs, wf->kpts[0]->bands[0]->num_waves, fftg);
	int* Gs = wf->kpts[0]->Gs;
	float complex* Cs = wf->kpts[0]->bands[0]->Cs;
	double inv_sqrt_vol = pow(determinant(wf->lattice), -0.5);
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
				if (cabs(x[ind] - temp) > 1e-5)
					return -1;
			}
		}
	}

	printf("FFTCHECK ASSERTS\n");
	float complex* CAs = (float complex*) calloc(wf->kpts[0]->num_waves, sizeof(float complex));
	fwd_fft3d(x, wf->G_bounds, wf->lattice, wf->kpts[0]->k, wf->kpts[0]->Gs,
		CAs, wf->kpts[0]->bands[0]->num_waves, fftg);
	for (int w = 0; w < wf->kpts[0]->num_waves; w++) {
		if (cabs(CAs[w] - wf->kpts[0]->bands[0]->Cs[w]) > 1e-5)
			return -2;
	}
	free(CAs);

	mkl_free(x);
	return 0;
}

void proj_check(int BAND_NUM, int KPOINT_NUM,
	pswf_t* wf, int* fftg, int* labels, double* coords) {

	ppot_t* pps = wf->pps;
	double complex* x = mkl_calloc(fftg[0]*fftg[1]*fftg[2], sizeof(double complex), 64);
	//printf("START FT\n");
	fft3d(x, wf->G_bounds, wf->lattice, wf->kpts[KPOINT_NUM]->k,
		wf->kpts[KPOINT_NUM]->Gs, wf->kpts[KPOINT_NUM]->bands[BAND_NUM]->Cs,
		wf->kpts[KPOINT_NUM]->bands[BAND_NUM]->num_waves, fftg);
	//printf("FINISH FT\n");
	double* lattice = wf->lattice;
	double vol = determinant(lattice);
	double dv = vol / (fftg[0]*fftg[1]*fftg[2]);
	for (int i = 0; i < fftg[0]; i++) {
		double frac[3] = {0,0,0};
		double kdotr = 0;
		for (int j = 0; j < fftg[1]; j++) {
			for (int k = 0; k < fftg[2]; k++) {
				frac[0] = (double) i / fftg[0];
				frac[1] = (double) j / fftg[1];
				frac[2] = (double) k / fftg[2];
				kdotr = dot(wf->kpts[KPOINT_NUM]->k, frac);
				x[i*fftg[1]*fftg[2] + j*fftg[2] + k] *= cexp(2*PI*I*kdotr);
			}
		}
	}

	double complex* y = (double complex*) malloc(fftg[0]*fftg[1]*fftg[2] * sizeof(double complex));
	memcpy(y, x, fftg[0]*fftg[1]*fftg[2] * sizeof(double complex));

	double err=0, err2=0;
	double normx=0, normy=0;
	int num_sites = wf->num_sites;
	#pragma omp parallel for
	for (int p = 0; p < num_sites; p++) {
		int ind;
		double serr=0, serr2=0;
		double snormx=0, snormy=0;
		projection_t pros = wf->kpts[KPOINT_NUM]->bands[BAND_NUM]->projections[p];
		//printf("READ PROJECTIONS\n");
		ppot_t pp = pps[labels[p]];
		double rmax = pp.wave_grid[pp.wave_gridsize-1];
		double res[3] = {0,0,0};
		vcross(res, lattice+3, lattice+6);
		int grid1 = (int) (mag(res) * rmax / vol * fftg[0]) + 1;
		vcross(res, lattice+0, lattice+6);
		int grid2 = (int) (mag(res) * rmax / vol * fftg[1]) + 1;
		vcross(res, lattice+0, lattice+3);
		int grid3 = (int) (mag(res) * rmax / vol * fftg[2]) + 1;
		int center1 = (int) round(coords[3*p+0] * fftg[0]);
		int center2 = (int) round(coords[3*p+1] * fftg[1]);
		int center3 = (int) round(coords[3*p+2] * fftg[2]);
		//printf("FINISH SETUP %d\n%d %d %d\n%d %d %d\n",p, center1, center2, center3, grid1, grid2, grid3);
		for (int i = -grid1 + center1; i <= grid1 + center1; i++) {
			double frac[3] = {0,0,0};
			double testcoord[3] = {0,0,0};
			int ii=0, jj=0, kk=0;
			double phasecoord[3] = {0,0,0};
			double phase = 0;
			for (int j = -grid2 + center2; j <= grid2 + center2; j++) {
				for (int k = -grid3 + center3; k <= grid3 + center3; k++) {
					testcoord[0] = (double) i / fftg[0] - coords[3*p+0];
					testcoord[1] = (double) j / fftg[1] - coords[3*p+1];
					testcoord[2] = (double) k / fftg[2] - coords[3*p+2];
					frac_to_cartesian(testcoord, lattice);
					if (mag(testcoord) < rmax) {
						
						ii = (i%fftg[0] + fftg[0]) % fftg[0];
						jj = (j%fftg[1] + fftg[1]) % fftg[1];
						kk = (k%fftg[2] + fftg[2]) % fftg[2];
						frac[0] = (double) ii / fftg[0];
						frac[1] = (double) jj / fftg[1];
						frac[2] = (double) kk / fftg[2];
						phasecoord[0] = coords[3*p+0] + ((ii-i) / fftg[0]);
						phasecoord[1] = coords[3*p+1] + ((jj-j) / fftg[1]);
						phasecoord[2] = coords[3*p+2] + ((kk-k) / fftg[2]);
						phase = dot(phasecoord, wf->kpts[KPOINT_NUM]->k);
						ind = ii*fftg[1]*fftg[2] + jj*fftg[2] + kk;
						x[ii*fftg[1]*fftg[2] + jj*fftg[2] + kk] = 0;
						for (int n = 0; n < pros.total_projs; n++) {
							x[ii*fftg[1]*fftg[2] + jj*fftg[2] + kk] +=
								wave_value2(pp.wave_grid,
								pp.funcs[pros.ns[n]].pswave,
								pp.funcs[pros.ns[n]].pswave_spline,
								pp.wave_gridsize,
								pros.ls[n], pros.ms[n],
								testcoord)
								* pros.overlaps[n] * cexp(2*PI*I*phase);
						}
						serr += pow(cabs(x[ind] - y[ind]), 2);
						serr2 += pow(cabs(x[ind]) - cabs(y[ind]), 2);
						snormx += pow(cabs(x[ind]), 2);
						snormy += pow(cabs(y[ind]), 2);
					}
				}
			}
		}
		#pragma omp critical
		{
			err += serr;
			err2 += serr2;
			normx += snormx;
			normy += snormy;
		}
	}

	printf("err magerr, normx normy %lf %lf %lf %lf\n", err/normy, err2/normy, normx, normy);

	mkl_free(x);
	free(y);
}
