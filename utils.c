#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <omp.h>
#include <time.h>
#include "utils.h"

#define PI 3.14159265359

void vcross(double* res, double* top, double* bottom) {
	res[0] = top[1] * bottom[2] - top[2] * bottom[1];
	res[1] = top[2] * bottom[0] - top[0] * bottom[2];
	res[2] = top[0] * bottom[1] - top[1] * bottom[0];
}

double dot(double* x1, double* x2) {
	return x1[0] * x2[0] + x1[1] * x2[1] + x1[2] * x2[2];
}

double mag(double* x1) {
	return pow(dot(x1, x1), 0.5);
}

double determinant(double* m) {
	return m[0] * m[4] * m[8]
		+  m[1] * m[5] * m[6]
		+  m[2] * m[3] * m[7]
		-  m[2] * m[4] * m[6]
		-  m[1] * m[3] * m[8]
		-  m[0] * m[5] * m[7];
}

double dist_from_frac(double* coords1, double* coords2, double* lattice) {
	double f1 = fmin(fabs(coords1[0]-coords2[0]), 1-fabs(coords1[0]-coords2[0]));
	double f2 = fmin(fabs(coords1[1]-coords2[1]), 1-fabs(coords1[1]-coords2[1]));
	double f3 = fmin(fabs(coords1[2]-coords2[2]), 1-fabs(coords1[2]-coords2[2]));
	return pow(pow(f1*lattice[0]+f2*lattice[3]+f3*lattice[6], 2)
		+ pow(f1*lattice[1]+f2*lattice[4]+f3*lattice[7], 2)
		+ pow(f1*lattice[2]+f2*lattice[5]+f3*lattice[8], 2), 0.5);
}

void free_kpoint(kpoint_t* kpt) {
	for (int b = 0; b < kpt->num_bands; b++) {
		band_t* curr_band = kpt->bands[b];
		free(curr_band->Cs);
		//free(curr_band->Gs);
		//free(curr_band->C_grid);
		free(curr_band);
	}
	//printf("ya");
	free(kpt->Gs);
	free(kpt->bands);
	//free(kpt->k);
	free(kpt);
}

void free_pswf(pswf_t* wf) {
	for (int i = 0; i < wf->nwk * wf->nspin; i++)
		free_kpoint(wf->kpts[i]);
}

double* get_occs(pswf_t* wf) {
	kpoint_t** kpts = wf->kpts;
	double* occs = (double*) malloc(wf->nwk*wf->nband*wf->nspin*sizeof(double));
	int NUM_KPTS = wf->nwk * wf->nspin;
	for (int kpt_num = 0; kpt_num < NUM_KPTS; kpt_num++) {
		for (int band_num = 0; band_num < wf->nband; band_num++) {
			occs[band_num*NUM_KPTS+kpt_num] = kpts[kpt_num]->bands[band_num]->occ;
		}
	}
}

int get_nband(pswf_t* wf) {
	return wf->nband;
}

int get_nwk(pswf_t* wf) {
	return wf->nwk;
}

int get_nspin(pswf_t* wf) {
	return wf->nspin;
}

double legendre(int l, double x) {
	double p0 = 1;
	double p1 = x;
	if (l == 0) return p0;
	if (l == 1) return p1;
	double pf;
	for (int currl = 1; currl < l; currl++) {
		pf = ((2*currl+1) * x * p1 - currl * p0) / (l+1);
		p1 = pf;
		p0 = p1;
	}
	return pf;
}

double fac(int n) {
	int m = 1;
	int t = 1;
	while (m <= n) {
		t *= m;
		m++;
	}
	return (double)t;
}

double complex Ylm(int l, int m, double theta, double phi) {
	return pow(-1, m) * pow((2*l+1)/(4*PI)*fac(l-m)/fac(l+m), 0.5) * legendre(l, cos(theta)) * cexp(I*m*phi);
}

void ALLOCATION_FAILED() {
	printf("ALLOCATION FAILED");
	exit(-1);
}
