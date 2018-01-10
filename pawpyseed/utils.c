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

void frac_to_cartesian(double* coord, double* lattice) {
	double temp[3] = {0,0,0};
	temp[0] = coord[0] * lattice[0] + coord[1] * lattice[3] + coord[2] * lattice[6];
	temp[1] = coord[0] * lattice[1] + coord[1] * lattice[4] + coord[2] * lattice[7];
	temp[2] = coord[0] * lattice[2] + coord[1] * lattice[5] + coord[2] * lattice[8];
	coord[0] = temp[0];
	coord[1] = temp[1];
	coord[2] = temp[2];
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
	return occs;
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

double legendre(int l, int m, double x) {
	double total = 0;
	if (m < 0) return pow(-1.0, m) * fac(l+m) / fac(l-m) * legendre(l, -m, x);
	for (int n = l; n >= 0 && 2*n-l-m >= 0; n--) {
		total += pow(x, 2*n-l-m) * fac(2*n) / fac(2*n-l-m) / fac(n) / fac(l-n) * pow(-1, l+m-n);
	}
	return total * pow(1 - x * x, m/2.0) / pow(2, l);
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
	//printf("%lf %lf %lf\n", pow((2*l+1)/(4*PI)*fac(l-m)/fac(l+m), 0.5), legendre(l, m, cos(theta)),
	//	creal(cexp(I*m*phi)));
	return pow(-1, m) * pow((2*l+1)/(4*PI)*fac(l-m)/fac(l+m), 0.5) *
		legendre(l, m, cos(theta)) * cexp(I*m*phi);
}

double complex proj_value(funcset_t funcs, int m, double rmax,
	double* ion_pos, double* pos, double* lattice) {

	double r = dist_from_frac(ion_pos, pos, lattice);
	double radial_val = funcs.proj[(int)(r/rmax*100)];
	if (r == 0) return Ylm(funcs.l, m, 0, 0) * radial_val;
	double theta = 0, phi = 0;
	double temp[3] = {0,0,0};
	for (int i = 0; i < 3; i++) {
		temp[i] = pos[i] - ion_pos[i];
		if (fabs(pos[i] - ion_pos[i]) > fabs(pos[i] + 1 - ion_pos[i]))
			temp[i] += 1;
		if (fabs(temp[i]) > fabs(pos[i] - 1 - ion_pos[i]))
			temp[i] = pos[i] - 1 - ion_pos[i];
	}
	frac_to_cartesian(temp, lattice);
	theta = acos(temp[2]/r);
	if (r - fabs(temp[2]) == 0) phi = 0;
	else phi = acos(temp[0] / pow(temp[0]*temp[0] + temp[1]*temp[1], 0.5));
	if (temp[1] < 0) phi = 2*PI - phi;
	//printf("inp %d %lf %d\n", m, rmax, funcs.l);
	double complex sph_val = Ylm(funcs.l, m, theta, phi);
	//printf("out %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", phi, theta, r, temp[0], temp[1], temp[2], radial_val, creal(sph_val), cimag(sph_val));
	return radial_val * sph_val;
}



void ALLOCATION_FAILED() {
	printf("ALLOCATION FAILED\n");
	exit(-1);
}
