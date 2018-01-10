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

void min_cart_path(double* coord, double* center, double* lattice, double* path, double* r) {
	*r = INFINITY;
	double testvec[3];
	double testdist;
	for (int i = -1; i <= 1; i++) {
		for (int j = -1; j <= 1; j++) {
			for (int k = -1; k <= 1; k++) {
				testvec[0] = coord[0] + i - center[0];
				testvec[1] = coord[1] + j - center[1];
				testvec[2] = coord[2] + k - center[2];
				frac_to_cartesian(testvec, lattice);
				testdist = mag(testvec);
				if (testdist < *r) {
					path[0] = testvec[0];
					path[1] = testvec[1];
					path[2] = testvec[2];
					*r = testdist;
				}
			}
		}
	}	
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

void free_ppot(ppot_t* pp) {
	for (int i = 0; i < pp->num_projs; i++) {
		free(pp->funcs[i].proj);
		free(pp->funcs[i].pswave);
		free(pp->funcs[i].aewave);
	}
	free(pp->funcs);
	free(pp->wave_grid);
	free(pp->proj_grid);
	free(pspw_overlap_matrix);
	free(aepw_overlap_matrix);
	free(diff_overlap_matrix);
}

void free_real_proj(real_proj_t* proj) {
	free(proj->values);
}

void free_real_proj_site(real_proj_site_t* site) {
	for (int i = 0; i < site->total_projs; i++) {
		free_real_proj(site->projs + i);
	}
	free(sites->projs);
	free(site->indices);
}

void free_pswf(pswf_t* wf) {
	for (int i = 0; i < wf->nwk * wf->nspin; i++)
		free_kpoint(wf->kpts[i]);
	free(wf->kpts);
	free(wf->G_bounds);
	free(wf->lattice);
	free(wf->reclattice);
	free(wf);
}

void free_ptr(void* ptr) {
	free(ptr);
}

void free_real_proj_site_list(real_proj_site_t* sites, int length) {
	for (int i = 0; i < length; i++) {
		free_real_proj_site(sites + i);
	}
	free(sites);
}

void free_ppot_list(ppot_t* pps, int length) {
	for (int i = 0; i < length; i++) {
		free_ppot(pps + i);
	}
	free(pps);
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
		total += pow(x, 2*n-l-m) * fac(2*n) / fac(2*n-l-m) / fac(n) / fac(l-n) * pow(-1, l-n);
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

	double temp[3] = {0,0,0};
	double r = 0;
	min_cart_path(pos, ion_pos, lattice, temp, &r);
	double radial_val = funcs.proj[(int)(r/rmax*100)];
	if (r == 0) return Ylm(funcs.l, m, 0, 0) * radial_val;
	double theta = 0, phi = 0;
	printf("ERROR %lf %lf\n", mag(temp), r);
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
