#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <omp.h>
#include <time.h>
#include "utils.h"

#define PI 3.14159265358979323846

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
	double testvec[3]= {0,0,0};
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

void frac_from_spherical(double* ion_frac, double r, double theta, double phi,
	double* lattice, double* reclattice, double* result) {

	double cart[3];
	cart[0] = r * sin(theta) * cos(phi);
	cart[1] = r * sin(theta) * sin(phi);
	cart[2] = r * cos(theta);
	cartesian_to_frac(cart, reclattice);
	result[0] = fmod(cart[0] + ion_frac[0], 1.00);
	result[1] = fmod(cart[1] + ion_frac[1], 1.00);
	result[2] = fmod(cart[2] + ion_frac[2], 1.00);
	if (result[0] < 0) result[0] += 1;
	if (result[1] < 0) result[1] += 1;
	if (result[2] < 0) result[2] += 1;
}

double complex trilinear_interpolate(double complex* c, double* frac, int* fftg) {
	//values: c000, c001, c010, c011, c100, c101, c110, c111
	double d[3];
	d[0] = fmod(frac[0] * fftg[0], 1.0);
	d[1] = fmod(frac[1] * fftg[1], 1.0);
	d[2] = fmod(frac[2] * fftg[2], 1.0);
	//printf("%lf %lf %lf\n", d[0], d[1], d[2]);
	double complex c00 = c[0] * (1-d[0]) + c[4] * d[0];
	double complex c01 = c[1] * (1-d[0]) + c[5] * d[0];
	double complex c10 = c[2] * (1-d[0]) + c[6] * d[0];
	double complex c11 = c[3] * (1-d[0]) + c[7] * d[0];

	double complex c0 = c00 * (1-d[1]) + c10 * d[1];
	double complex c1 = c01 * (1-d[1]) + c11 * d[1];

	return c0 * (1-d[2]) + c1 * d[2];
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

void cartesian_to_frac(double* coord, double* reclattice) {
	double temp[3] = {0,0,0};
	temp[0] = coord[0] * reclattice[0] + coord[1] * reclattice[1] + coord[2] * reclattice[2];
	temp[1] = coord[0] * reclattice[3] + coord[1] * reclattice[4] + coord[2] * reclattice[5];
	temp[2] = coord[0] * reclattice[6] + coord[1] * reclattice[7] + coord[2] * reclattice[8];
	coord[0] = temp[0] / 2 / PI;
	coord[1] = temp[1] / 2 / PI;
	coord[2] = temp[2] / 2 / PI;
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
	free(kpt->k);
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
	if (pp->pspw_overlap_matrix != NULL) free(pp->pspw_overlap_matrix);
	if (pp->aepw_overlap_matrix != NULL) free(pp->aepw_overlap_matrix);
	if (pp->diff_overlap_matrix != NULL) free(pp->diff_overlap_matrix);
}

void free_real_proj(real_proj_t* proj) {
	free(proj->values);
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

void free_real_proj_site(real_proj_site_t* site) {
	for (int i = 0; i < site->total_projs; i++) {
		free_real_proj(site->projs + i);
	}
	free(site->projs);
	free(site->indices);
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

int min(int a, int b) {
	if (a > b)
		return b;
	else
		return a;
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
	return total * pow(-1, m) * pow(1 - x * x, m/2.0) / pow(2, l);
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
	double complex multiplier = 0;
	if (m == 0) multiplier = 1;
	else if (m < 0) multiplier = pow(2.0, 0.5) * cos(-m*phi);
	else multiplier = pow(-1,m) * pow(2.0, 0.5) * sin(m*phi);
	return pow((2*l+1)/(4*PI)*fac(l-m)/fac(l+m), 0.5) *
		legendre(l, m, cos(theta))*multiplier;// cexp(I*m*phi);
}

double complex Ylm2(int l, int m, double costheta, double phi) {
	//printf("%lf %lf %lf\n", pow((2*l+1)/(4*PI)*fac(l-m)/fac(l+m), 0.5), legendre(l, m, cos(theta)),
	//	creal(cexp(I*m*phi)));
	double complex multiplier = 0;
        if (m == 0) multiplier = 1;
        else if (m < 0) multiplier = pow(2.0, 0.5) * cos(-m*phi);
        else multiplier = pow(-1,m) * pow(2.0, 0.5) * sin(m*phi);
	return pow((2*l+1)/(4*PI)*fac(l-m)/fac(l+m), 0.5) *
		legendre(l, m, costheta) *multiplier;//* cexp(I*m*phi);
}

double proj_interpolate(double r, double rmax, double* x, double* proj, double** proj_spline) {
	int ind = min((int)(r/rmax*100), 98);
	double rem = r - x[ind];
	double radval = proj[ind] + rem * (proj_spline[0][ind] +
						rem * (proj_spline[1][ind] +
						rem * proj_spline[2][ind]));
	return radval;
}

double wave_interpolate(double r, int size, double* x, double* f, double** wave_spline) {
	int ind = min((int) (log(r/x[0]) / log(x[1]/x[0])), size-2);
	double rem = r - x[ind];
	return f[ind] + rem * (wave_spline[0][ind] + 
				rem * (wave_spline[1][ind] +
				rem * wave_spline[2][ind]));
}

double complex proj_value(funcset_t funcs, double* x, int m, double rmax,
	double* ion_pos, double* pos, double* lattice) {

	double temp[3] = {0,0,0};
	double r = 0;
	min_cart_path(pos, ion_pos, lattice, temp, &r);

	double radial_val = proj_interpolate(r, rmax, x, funcs.proj, funcs.proj_spline);
	if (r == 0) return Ylm(funcs.l, m, 0, 0) * radial_val;
	double theta = 0, phi = 0;
	theta = acos(temp[2]/r);
	if (r - fabs(temp[2]) == 0) phi = 0;
	else phi = acos(temp[0] / pow(temp[0]*temp[0] + temp[1]*temp[1], 0.5));
	if (temp[1] < 0) phi = 2*PI - phi;
	//printf("inp %d %lf %d\n", m, rmax, funcs.l);
	double complex sph_val = Ylm(funcs.l, m, theta, phi);
	//printf("out %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", phi, theta, r, temp[0], temp[1], temp[2], radial_val, creal(sph_val), cimag(sph_val));
	return radial_val * sph_val;
}

double complex onto_partial_wave(double* pvec, double* reclattice, double complex* kwave,
	double kappamin, double dkappa, int N, double** wave_spline) {

	double temp[3] = {0,0,0};
	temp[0] = reclattice[0] * pvec[0] + reclattice[3] * pvec[1] + reclattice[6] * pvec[2];
	temp[0] = reclattice[1] * pvec[0] + reclattice[4] * pvec[1] + reclattice[7] * pvec[2];
	temp[0] = reclattice[2] * pvec[0] + reclattice[5] * pvec[1] + reclattice[8] * pvec[2];

	int n = (int) ((log(mag(temp)) - kappamin) / dkappa);
	double rem = mag(temp) - exp(kappamin + n * dkappa);
	double val = kwave[n] + rem * (wave_spline[0][n] +
				rem * (wave_spline[1][n] + 
				rem * wave_spline[2][n]));
	return val;

}

//adapted from VASP source code
double** spline_coeff(double* x, double* y, int N) {
	double** coeff = (double**) malloc(3 * sizeof(double*));
	coeff[0] = (double*) malloc(N * sizeof(double));
	coeff[1] = (double*) malloc(N * sizeof(double));
	coeff[2] = (double*) malloc(N * sizeof(double));

	printf("pl %d\n", N);
	double d1p1 = (y[1] - y[0]) / (x[1] - x[0]);
	if (d1p1 > 0.99E30) {
		coeff[1][0] = 0;
		coeff[0][0] = 0;
	}
	else {
		coeff[1][0] = -0.5;
		coeff[0][0] = (3 / (x[1] - x[0])) * ((y[1] - y[0]) / (x[1] - x[0]) - d1p1);
	}

	double s = 0, r = 0;

	for (int i = 1; i < N - 1; i++) {
		s = (x[i] - x[i-1]) / (x[i+1] - x[i-1]);
		r = s * coeff[1][i-1] + 2;
		coeff[1][i] = (s - 1) / r;
		coeff[0][i] = (6 * ( (y[i+1] - y[i]) / (x[i+1] - x[i]) -
			(y[i] - y[i-1]) / (x[i] - x[i-1])) /
			(x[i+1] - x[i-1]) - s*coeff[0][i-1]) / r;
	}

	coeff[0][N-1] = 0;
	coeff[1][N-1] = 0;
	coeff[2][N-1] = 0;

	for (int i = N-2; i >= 0; i--) {
		coeff[1][i] = coeff[1][i] * coeff[1][i+1] + coeff[0][i];
	}

	for (int i = 0; i < N-1; i++) {
		s = x[i+1] - x[i];
		r = (coeff[1][i+1] - coeff[1][i]) / 6;
		coeff[2][i] = r / s;
		coeff[1][i] = coeff[1][i] / 2;
		coeff[0][i] = (y[i+1]-y[i]) / s - (coeff[1][i] + r) * s;
	}

	return coeff;
}

void frac_from_index(int index, double* coord, int* fftg) {
	int t1 = index / (fftg[1] * fftg[2]);
	int t2 = index % (fftg[1] * fftg[2]);
	int t3 = t2 % fftg[2];
	t2 /= fftg[2];
	coord[0] = ((double) t1) / fftg[0];
	coord[1] = ((double) t2) / fftg[1];
	coord[2] = ((double) t3) / fftg[2];
}

void direction(double* cart, double* dir) {
	double theta = 0, phi = 0;
	double r = mag(cart);
	theta = acos(cart[2]/r);
	if (r - fabs(cart[2]) == 0) phi = 0;
	else phi = acos(cart[0] / pow(cart[0]*cart[0] + cart[1]*cart[1], 0.5));
	if (cart[1] < 0) phi = 2*PI - phi;
	dir[0] = theta;
	dir[1] = phi;
}

double sph_bessel(double k, double r, int l) {
	double x = k * r;
	if (l == 0)
		return sin(x) / x;
	else if (l == 1)
		return sin(x) / (x*x) - cos(x) / x;
	else if (l == 2)
		return (3 / (x*x) -1) * sin(x) / x - 3 * cos(x) / (x*x);
	else if (l == 3)
		return (15 / (x*x*x) - 6 / x) * sin(x) / x - (15 / (x*x) -1) * cos(x) / x;
	else
		printf("ERROR: sph_bessel l too high");
		return 0;
}

double complex rayexp(double* kpt, int* Gs, float complex* Cs, int l, int m,
	int num_waves, double* grid, double* wave, double* spline, double* ionp) {

	double complex result = 0;
	double pvec[3] = {0,0,0};
	double phase = 0;
	for (int w = 0; w < num_waves; w++) {
		pvec[0] = kpt[0] + Gs[3*w+0];
		pvec[1] = kpt[1] + Gs[3*w+1];
		pvec[2] = kpt[2] + Gs[3*w+2];
		phase = cexp(2*PI*I*dot(ionp, pvec));
		result += phase * Cs[w] * sum_terms[(2*l+1)*w+l+m];
	}

	return result * 4 * PI * cpow(I, l);
}

double complex* rayexp_terms(double* kpt, int* Gs, int num_waves,
	int l, int wave_gridsize, double* grid,
	double* wave, double* spline, double* reclattice) {

	double complex ylmdir = 0;
	double k = 0;
	double complex* terms = (double complex*) malloc((2*l+1) * num_waves * sizeof(double complex));

	double pvec[3] = {0,0,0};
	double phat[2] = {0,0};
	for (int w = 0; w < num_waves; w++) {
		pvec[0] = kpt[0] + Gs[3*w+0];
		pvec[1] = kpt[1] + Gs[3*w+1];
		pvec[2] = kpt[2] + Gs[3*w+2];
		frac_to_cartesian(pvec, reclattice);
		k = mag(pvec);
		direction(pvec, phat);
		ylmdir = conj(Ylm(l, m, phat[0], phat[1]));
		overlap = wave_interpolate(k, wave_gridsize, grid, wave, spline);
		terms[(2*l+1)*w+l+m] = ylmdir * overlap;
	}
	return terms;
}

void generate_rayleigh_expansion_terms(pswf_t* wf, ppot_t* pps, int num_elems) {
	for (int k_num = 0; k_num < wf->nwk * wf->nspin; k_num++) {
		kpoint_t* kpt = wf->kpts[k_num];
		kpt->expansion = (rayleigh_set_t**) malloc(num_elems * sizeof(rayleigh_set_t*));
		for (int i = 0; i < num_elems; i++) {
			ppot_t pp = pps[i];
			kpt->expansion[i] = (rayleigh_set_t*) malloc(pp.num_projs * sizeof(rayleigh_set_t));
			for (int j = 0; j < pp.num_projs; j++) {
				double complex* terms = rayexp_terms(kpt->k, kpt->Gs, kpt->num_waves,
					pp.funcs[j].l, pp.wave_gridsize, pp.wave_grid,
					pp.funcs[j].aewave, pp.funcs[j].pswave, wf->reclattice);
				kpt->expansion[i][j].terms = terms;
				kpt->expansion[i][j].l = pp.funcs[j].l;
			}
		}
	}
}

void ALLOCATION_FAILED() {
	printf("ALLOCATION FAILED\n");
	exit(-1);
}
