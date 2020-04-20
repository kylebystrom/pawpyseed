#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <math.h>
#include <omp.h>
#include <time.h>
#include <mkl.h>
#include <mkl_types.h>
#include "utils.h"

#define PI 3.14159265358979323846
#define CCONST 0.262465831
#define OPSIZE 9

void affine_transform(double* out, double* op, double* inv) {
	//0 1 2 3
	//4 5 6 7
	//8 9 10 11
	//12 13 14 15
	double in[4] = {inv[0], inv[1], inv[2], 1};
	out[0] = op[0]*in[0] + op[1]*in[1] + op[2]*in[2] + op[3]*in[3];
	out[1] = op[4]*in[0] + op[5]*in[1] + op[6]*in[2] + op[7]*in[3];
	out[2] = op[8]*in[0] + op[9]*in[1] + op[10]*in[2] + op[11]*in[3];
}

void rotation_transform(double* out, double* op, double* inv) {
	double in[3] = {inv[0], inv[1], inv[2]};
	out[0] = op[0]*in[0] + op[1]*in[1] + op[2]*in[2];
	out[1] = op[3]*in[0] + op[4]*in[1] + op[5]*in[2];
	out[2] = op[6]*in[0] + op[7]*in[1] + op[8]*in[2];
}

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

void trilinear_interpolate_values(double complex* x, double* frac, int* fftg, double complex* values) {
	//values: c000, c001, c010, c011, c100, c101, c110, c111
	int i = (int) (frac[0] * fftg[0]);
	int j = (int) (frac[1] * fftg[1]);
	int k = (int) (frac[2] * fftg[2]);
	int ip = (i+1)%fftg[0];
	int jp = (j+1)%fftg[1];
	int kp = (k+1)%fftg[2];

	int ind[8];
	ind[0] = i*fftg[1]*fftg[2] + j*fftg[2] + k;
	ind[1] = i*fftg[1]*fftg[2] + j*fftg[2] + kp;
	ind[2] = i*fftg[1]*fftg[2] + jp*fftg[2] + k;
	ind[3] = i*fftg[1]*fftg[2] + jp*fftg[2] + kp;
	ind[4] = ip*fftg[1]*fftg[2] + j*fftg[2] + k;
	ind[5] = ip*fftg[1]*fftg[2] + j*fftg[2] + kp;
	ind[6] = ip*fftg[1]*fftg[2] + jp*fftg[2] + k;
	ind[7] = ip*fftg[1]*fftg[2] + jp*fftg[2] + kp;
	for (int n = 0; n < 8; n++) values[n] = x[ind[n]];
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

void free_rayleigh_set_list(rayleigh_set_t* sets, int num_projs) {
	for (int i = 0; i < num_projs; i++)
		free(sets[i].terms);
	free(sets);
}

void free_projection_list(projection_t* projlist, int num) {
	for (int i = 0; i < num; i++) {
		free(projlist[i].ms);
		free(projlist[i].ls);
		free(projlist[i].ns);
		free(projlist[i].overlaps);
	}
	free(projlist);
}

void clean_wave_projections(pswf_t* wf) {

	for (int i = 0; i < wf->nwk * wf->nspin; i++) {
		kpoint_t* kpt = wf->kpts[i];
		for (int b = 0; b < kpt->num_bands; b++) {
			if (kpt->bands[b]->wave_projections != NULL) {
				free_projection_list(kpt->bands[b]->wave_projections, wf->wp_num);
				kpt->bands[b]->wave_projections = NULL;
			}
		}
	}

}

void free_kpoint(kpoint_t* kpt, int num_elems, int num_sites, int wp_num, int* num_projs) {
	for (int b = 0; b < kpt->num_bands; b++) {
		band_t* curr_band = kpt->bands[b];
		free(curr_band->Cs);
		if (curr_band->projections != NULL) {
			free_projection_list(curr_band->projections, num_sites);
		}
		if (curr_band->wave_projections != NULL) {
			free_projection_list(curr_band->wave_projections, wp_num);
		}
		if (curr_band->up_projections != NULL) {
			free_projection_list(curr_band->up_projections, num_sites);
		}
		if (curr_band->down_projections != NULL) {
			free_projection_list(curr_band->down_projections, num_sites);
		}
		if (curr_band->CRs != NULL) {
			mkl_free(curr_band->CRs);
		}
		if (curr_band->CAs != NULL) {
			mkl_free(curr_band->CAs);
		}
		free(curr_band);
	}
	if (kpt->expansion != NULL) {
		for (int i = 0; i < num_elems; i++)
			free_rayleigh_set_list(kpt->expansion[i], num_projs[i]);
		free(kpt->expansion);
	}
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
		free(pp->funcs[i].diffwave);
		free(pp->funcs[i].kwave);
		for (int j = 0; j < 3; j++) {
			free(pp->funcs[i].proj_spline[j]);
			free(pp->funcs[i].aewave_spline[j]);
			free(pp->funcs[i].pswave_spline[j]);
			free(pp->funcs[i].diffwave_spline[j]);
			free(pp->funcs[i].kwave_spline[j]);
		}
		free(pp->funcs[i].proj_spline);
		free(pp->funcs[i].aewave_spline);
		free(pp->funcs[i].pswave_spline);
		free(pp->funcs[i].diffwave_spline);
		free(pp->funcs[i].kwave_spline);
	}
	free(pp->funcs);
	free(pp->wave_grid);
	free(pp->proj_grid);
	free(pp->pspw_overlap_matrix);
	free(pp->aepw_overlap_matrix);
	free(pp->diff_overlap_matrix);
}

void free_real_proj(real_proj_t* proj) {
	free(proj->values);
}

void free_pswf(pswf_t* wf) {
	for (int i = 0; i < wf->nwk * wf->nspin; i++)
		free_kpoint(wf->kpts[i], wf->num_elems, wf->num_sites, wf->wp_num, wf->num_projs);
	if (wf->overlaps != NULL) {
		for (int i = 0; i < wf->num_aug_overlap_sites; i++)
			free(wf->overlaps[i]);
		free(wf->overlaps);
	}
	if (wf->num_projs != NULL) {
		free(wf->num_projs);
	}
	free(wf->kpts);
	free(wf->G_bounds);
	free(wf->lattice);
	free(wf->reclattice);
	if (wf->pps != NULL) {
		free_ppot_list(wf->pps, wf->num_elems);
	}
	free(wf);
}

void free_real_proj_site(real_proj_site_t* site) {
	for (int i = 0; i < site->total_projs; i++) {
		free_real_proj(site->projs + i);
	}
	free(site->projs);
	free(site->indices);
	free(site->coord);
	free(site->paths);
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

int max(int a, int b) {
	if (a < b)
		return b;
	else
		return a;
}

double* get_occs(pswf_t* wf) {
	kpoint_t** kpts = wf->kpts;
	double* occs = (double*) malloc(wf->nwk*wf->nband*wf->nspin*sizeof(double));
	CHECK_ALLOCATION(occs);
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

double get_encut(pswf_t* wf) {
	return wf->encut;
}

int is_ncl(pswf_t* wf) {
	return wf->is_ncl;
}

double get_energy(pswf_t* wf, int band, int kpt, int spin) {
	return wf->kpts[kpt+spin*wf->nwk]->bands[band]->energy;
}

double get_occ(pswf_t* wf, int band, int kpt, int spin) {
	return wf->kpts[kpt+spin*wf->nwk]->bands[band]->occ;
}

void set_num_sites(pswf_t* wf, int nsites) {
	wf->num_sites = nsites;
}

double legendre(int l, int m, double x) {
	double total = 0;
	if (m < 0) return pow(-1.0, m) * fac(l+m) / fac(l-m) * legendre(l, -m, x);
	for (int n = l; n >= 0 && 2*n-l-m >= 0; n--) {
		total += pow(x, 2*n-l-m) * fac(2*n) / fac(2*n-l-m) / fac(n) / fac(l-n) * pow(-1, l-n);
	}
	return total * pow(-1, m) * pow(1 - x * x, m/2.0) / pow(2, l);
}

void legendre_coeff(double* ptr, int l, int m) {
	// assumes ptr is cleared
	double prefac = pow(-1, m) / pow(2, l);
	if (m < 0) {
		prefac *= pow(-1.0, m) * fac(l+m) / fac(l-m);
	}
	m = abs(m);
	for (int n = l; n >= 0 && 2*n-l-m >= 0; n--) {
		ptr[n] = fac(2*n) / fac(2*n-l-m) / fac(n) / fac(l-n) * pow(-1, l-n) * prefac;
	}
}

double* legendre_product(int l1, int l2, int m1, int m2) {
	int m = m2 - m1;
	int maxl = l1 + l2;
	double* lp1 = (double*) calloc((l1+1), sizeof(double));
	double* lp2 = (double*) calloc((l2+1), sizeof(double));
	double* polynomial = (double*) calloc((maxl+1), sizeof(double));
	double* test = (double*) calloc((maxl+1), sizeof(double));
	double* coeff = (double*) calloc((maxl+1), sizeof(double));
	legendre_coeff(lp1, l1, m1);
	legendre_coeff(lp2, l2, m2);
	for (int n1 = 0; n1 <= l1; n1++) {
		for (int n2 = 0; n2 <= l2; n2++) {
			polynomial[n1+n2] += lp1[n1] * lp2[n2];
		}
	}
	for (int l = maxl; l >= abs(m); l--) {
		legendre_coeff(test, l, m);
		coeff[l] = polynomial[l] / test[l];
		for (int lp = abs(m); lp <= maxl; lp++) {
			polynomial[lp] -= coeff[l] * test[lp];
			test[lp] = 0;
		}
	}
	free(lp1);
	free(lp2);
	free(polynomial);
	free(test);
	return coeff;
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
	//double complex multiplier = 0;
	//if (m == 0) multiplier = 1;
	//else if (m < 0) multiplier = pow(2.0, 0.5) * cos(-m*phi);
	//else multiplier = pow(-1,m) * pow(2.0, 0.5) * sin(m*phi);
	return pow((2*l+1)/(4*PI)*fac(l-m)/fac(l+m), 0.5) *
		legendre(l, m, cos(theta))* cexp(I*m*phi);
}

double complex Ylm2(int l, int m, double costheta, double phi) {
	//printf("%lf %lf %lf\n", pow((2*l+1)/(4*PI)*fac(l-m)/fac(l+m), 0.5), legendre(l, m, cos(theta)),
	//	creal(cexp(I*m*phi)));
	//double complex multiplier = 0;
	//if (m == 0) multiplier = 1;
	//else if (m < 0) multiplier = pow(2.0, 0.5) * cos(-m*phi);
	//else multiplier = pow(-1,m) * pow(2.0, 0.5) * sin(m*phi);
	return pow((2*l+1)/(4*PI)*fac(l-m)/fac(l+m), 0.5) *
		legendre(l, m, costheta) * cexp(I*m*phi);
}

double proj_interpolate(double r, double rmax, int size, double* x,
	double* proj, double** proj_spline) {
	
	if (r > x[size-1]) return 0;
	if (r < x[0]) return proj[0];
	int ind = min((int)(r/rmax*size), size-2);
	double rem = r - x[ind];
	double radval = proj[ind] + rem * (proj_spline[0][ind] +
						rem * (proj_spline[1][ind] +
						rem * proj_spline[2][ind]));
	return radval;
}

double wave_interpolate(double r, int size, double* x, double* f,
	double** wave_spline) {

	if (r > x[size-1]) return 0;
	if (r < x[0]) return f[0];
	int ind = min((int) (log(r/x[0]) / log(x[1]/x[0])), size-2);
	double rem = r - x[ind];
	return f[ind] + rem * (wave_spline[0][ind] + 
				rem * (wave_spline[1][ind] +
				rem * wave_spline[2][ind]));
}

double complex wave_value(funcset_t funcs, int size, double* x, int m,
	double* ion_pos, double* pos, double* lattice) {

	double temp[3] = {0,0,0};
	double r = 0;
	min_cart_path(pos, ion_pos, lattice, temp, &r);

	double ae_radial_val = wave_interpolate(r, size, x, funcs.aewave, funcs.aewave_spline);
	double ps_radial_val = wave_interpolate(r, size, x, funcs.pswave, funcs.pswave_spline);
	double radial_val = (ae_radial_val - ps_radial_val);
	if (r < x[0])
		radial_val /= x[0];
	else
		radial_val /= r;

	if (r==0) return Ylm(funcs.l, m, 0, 0) * radial_val;
	double theta = 0, phi = 0;
	theta = acos(temp[2]/r);
	if (r - fabs(temp[2]) == 0) phi = 0;
	else phi = acos(temp[0] / pow(temp[0]*temp[0] + temp[1]*temp[1], 0.5));
	if (temp[1] < 0) phi = 2*PI - phi;
	double complex sph_val = Ylm(funcs.l, m, theta, phi);
	return radial_val * sph_val;
}

double complex wave_value2(double* x, double* wave, double** spline, int size,
	int l, int m, double* pos) {

	double r = mag(pos);

	double radial_val = wave_interpolate(r, size, x, wave, spline);
	if (r < x[0])
		radial_val /= x[0];
	else
		radial_val /= r;

	if (r==0) return Ylm(l, m, 0, 0) * radial_val;
	double theta = 0, phi = 0;
	theta = acos(pos[2]/r);
	if (r - fabs(pos[2]) == 0) phi = 0;
	else phi = acos(pos[0] / pow(pos[0]*pos[0] + pos[1]*pos[1], 0.5));
	if (pos[1] < 0) phi = 2*PI - phi;
	double complex sph_val = Ylm(l, m, theta, phi);
	return radial_val * sph_val;
}

double complex proj_value_helper(double r, double rmax, int size,
	double* temp, double* x, double* f, double** s, int l, int m) {

	double radial_val = proj_interpolate(r, rmax, size, x, f, s);
	if (r == 0) return Ylm(l, m, 0, 0) * radial_val;
	double theta = 0, phi = 0;
	theta = acos(temp[2]/r);
	if (r - fabs(temp[2]) == 0) phi = 0;
	else phi = acos(temp[0] / pow(temp[0]*temp[0] + temp[1]*temp[1], 0.5));
	if (temp[1] < 0) phi = 2*PI - phi;
	double complex sph_val = Ylm(l, m, theta, phi);
	return radial_val * sph_val;
}

double complex proj_value(funcset_t funcs, double* x, int m, double rmax,
	int size, double* ion_pos, double* pos, double* lattice) {

	double temp[3] = {0,0,0};
	double r = 0;
	min_cart_path(pos, ion_pos, lattice, temp, &r);

	return proj_value_helper(r, rmax, size, temp, x, funcs.proj,
		funcs.proj_spline, funcs.l, m);
}

double complex smooth_wave_value(funcset_t funcs, double* x, int m, double rmax,
	int size, double* ion_pos, double* pos, double* lattice) {

	double temp[3] = {0,0,0};
	double r = 0;
	min_cart_path(pos, ion_pos, lattice, temp, &r);

	return proj_value_helper(r, rmax, size, temp, x, funcs.smooth_diffwave,
		funcs.smooth_diffwave_spline, funcs.l, m);
}

void setup_site(real_proj_site_t* sites, ppot_t* pps, int num_sites, int* site_nums,
	int* labels, double* coords, double* lattice, int* fftg, int pr0_pw1) {
	
	double vol = determinant(lattice);

	for (int s = 0; s < num_sites; s++) {
		int i = site_nums[s];
		sites[s].index = i;
		sites[s].elem = labels[i];
		sites[s].gridsize = pps[labels[i]].proj_gridsize;
		sites[s].num_projs = pps[labels[i]].num_projs;
		if (pr0_pw1) sites[s].rmax = pps[labels[i]].wave_rmax;
		else sites[s].rmax = pps[labels[i]].rmax;
		sites[s].total_projs = pps[labels[i]].total_projs;
		sites[s].num_indices = 0;
		sites[s].coord = malloc(3 * sizeof(double));
		CHECK_ALLOCATION(sites[s].coord);
		sites[s].coord[0] = coords[3*i+0];
		sites[s].coord[1] = coords[3*i+1];
		sites[s].coord[2] = coords[3*i+2];
		sites[s].indices = calloc(pps[labels[i]].num_cart_gridpts, sizeof(int));
		CHECK_ALLOCATION(sites[s].indices);
		sites[s].projs = (real_proj_t*) malloc(sites[s].total_projs * sizeof(real_proj_t));
		int p = 0;
		sites[s].paths = malloc(3*pps[labels[i]].num_cart_gridpts * sizeof(double));
		CHECK_ALLOCATION(sites[s].paths);
		for (int j = 0; j < sites[s].num_projs; j++) {
			for (int m = -pps[labels[i]].funcs[j].l; m <= pps[labels[i]].funcs[j].l; m++) {
				sites[s].projs[p].l = pps[labels[i]].funcs[j].l;
				sites[s].projs[p].m = m;
				sites[s].projs[p].func_num = j;
				sites[s].projs[p].values = malloc(pps[labels[i]].num_cart_gridpts * sizeof(double complex));
				CHECK_ALLOCATION(sites[s].projs[p].values);
				p++;
			}
		}
	}

	//#pragma omp parallel for
	for (int s = 0; s < num_sites; s++) {
		int p = site_nums[s];
		double res[3] = {0,0,0};
		double frac[3] = {0,0,0};
		double testcoord[3] = {0,0,0};
		vcross(res, lattice+3, lattice+6);
		int grid1 = (int) (mag(res) * sites[s].rmax / vol * fftg[0]) + 1;
		vcross(res, lattice+0, lattice+6);
		int grid2 = (int) (mag(res) * sites[s].rmax / vol * fftg[1]) + 1;
		vcross(res, lattice+0, lattice+3);
		int grid3 = (int) (mag(res) * sites[s].rmax / vol * fftg[2]) + 1;
		int center1 = (int) round(coords[3*p+0] * fftg[0]);
		int center2 = (int) round(coords[3*p+1] * fftg[1]);
		int center3 = (int) round(coords[3*p+2] * fftg[2]);
		int ii=0, jj=0, kk=0;
		double R0 = (pps[labels[p]].proj_gridsize-1) * sites[s].rmax
			/ pps[labels[s]].proj_gridsize;
		for (int i = -grid1 + center1; i <= grid1 + center1; i++) {
			for (int j = -grid2 + center2; j <= grid2 + center2; j++) {
				for (int k = -grid3 + center3; k <= grid3 + center3; k++) {
					testcoord[0] = (double) i / fftg[0] - coords[3*p+0];
					testcoord[1] = (double) j / fftg[1] - coords[3*p+1];
					testcoord[2] = (double) k / fftg[2] - coords[3*p+2];
					frac_to_cartesian(testcoord, lattice);
					if (mag(testcoord) < R0) {
						ii = (i%fftg[0] + fftg[0]) % fftg[0];
						jj = (j%fftg[1] + fftg[1]) % fftg[1];
						kk = (k%fftg[2] + fftg[2]) % fftg[2];
						frac[0] = (double) ii / fftg[0];
						frac[1] = (double) jj / fftg[1];
						frac[2] = (double) kk / fftg[2];
						sites[s].indices[sites[s].num_indices] = ii*fftg[1]*fftg[2] + jj*fftg[2] + kk;
						sites[s].paths[3*sites[s].num_indices+0] = testcoord[0];
						sites[s].paths[3*sites[s].num_indices+1] = testcoord[1];
						sites[s].paths[3*sites[s].num_indices+2] = testcoord[2];
						for (int n = 0; n < sites[s].total_projs; n++) {
							if (pr0_pw1)
								sites[s].projs[n].values[sites[s].num_indices] = smooth_wave_value(
									pps[labels[p]].funcs[sites[s].projs[n].func_num],
									pps[labels[p]].smooth_grid, sites[s].projs[n].m, sites[s].rmax,
									pps[labels[p]].proj_gridsize, coords+3*p, frac, lattice);
							else
								sites[s].projs[n].values[sites[s].num_indices] = proj_value(
									pps[labels[p]].funcs[sites[s].projs[n].func_num],
									pps[labels[p]].proj_grid, sites[s].projs[n].m, sites[s].rmax,
									pps[labels[p]].proj_gridsize, coords+3*p, frac, lattice);
						}
						sites[s].num_indices++;
						//if (sites[s].num_indices >= pps[labels[p]].num_cart_gridpts)
						//	printf("SETUP_SITE ERROR %d %d %d\n", sites[s].num_indices, p, pr0_pw1);
					}
				}
			}
		}
	}
}

//adapted from VASP source code
double** spline_coeff(double* x, double* y, int N) {
	double** coeff = (double**) malloc(3 * sizeof(double*));
	CHECK_ALLOCATION(coeff);
	coeff[0] = (double*) malloc(N * sizeof(double));
	coeff[1] = (double*) malloc(N * sizeof(double));
	coeff[2] = (double*) malloc(N * sizeof(double));
	CHECK_ALLOCATION(coeff[0]);
	CHECK_ALLOCATION(coeff[1]);
	CHECK_ALLOCATION(coeff[2]);

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

double spline_integral(double* x, double* a, double** s, int size) {
	double* b = s[0];
	double* c = s[1];
	double* d = s[2];
	double dx = 0;
	double integral = 0;

	for (int i = 0; i < size - 1; i++) {
		dx = x[i+1] - x[i];
		integral += dx * (a[i] + dx * (b[i]/2 + dx * (c[i]/3 + d[i]*dx/4)));
	}

	return integral;

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
	else {
		printf("ERROR: sph_bessel l too high");
		return 0;
	}
}

double sbf(double x, int l) {
	if (x < 10e-6) {
		if (l==0) return 1;
		else return 0;
	}
	double jlm1 = sin(x) / x;
	double jl = sin(x) / (x*x) - cos(x) / x;
	double jlp1 = 0;
	if (l == 0) return jlm1;
	if (l == 1) return jl;
	for (int ll = 1; ll < l; ll++) {
		jlp1 = (2*ll+1)/x*jl - jlm1;
		jlm1 = jl;
		jl = jlp1;
	}
	return jlp1;
}

pswf_t* expand_symm_wf(pswf_t* rwf, int num_kpts, int* maps,
	double* ops, double* drs, double* kws, int* trs) {

	double* lattice = rwf->lattice;
	double* reclattice = rwf->reclattice;

	pswf_t* wf = (pswf_t*) malloc(sizeof(pswf_t));

	wf->num_elems = rwf->num_elems;
	wf->num_sites = rwf->num_sites;
	wf->pps = NULL;
	wf->G_bounds = (int*) malloc(6*sizeof(int));
	for (int i = 0; i < 6; i++) {
		wf->G_bounds[i] = 0;//rwf->G_bounds[i];
	}

	wf->kpts = (kpoint_t**) malloc(num_kpts * rwf->nspin * sizeof(kpoint_t*));
	wf->nspin = rwf->nspin;
	wf->nband = rwf->nband;
	wf->nwk = num_kpts;
	wf->lattice = (double*) malloc(9*sizeof(double));
	wf->reclattice = (double*) malloc(9*sizeof(double));
	for (int i = 0; i < 9; i++) {
		wf->lattice[i] = lattice[i];
		wf->reclattice[i] = reclattice[i];
	}
	wf->fftg = NULL;

	wf->is_ncl = rwf->is_ncl;

	wf->num_aug_overlap_sites = 0;
	wf->dcoords = NULL;
	wf->overlaps = NULL;
	wf->num_projs = NULL;
	wf->wp_num = 0;

	//#pragma omp parallel for
	for (int knum = 0; knum < num_kpts * wf->nspin; knum++) {
		wf->kpts[knum] = (kpoint_t*) malloc(sizeof(kpoint_t));

		double pw[3] = {0,0,0};
		int rnum = maps[knum%num_kpts];
		int tr = trs[knum%num_kpts];
		if (knum >= num_kpts && rwf->nspin ==2) {
			rnum += rwf->nwk;
		}

		kpoint_t* kpt = wf->kpts[knum];
		kpoint_t* rkpt = rwf->kpts[rnum];

		kpt->up = rkpt->up;
		kpt->num_waves = rkpt->num_waves;
		kpt->k = (double*) malloc(3 * sizeof(double));
		rotation_transform(kpt->k, ops+OPSIZE*(knum%num_kpts), rkpt->k);
		if (tr == 1) {
			kpt->k[0] *= -1;
			kpt->k[1] *= -1;
			kpt->k[2] *= -1;
		}
		double kdiff[3] = {round(kpt->k[0]), round(kpt->k[1]), round(kpt->k[2])};
		kpt->k[0] -= kdiff[0];
		kpt->k[1] -= kdiff[1];
		kpt->k[2] -= kdiff[2];
		//printf("OLD KPT %lf %lf %lf\n", rkpt->k[0], rkpt->k[1], rkpt->k[2]);
		//printf("NEW KPT %lf %lf %lf\n", kpt->k[0], kpt->k[1], kpt->k[2]);
		//kpt->Gs = (int*) malloc(3 * kpt->num_waves * sizeof(int));

		kpt->weight = kws[knum%num_kpts];
		kpt->num_bands = rkpt->num_bands;
		kpt->bands = (band_t**) malloc(kpt->num_bands * sizeof(band_t*));
		kpt->expansion = NULL;

		int* igall = malloc(3*kpt->num_waves*sizeof(int));
		if (igall == NULL) {
		    	ALLOCATION_FAILED();
		}
		int nb1max = rwf->G_bounds[1] - rwf->G_bounds[0] + 2;
		int nb2max = rwf->G_bounds[3] - rwf->G_bounds[2] + 2;
		int nb3max = rwf->G_bounds[5] - rwf->G_bounds[4] + 2;
		double encut = rwf->encut;
		double* b1 = reclattice;
		double* b2 = reclattice+3;
		double* b3 = reclattice+6;
		int ncnt = -1;
		for (int ig3 = 0; ig3 <= 2 * nb3max; ig3++) {
			int ig3p = ig3;
			if (ig3 > nb3max) ig3p = ig3 - 2 * nb3max - 1;
			for (int ig2 = 0; ig2 <= 2 * nb2max; ig2++) {
				int ig2p = ig2;
				if (ig2 > nb2max) ig2p = ig2 - 2 * nb2max - 1;
				for (int ig1 = 0; ig1 <= 2 * nb1max; ig1++) {
					int ig1p = ig1;
					if (ig1 > nb1max) ig1p = ig1 - 2 * nb1max - 1;
					double sumkg[3];
					for (int j = 0; j < 3; j++) {
						sumkg[j] = (kpt->k[0]+ig1p) * b1[j]
									+ (kpt->k[1]+ig2p) * b2[j]
									+ (kpt->k[2]+ig3p) * b3[j];
					}
					double gtot = mag(sumkg);
					double etot = pow(gtot,2.0) / CCONST;
					//printf("%lf %lf\n", etot, gtot);
					if (etot <= encut) {
						ncnt++;
						igall[ncnt*3+0] = ig1p;
						igall[ncnt*3+1] = ig2p;
						igall[ncnt*3+2] = ig3p;
						if (ig1p < wf->G_bounds[0]) wf->G_bounds[0] = ig1p;
						else if (ig1p > wf->G_bounds[1]) wf->G_bounds[1] = ig1p;
						if (ig2p < wf->G_bounds[2]) wf->G_bounds[2] = ig2p;
						else if (ig2p > wf->G_bounds[3]) wf->G_bounds[3] = ig2p;
						if (ig3p < wf->G_bounds[4]) wf->G_bounds[4] = ig3p;
						else if (ig3p > wf->G_bounds[5]) wf->G_bounds[5] = ig3p;
					}
				}
			}
		}
		ncnt++;

		if (ncnt * 2 == rkpt->num_waves) {
			printf("This is an NCL wavefunction!\n");
			wf->is_ncl = 1;
			for (int iplane = 0; iplane < rkpt->num_waves/2; iplane++) {
				igall[3*(rkpt->num_waves/2+iplane)+0] = igall[3*iplane+0];
				igall[3*(rkpt->num_waves/2+iplane)+1] = igall[3*iplane+1];
				igall[3*(rkpt->num_waves/2+iplane)+2] = igall[3*iplane+2];
			}
		} else if (ncnt != rkpt->num_waves) {
			printf("ERROR %d %d %lf %lf %lf %lf\n", ncnt, kpt->num_waves,
				kpt->k[0], kpt->k[1], kpt->k[2], CCONST);
		}
		kpt->Gs = igall;

		int ngx = wf->G_bounds[1] - wf->G_bounds[0] + 1;
		int gxmin = wf->G_bounds[0];
		int ngy = wf->G_bounds[3] - wf->G_bounds[2] + 1;
		int gymin = wf->G_bounds[2];
		int ngz = wf->G_bounds[5] - wf->G_bounds[4] + 1;
		int gzmin = wf->G_bounds[4];
		int* kptinds = (int*) malloc(ngx*ngy*ngz * sizeof(int));
		for (int w = 0; w < ngx*ngy*ngz; w++) kptinds[w] = -1;

		int gx, gy, gz;
		int* gmaps = (int*) malloc(kpt->num_waves * sizeof(int));
		float complex* factors = (float complex*) malloc(
			kpt->num_waves * sizeof(float complex));
		for (int w = 0; w < kpt->num_waves; w++) gmaps[w] = -1;
		for (int w = 0; w < rkpt->num_waves; w++) {
			gx = kpt->Gs[3*w+0];
			gy = kpt->Gs[3*w+1];
			gz = kpt->Gs[3*w+2];
			kptinds[(gx-gxmin)*ngy*ngz + (gy-gymin)*ngz + (gz-gzmin)] = w;
		}

		double* dr = drs + 3 * (knum%num_kpts);
		double* op = ops+OPSIZE*(knum%num_kpts);

		for (int g = 0; g < kpt->num_waves; g++) {

			//pw[0] = rkpt->k[0] + rkpt->Gs[3*g+0];
			//pw[1] = rkpt->k[1] + rkpt->Gs[3*g+1];
			//pw[2] = rkpt->k[2] + rkpt->Gs[3*g+2];
			pw[0] = rkpt->Gs[3*g+0];
			pw[1] = rkpt->Gs[3*g+1];
			pw[2] = rkpt->Gs[3*g+2];

			rotation_transform(pw, ops+OPSIZE*(knum%num_kpts), pw);
			if (tr == 1) {
				pw[0] *= -1;
				pw[1] *= -1;
				pw[2] *= -1;
			}

			//pw[0] -= kpt->k[0];
			//pw[1] -= kpt->k[1];
			//pw[2] -= kpt->k[2];
			pw[0] += kdiff[0];
			pw[1] += kdiff[1];
			pw[2] += kdiff[2];

			gx = (int) round(pw[0]);
			gy = (int) round(pw[1]);
			gz = (int) round(pw[2]);
			gmaps[kptinds[(gx-gxmin)*ngy*ngz + (gy-gymin)*ngz + (gz-gzmin)]] = g;
			if (tr == 0) {
				factors[kptinds[(gx-gxmin)*ngy*ngz + (gy-gymin)*ngz + (gz-gzmin)]] = cexpf(
					-I * 2 * PI * (dot(kpt->k, dr) + dot(pw, dr)) );
			} else {
				factors[kptinds[(gx-gxmin)*ngy*ngz + (gy-gymin)*ngz + (gz-gzmin)]] = cexpf(
					I * 2 * PI * (dot(kpt->k, dr) + dot(pw, dr)) );
			}

			if (kptinds[(gx-gxmin)*ngy*ngz + (gy-gymin)*ngz + (gz-gzmin)] < 0) {
				printf("ERROR, BAD PLANE WAVE MAPPING %d %d %d %d %d %d %lf %lf %lf\n %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
						rkpt->Gs[3*g+0], rkpt->Gs[3*g+1], rkpt->Gs[3*g+2], gx, gy, gz,
						pw[0], pw[1], pw[2], op[0],op[1],op[2],op[3],op[4],op[5],op[6],op[7],op[8]);
			}

			//printf("OLD G %d %d %d\n", okpt->Gs[3*g+0], okpt->Gs[3*g+1], okpt->Gs[3*g+2]);
			//printf("NEW G %d %d %d\n", kpt->Gs[3*g+0],  kpt->Gs[3*g+1],  kpt->Gs[3*g+2]);
		}

		for (int b = 0; b < kpt->num_bands; b++) {
			kpt->bands[b] = (band_t*) malloc(sizeof(band_t));
			kpt->bands[b]->n = rkpt->bands[b]->n;
			kpt->bands[b]->num_waves = rkpt->bands[b]->num_waves;
			kpt->bands[b]->occ = rkpt->bands[b]->occ;
			kpt->bands[b]->energy = rkpt->bands[b]->energy;
			kpt->bands[b]->Cs = (float complex*) malloc(
				kpt->num_waves * sizeof(float complex));
			kpt->bands[b]->CRs = NULL;
			kpt->bands[b]->CAs = NULL;
			kpt->bands[b]->projections = NULL;
			kpt->bands[b]->up_projections = NULL;
			kpt->bands[b]->down_projections = NULL;
			kpt->bands[b]->wave_projections = NULL;
			//double total = 0;
			for (int w = 0; w < kpt->num_waves; w++) {
				if (gmaps[w] < 0) {
					printf("ERROR, INCOMPLETE PLANE WAVE MAPPING\n");
				}
				kpt->bands[b]->Cs[w] = factors[w] * rkpt->bands[b]->Cs[gmaps[w]];
				if (tr == 1) {
					kpt->bands[b]->Cs[w] = conj(kpt->bands[b]->Cs[w]);
				}
				//total += cabs(cabs(kpt->bands[b]->Cs[w]) - cabs(okpt->bands[b]->Cs[w]));
			}
			//printf("energies %lf %lf %e\n", kpt->bands[b]->energy, okpt->bands[b]->energy, total);
		}
	
		free(gmaps);
		free(kptinds);
		free(factors);
	}

	wf->encut = rwf->encut;
	return wf;
}

void CHECK_ALLOCATION(void* ptr) {
	if (ptr == NULL) {
		ALLOCATION_FAILED();
	}
}

void ALLOCATION_FAILED(void) {
	printf("ALLOCATION FAILED\n");
	exit(-1);
}

void CHECK_STATUS(int status) {
	if (status != 0) {
		printf("ROUTINE FAILED WITH STATUS %d:\n", status);
		char* message = DftiErrorMessage(status);
		printf("%s\n", message);
		exit(-1);
	}
}

double* get_smooth_wave(ppot_t* lst, int num) {
	//return lst[num].funcs[0].smooth_diffwave;
	return lst[num].smooth_grid;
}

