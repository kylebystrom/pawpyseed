#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <omp.h>
#include <time.h>
#include "utils.h"
#include "projector.h"
#include <mkl.h>
#include <mkl_types.h>
#include "linalg.h"
#include "quadrature.h"
#include "radial.h"
#include "sbt.h"

#define c 0.262465831
#define PI 3.14159265358979323846

ppot_t* get_projector_list(int num_els, int* labels, int* ls, double* proj_grids, double* wave_grids,
	double* projectors, double* aewaves, double* pswaves, double* rmaxs) {

	setbuf(stdout,NULL);	
	ppot_t* pps = (ppot_t*) malloc(num_els * sizeof(ppot_t));
	CHECK_ALLOCATION(pps);
	int wt = 0;
	int pt = 0;
	int wgt = 0;
	int pgt = 0;
	int l_num = 0;
	for (int i = 0; i < num_els; i++) {
		pps[i].num_projs = labels[4*i+1];
		pps[i].rmax = rmaxs[i];
		//sscanf(rmaxs[i], "%lf", &(pps[i].rmax));
		pps[i].proj_gridsize = labels[4*i+2];
		pps[i].wave_gridsize = labels[4*i+3];
		printf("vals %d %d %d\n", pps[i].num_projs, pps[i].proj_gridsize, pps[i].wave_gridsize);
		pps[i].total_projs = 0;
		pps[i].wave_grid = (double*) malloc((pps[i].wave_gridsize)*sizeof(double));
		pps[i].kwave_grid = (double*) malloc((pps[i].wave_gridsize)*sizeof(double));
		CHECK_ALLOCATION(pps[i].wave_grid);
		CHECK_ALLOCATION(pps[i].kwave_grid);
		pps[i].lmax = 0;
		pps[i].pspw_overlap_matrix = NULL;
		pps[i].aepw_overlap_matrix = NULL;
		pps[i].diff_overlap_matrix = NULL;
		for (int j = 0; j < pps[i].wave_gridsize; j++) {
			pps[i].wave_grid[j] = wave_grids[wgt];
			wgt++;
		}
		pps[i].proj_grid = (double*) malloc(pps[i].proj_gridsize*sizeof(double));
		CHECK_ALLOCATION(pps[i].proj_grid);
		for (int j = 0; j < pps[i].proj_gridsize; j++) {
			pps[i].proj_grid[j] = pps[i].rmax / pps[i].proj_gridsize * j;
			pgt++;
		}
		funcset_t* funcs = (funcset_t*) malloc(pps[i].num_projs*sizeof(funcset_t));
		CHECK_ALLOCATION(funcs);
		for (int k = 0; k < pps[i].num_projs; k++) {
			funcs[k].proj = (double*) malloc(sizeof(double)*pps[i].proj_gridsize);
			funcs[k].aewave = (double*) malloc(sizeof(double)*pps[i].wave_gridsize);
			funcs[k].pswave = (double*) malloc(sizeof(double)*pps[i].wave_gridsize);
			funcs[k].diffwave = (double*) malloc(sizeof(double)*pps[i].wave_gridsize);
			CHECK_ALLOCATION(funcs[k].proj);
			CHECK_ALLOCATION(funcs[k].aewave);
			CHECK_ALLOCATION(funcs[k].pswave);
			CHECK_ALLOCATION(funcs[k].diffwave);
			funcs[k].l = ls[l_num];
			if (funcs[k].l > pps[i].lmax)
				pps[i].lmax = funcs[k].l;
			pps[i].total_projs += 2 * ls[l_num] + 1;
			l_num++;
			for (int j = 0; j < pps[i].wave_gridsize; j++) {
				funcs[k].aewave[j] = aewaves[wt];
				funcs[k].pswave[j] = pswaves[wt];
				funcs[k].diffwave[j] = aewaves[wt] - pswaves[wt];
				wt++;
			}
			for (int j = 0; j < pps[i].proj_gridsize; j++) {
				funcs[k].proj[j] = projectors[pt];
				pt++;
			}
			funcs[k].proj_spline = spline_coeff(pps[i].proj_grid, funcs[k].proj, pps[i].proj_gridsize);
			funcs[k].aewave_spline = spline_coeff(pps[i].wave_grid, funcs[k].aewave, pps[i].wave_gridsize);
			funcs[k].pswave_spline = spline_coeff(pps[i].wave_grid, funcs[k].pswave, pps[i].wave_gridsize);
			funcs[k].diffwave_spline = spline_coeff(pps[i].wave_grid, funcs[k].diffwave, pps[i].wave_gridsize);
		}
		sbt_descriptor_t* d = spherical_bessel_transform_setup(520.0, 5000.0, pps[i].lmax,
			pps[i].wave_gridsize, pps[i].wave_grid, pps[i].kwave_grid);
		for (int k = 0; k < pps[i].num_projs; k++) {
			funcs[k].kwave = wave_spherical_bessel_transform(d, funcs[k].diffwave, funcs[k].l);
			//funcs[k].kwave = besselt(pps[i].wave_grid, pps[i].kwave_grid, funcs[k].diffwave, 520.0, pps[i].wave_gridsize, funcs[k].l);
			funcs[k].kwave_spline = spline_coeff(pps[i].kwave_grid, funcs[k].kwave, pps[i].wave_gridsize);
		}
		for (int l = 0; l <= pps[i].lmax; l++) {
			free(d->mult_table[l]);
		}
		free(d->ks);
		free(d->rs);
		free(d->mult_table);
		free(d);
		pps[i].funcs = funcs;
		make_pwave_overlap_matrices(pps+i);
	}
	mkl_free_buffers();
	printf("finished making projector list\n");
	return pps;
}

double* besselt(double* r, double* k, double* f, double encut, int N, int l) {
	
	double kmax = pow(encut*c, 0.5);
	double drho = log(r[1]/r[0]);
	double kmin = kmax * exp((1-N)*drho);
	double* g = (double*) malloc(N*sizeof(double));
	CHECK_ALLOCATION(g);
	for (int i = 0; i < N; i++) {
		double dr = r[0];
		g[i] = 0;
		k[i] = kmin * exp(i*drho);
		for (int j = 0; j < N; j++) {
			g[i] += f[j] * r[j] * sph_bessel(k[i], r[j], l) * dr; 
			if (j != N-1) dr = r[j+1]-r[j];
		}
	}
	return g;

}

real_proj_site_t* projector_values(int num_sites, int* labels, double* coords,
	double* lattice, double* reclattice, ppot_t* pps, int* fftg) {

	double intervals[3] = {mag(lattice)/fftg[0], mag(lattice+3)/fftg[1], mag(lattice+6)/fftg[2]};
	double vol = determinant(lattice);
	int num_pts = fftg[0] * fftg[1] * fftg[2];

	real_proj_site_t* sites = (real_proj_site_t*) malloc(num_sites * sizeof(real_proj_site_t));
	CHECK_ALLOCATION(sites);
	for (int i = 0; i < num_sites; i++) {
		sites[i].index = i;
		sites[i].elem = labels[i];
		sites[i].num_projs = pps[labels[i]].num_projs;
		sites[i].rmax = pps[labels[i]].rmax;
		sites[i].total_projs = pps[labels[i]].total_projs;
		sites[i].num_indices = 0;
		sites[i].coord = malloc(3 * sizeof(double));
		CHECK_ALLOCATION(sites[i].coord);
		sites[i].coord[0] = coords[3*i+0];
		sites[i].coord[1] = coords[3*i+1];
		sites[i].coord[2] = coords[3*i+2];
		sites[i].indices = calloc(pps[labels[i]].num_cart_gridpts, sizeof(int));
		CHECK_ALLOCATION(sites[i].indices);
		sites[i].projs = (real_proj_t*) malloc(sites[i].total_projs * sizeof(real_proj_t));
		int p = 0;
		for (int j = 0; j < sites[i].num_projs; j++) {
			for (int m = -pps[labels[i]].funcs[j].l; m <= pps[labels[i]].funcs[j].l; m++) {
				sites[i].projs[p].l = pps[labels[i]].funcs[j].l;
				sites[i].projs[p].m = m;
				sites[i].projs[p].func_num = j;
				sites[i].projs[p].values = malloc(pps[labels[i]].num_cart_gridpts * sizeof(double complex));
				sites[i].projs[p].paths = malloc(3*pps[labels[i]].num_cart_gridpts * sizeof(double));
				CHECK_ALLOCATION(sites[i].projs[p].values);
				CHECK_ALLOCATION(sites[i].projs[p].paths);
				p++;
			}
		}
	}

	for (int i = 0; i < fftg[0]; i++) {
		double path[3] = {0,0,0};
		double r = 0;
		for (int j = 0; j < fftg[1]; j++) {
			for (int k = 0; k  < fftg[2]; k++) {
				double frac[3] = {(double)i/fftg[0], (double)j/fftg[1], (double)k/fftg[2]};
				for (int p = 0; p < num_sites; p++) {
					min_cart_path(frac, coords+3*p, lattice, path, &r);
					if (r < 0.99 * sites[p].rmax) {
						sites[p].indices[sites[p].num_indices] = i*fftg[1]*fftg[2] + j*fftg[2] + k;
						for (int n = 0; n < sites[p].total_projs; n++) {
							sites[p].projs[n].paths[3*sites[p].num_indices] = path[0];
							sites[p].projs[n].paths[3*sites[p].num_indices+1] = path[1];
							sites[p].projs[n].paths[3*sites[p].num_indices+2] = path[2];
							sites[p].projs[n].values[sites[p].num_indices] = proj_value(pps[labels[p]].funcs[sites[p].projs[n].func_num],
								pps[labels[p]].proj_grid,
								sites[p].projs[n].m, sites[p].rmax, coords+3*p, frac, lattice);
						//if (p == 0) printf("coordandval %lf %lf %lf %d %d\n", r, creal(sites[p].projs[n].values[sites[p].num_indices])*pow(vol,0.5), cimag(sites[p].projs[n].values[sites[p].num_indices])*pow(vol, 0.5), sites[p].num_indices,n);
						}
						sites[p].num_indices++;
					}
				}
			}
		}
	}
	//for (int i = 0; i < num_sites; i++) {
	//	printf("looking for nan %d %e\n", sites[0].num_indices, creal(sites[i].projs[0].values[0]));
	//}

	return sites;
}

void onto_projector_helper(band_t* band, fft_complex* x, real_proj_site_t* sites,
	int num_sites, int* labels, double* lattice, double* reclattice, double* kpt, ppot_t* pps, int* fftg) {

	double dv = determinant(lattice) / fftg[0] / fftg[1] / fftg[2];

	band->projections = (projection_t*) malloc(num_sites * sizeof(projection_t));
	CHECK_ALLOCATION (band->projections);

	double path[3] = {0,0,0};
	double kdotr = 0;

	double kpt_cart[3] = {0,0,0};
	kpt_cart[0] = kpt[0];
	kpt_cart[1] = kpt[1];
	kpt_cart[2] = kpt[2];
	frac_to_cartesian(kpt_cart, reclattice);

	int num_indices, index, t=0;
	for (int s = 0; s < num_sites; s++) {
		//printf("checking coord %lf %lf %lf\n", sites[s].coord[0], sites[s].coord[1], sites[s].coord[2]);
		num_indices = sites[s].num_indices;
		int* indices = sites[s].indices;
		projection_t* projections = band->projections;
		projections[s].num_projs = sites[s].num_projs;
		projections[s].total_projs = sites[s].total_projs;
		projections[s].ns = malloc(sites[s].total_projs * sizeof(int));
		projections[s].ls = malloc(sites[s].total_projs * sizeof(int));
		projections[s].ms = malloc(sites[s].total_projs * sizeof(int));
		projections[s].overlaps = (double complex*) malloc(sites[s].total_projs * sizeof(double complex));
		CHECK_ALLOCATION(projections[s].ns);
		CHECK_ALLOCATION(projections[s].ls);
		CHECK_ALLOCATION(projections[s].ms);
		CHECK_ALLOCATION(projections[s].overlaps);
		for (int p = 0; p < sites[s].total_projs; p++) {
			projections[s].ns[p] = sites[s].projs[p].func_num;
			projections[s].ls[p] = sites[s].projs[p].l;
			projections[s].ms[p] = sites[s].projs[p].m;
			double complex* values = sites[s].projs[p].values;
			double complex total = 0 + 0 * I;
			for (int i = 0; i < num_indices; i++) {
				index = indices[i];
				kdotr = dot(kpt_cart, sites[s].projs[p].paths+i*3);
				total += fft_mult(index, x, conj(values[i]) * dv * cexp(I * kdotr));
			}
			projections[s].overlaps[p] = total;
		}
	}
}

void onto_projector(kpoint_t* kpt, int band_num, real_proj_site_t* sites, int num_sites, int* labels,
	int* G_bounds, double* lattice, double* reclattice, ppot_t* pps, int* fftg) {

	double* k = kpt->k;
	int* Gs = kpt->Gs;
	float complex* Cs = kpt->bands[band_num]->Cs;
	int num_waves = kpt->num_waves;
	
	fft_complex* x = (fft_complex*) mkl_calloc(fftg[0]*fftg[1]*fftg[2], sizeof(fft_complex), 64);
	CHECK_ALLOCATION(x);
	//printf("integrating params %e %e %e %e %e\n", dv, inv_sqrt_vol, kmins[0], kmins[1], kmins[2]);
	//printf("determinant %lf\n", determinant(lattice));
	fft3d(x, G_bounds, lattice, k, Gs, Cs, num_waves, fftg);

	onto_projector_helper(kpt->bands[band_num], x, sites, num_sites, labels, lattice, reclattice, k, pps, fftg);

	mkl_free(x);
}

void add_num_cart_gridpts(ppot_t* pp_ptr, double* lattice, int* fftg) {

	ppot_t pp = *pp_ptr;

	double maga1 = mag(lattice+0);
	double maga2 = mag(lattice+3);
	double maga3 = mag(lattice+6);

	double vtemp[3];
	double vmag, sinphi123;
	
	double phi12 = acos(dot(lattice+0, lattice+3) / (maga1 * maga2));
	vcross(vtemp, lattice+0, lattice+3);
	vmag = mag(vtemp);
	sinphi123 = dot(lattice+6, vtemp) / (vmag * maga3);
	double na1maxA = pp.rmax * fftg[0] / (maga1 * fabs(sin(phi12))) + 1;
	double na2maxA = pp.rmax * fftg[1] / (maga2 * fabs(sin(phi12))) + 1;
	double na3maxA = pp.rmax * fftg[2] / (maga3 * fabs(sinphi123)) + 1;
	int npmaxA = (int)(4.0/3.0*PI*na1maxA*na2maxA*na3maxA) + 1;

	double phi13 = acos(dot(lattice+0, lattice+6) / (maga1 * maga3));
	vcross(vtemp, lattice+0, lattice+6);
	vmag = mag(vtemp);
	sinphi123 = dot(lattice+3, vtemp) / (vmag * maga2);
	double na1maxB = pp.rmax * fftg[0] / (maga1 * fabs(sin(phi13))) + 1;
	double na2maxB = pp.rmax * fftg[1] / (maga2 * fabs(sinphi123)) + 1;
	double na3maxB = pp.rmax * fftg[2] / (maga3 * fabs(sin(phi13))) + 1;
	int npmaxB = (int)(4.0/3.0*PI*na1maxB*na2maxB*na3maxB) + 1;

	double phi23 = acos(dot(lattice+3, lattice+6) / (maga2 * maga3));
	vcross(vtemp, lattice+3, lattice+6);
	vmag = mag(vtemp);
	sinphi123 = dot(lattice, vtemp) / (vmag * maga1);
	double na1maxC = pp.rmax * fftg[0] / (maga1 * fabs(sinphi123)) + 1;
	double na2maxC = pp.rmax * fftg[1] / (maga2 * fabs(sin(phi23))) + 1;
	double na3maxC = pp.rmax * fftg[2] / (maga3 * fabs(sin(phi23))) + 1;
	int npmaxC = (int)(4.0/3.0*PI*na1maxC*na2maxC*na3maxC) + 1;

	printf("ancg %lf %lf %lf %lf %lf %d %d %d\n", maga1, maga2, maga3, na2maxA, na3maxA, npmaxA, npmaxB, npmaxC);

	int npmax = npmaxA;
	if (npmaxB > npmax) npmax = npmaxB;
	if (npmaxC > npmax) npmax = npmaxC;

	pp_ptr->num_cart_gridpts = npmax;
}

void make_pwave_overlap_matrices(ppot_t* pp_ptr) {
	ppot_t pp = *pp_ptr;
	int size = pp.num_projs * pp.num_projs;
	double* psov = (double*) calloc(size, sizeof(double));
	double* aeov = (double*) calloc(size, sizeof(double));
	double* diov = (double*) calloc(size, sizeof(double));
	CHECK_ALLOCATION(psov);
	CHECK_ALLOCATION(aeov);
	CHECK_ALLOCATION(diov);

	for (int i = 0; i < pp.num_projs; i++) {
		for (int j = i; j < pp.num_projs; j++) {
			if (pp.funcs[i].l == pp.funcs[j].l) {
				double* ps1 = pp.funcs[i].pswave;
				double* ps2 = pp.funcs[j].pswave;
				double* ae1 = pp.funcs[i].aewave;
				double* ae2 = pp.funcs[j].aewave;
				double* psprod = (double*) malloc(pp.wave_gridsize*sizeof(double));
				double* aeprod = (double*) malloc(pp.wave_gridsize*sizeof(double));
				double* diprod = (double*) malloc(pp.wave_gridsize*sizeof(double));
				for (int k = 0; k < pp.wave_gridsize; k++) {
					psprod[k] = ps1[k] * ps2[k];
					aeprod[k] = ae1[k] * ae2[k];
					diprod[k] = (ae1[k]-ps1[k]) * (ae2[k]-ps2[k]);
				}
				double** psspline = spline_coeff(pp.wave_grid, psprod, pp.wave_gridsize);
				double** aespline = spline_coeff(pp.wave_grid, aeprod, pp.wave_gridsize);
				double** displine = spline_coeff(pp.wave_grid, diprod, pp.wave_gridsize);
				psov[i*pp.num_projs+j] = spline_integral(pp.wave_grid, psprod, psspline, pp.wave_gridsize);
				aeov[i*pp.num_projs+j] = spline_integral(pp.wave_grid, aeprod, aespline, pp.wave_gridsize);
				diov[i*pp.num_projs+j] = spline_integral(pp.wave_grid, diprod, displine, pp.wave_gridsize);
			}
		}
	}
	for (int i = 1; i < pp.num_projs; i++) {
		for (int j = 0; j < i; j++) {
			psov[pp.num_projs*i+j] = psov[pp.num_projs*j+i];
			aeov[pp.num_projs*i+j] = aeov[pp.num_projs*j+i];
			diov[pp.num_projs*i+j] = diov[pp.num_projs*j+i];
		}
	}

	pp_ptr->pspw_overlap_matrix = psov;
	pp_ptr->aepw_overlap_matrix = aeov;
	pp_ptr->diff_overlap_matrix = diov;
}

void setup_projections(pswf_t* wf, ppot_t* pps, int num_elems,
		int num_sites, int* fftg, int* labels, double* coords) {

	printf("ptrvals %p %p\n", wf, pps);
	wf->num_elems = num_elems;
	wf->pps = pps;
	printf("started setup_proj\n");
	#pragma omp parallel for 
	for (int p = 0; p < num_elems; p++) {
		add_num_cart_gridpts(pps+p, wf->lattice, fftg);
	}
	int NUM_KPTS = wf->nwk * wf->nspin;
	int NUM_BANDS = wf->nband;
	printf("calculating projector_values\n");
	real_proj_site_t* sites = projector_values(num_sites, labels, coords,
		wf->lattice, wf->reclattice, pps, fftg);
	printf("onto_projector calcs\n");
	#pragma omp parallel for 
	for (int w = 0; w < NUM_BANDS * NUM_KPTS; w++) {
		kpoint_t* kpt = wf->kpts[w % NUM_KPTS];
		int band_num = w / NUM_KPTS;
		onto_projector(kpt, band_num, sites, num_sites, labels,
			wf->G_bounds, wf->lattice, wf->reclattice, pps, fftg);
	}
	free_real_proj_site_list(sites, num_sites);	
	generate_rayleigh_expansion_terms(wf, pps, num_elems);
}

void setup_projections_copy_rayleigh(pswf_t* wf, pswf_t* wf_R, ppot_t* pps, int num_elems,
		int num_sites, int* fftg, int* labels, double* coords) {

	printf("ptrvals %p %p\n", wf, pps);
	wf->num_elems = num_elems;
	wf->pps = pps;
	printf("started setup_proj\n");
	#pragma omp parallel for 
	for (int p = 0; p < num_elems; p++) {
		add_num_cart_gridpts(pps+p, wf->lattice, fftg);
	}
	int NUM_KPTS = wf->nwk * wf->nspin;
	int NUM_BANDS = wf->nband;
	printf("calculating projector_values\n");
	real_proj_site_t* sites = projector_values(num_sites, labels, coords,
		wf->lattice, wf->reclattice, pps, fftg);
	printf("onto_projector calcs\n");
	#pragma omp parallel for 
	for (int w = 0; w < NUM_BANDS * NUM_KPTS; w++) {
		kpoint_t* kpt = wf->kpts[w % NUM_KPTS];
		int band_num = w / NUM_KPTS;
		onto_projector(kpt, band_num, sites, num_sites, labels,
			wf->G_bounds, wf->lattice, wf->reclattice, pps, fftg);
	}
	free_real_proj_site_list(sites, num_sites);	
	copy_rayleigh_expansion_terms(wf, pps, num_elems, wf_R);
}

void overlap_setup(pswf_t* wf_R, pswf_t* wf_S, ppot_t* pps,
	int* labels_R, int* labels_S, double* coords_R, double* coords_S,
	int* N_R, int* N_S, int* N_RS_R, int* N_RS_S, int num_N_R, int num_N_S, int num_N_RS) {

	double complex** overlaps = (double complex**) malloc(num_N_RS * sizeof(double complex*));
	CHECK_ALLOCATION(overlaps);

	printf("STARTING OVERLAP_SETUP\n");
	int NUM_KPTS = wf_R->nwk * wf_R->nspin;
	int NUM_BANDS = wf_R->nband;
	double inv_sqrt_vol = pow(determinant(wf_R->lattice), -0.5);
	#pragma omp parallel for
	for (int w = 0; w < NUM_BANDS * NUM_KPTS; w++) {
		kpoint_t* kpt_R = wf_R->kpts[w%NUM_KPTS];
		kpoint_t* kpt_S = wf_S->kpts[w%NUM_KPTS];
		band_t* band_R = kpt_R->bands[w/NUM_KPTS];
		band_t* band_S = kpt_S->bands[w/NUM_KPTS];

		projection_t* wps = (projection_t*) malloc(num_N_R * sizeof(projection_t));
		CHECK_ALLOCATION(wps);
		for (int n = 0; n < num_N_R; n++) {
			int s = N_R[n];
			ppot_t pp = pps[labels_R[s]];
			wps[n].overlaps = (double complex*) malloc(pp.total_projs * sizeof(double complex));
			CHECK_ALLOCATION(wps[n].overlaps);
			int t = 0;
			for (int i = 0; i < pp.num_projs; i++) {
				for (int m = -pp.funcs[i].l; m <= pp.funcs[i].l; m++) {
					wps[n].overlaps[t] = rayexp(kpt_S->k, kpt_S->Gs,
						band_S->Cs, pp.funcs[i].l, m,
						kpt_S->num_waves,
						kpt_R->expansion[labels_R[s]][i].terms,
						coords_R + s*3) * inv_sqrt_vol;
					t++;
				}
			}
		}
		band_S->wave_projections = wps;
	}
	printf("ONE THIRD DONE\n");
	#pragma omp parallel for
	for (int w = 0; w < NUM_BANDS * NUM_KPTS; w++) {
		kpoint_t* kpt_R = wf_R->kpts[w%NUM_KPTS];
		kpoint_t* kpt_S = wf_S->kpts[w%NUM_KPTS];
		band_t* band_R = kpt_R->bands[w/NUM_KPTS];
		band_t* band_S = kpt_S->bands[w/NUM_KPTS];

		projection_t* wps = (projection_t*) malloc(num_N_S * sizeof(projection_t));
		for (int n = 0; n < num_N_S; n++) {
			int s = N_S[n];
			ppot_t pp = pps[labels_S[s]];
			wps[n].overlaps = (double complex*) malloc(pp.total_projs * sizeof(double complex));
			int t = 0;
			for (int i = 0; i < pp.num_projs; i++) {
				for (int m = -pp.funcs[i].l; m <= pp.funcs[i].l; m++) {
					wps[n].overlaps[t] = rayexp(kpt_R->k, kpt_R->Gs,
						band_R->Cs, pp.funcs[i].l, m,
						kpt_R->num_waves,
						kpt_S->expansion[labels_S[s]][i].terms,
						coords_S + s*3) * inv_sqrt_vol;
					t++;
				}
			}
		}
		band_R->wave_projections = wps;
	}
	printf("TWO THIRDS DONE\n");
	
	double* dcoords = malloc(3 * num_N_RS * sizeof(double));
	CHECK_ALLOCATION(dcoords);
	#pragma omp parallel for
	for (int i = 0; i < num_N_RS; i++) {
		double R = 0;
		int l1, l2;
		int s1 = N_RS_R[i];
		int s2 = N_RS_S[i];
		printf("%d %d\n", N_RS_R[i], N_RS_S[i]);
		ppot_t pp1 = pps[labels_R[s1]];
		ppot_t pp2 = pps[labels_S[s2]];
		// CALCULATE THE DIFF COORD HERE, PASS TO offsite_wave_overlap AND SAVE IT FOR USE IN compensation_terms
		overlaps[i] = calloc(pp1.total_projs * pp2.total_projs, sizeof(double complex));
		CHECK_ALLOCATION(overlaps[i]);
		double* coord1 = coords_R + 3 * s1;
		double* coord2 = coords_S + 3 * s2;
		min_cart_path(coord2, coord1, wf_R->lattice, dcoords + 3*i, &R);
		printf("herewego %lf %lf %lf\nherewego%lf %lf %lf\n", coord1[0], coord1[1], coord1[2], coord2[0], coord2[1], coord2[2]);
		int tj = 0;
		for (int j = 0; j < pp1.num_projs; j++) {
			l1 = pp1.funcs[j].l;
			for (int m1 = -l1; m1 <= l1; m1++) {
				int tk = 0;
				for (int k = 0; k < pp2.num_projs; k++) {
					l2 = pp2.funcs[k].l;
					for (int m2 = -l2; m2 <= l2; m2++) {
						if (1) {//R > 0.001) {
							overlaps[i][tj*pp2.total_projs+tk] =
								offsite_wave_overlap(dcoords + 3*i, pp1.wave_grid,
								pp1.funcs[j].diffwave,
								pp1.funcs[j].diffwave_spline, pp1.wave_gridsize,
								pp2.wave_grid, pp2.funcs[k].diffwave,
								pp2.funcs[k].diffwave_spline, pp2.wave_gridsize,
								wf_R->lattice, l1, m1, l2, m2);
						} else if (l1 == l2 && m1 == m2) {
							overlaps[i][tj*pp2.total_projs+tk] = pp2.diff_overlap_matrix[j*pp2.num_projs+k];
						}
						//printf("%d %d %d %d\n", N_RS_R[i], N_RS_S[i], tj, tk);
						tk++;
					}
				}
				tj++;
			}
		}
	}
	printf("%d %d %d %d stff\n", num_N_R, num_N_S, num_N_RS, 0);
	//printf("%p %p %p %p\n", overlaps[0], overlaps[1], overlaps[2], overlaps[3]);
	wf_S->overlaps = overlaps;
	wf_S->dcoords = dcoords;
	wf_S->num_aug_overlap_sites = num_N_RS;
	wf_R->num_aug_overlap_sites = num_N_RS;
	printf("finished overlap setup\n");
}

double* compensation_terms(int BAND_NUM, pswf_t* wf_S, pswf_t* wf_R, ppot_t* pps,
	int num_elems, int num_M, int num_N_R, int num_N_S, int num_N_RS,
	int* M_R, int* M_S, int* N_R, int* N_S, int* N_RS_R, int* N_RS_S,
	int* proj_labels, double* proj_coords, int* ref_labels, double* ref_coords,
	int* fft_grid) {

	printf("fftg %d\n", fft_grid[0]);
	printf("fftg %d\n", fft_grid[1]);
	printf("fftg %d\n", fft_grid[2]);
	//freopen("tst.out", "w", stdout);
	setbuf(stdout, NULL);

	printf("%d %d %d %d %d %d\n", BAND_NUM, num_elems, num_M, num_N_R, num_N_S, num_N_RS);
	printf("%d %lf %d %lf %d\n", proj_labels[0], proj_coords[0], ref_labels[0], ref_coords[0], fft_grid[1]);
	
	int NUM_KPTS = wf_S->nwk * wf_S->nspin;
	int NUM_BANDS = wf_S->nband;

	double* overlap = (double*) calloc(2 * NUM_KPTS * NUM_BANDS, sizeof(double));
	CHECK_ALLOCATION(overlap);

	double complex** N_RS_overlaps = wf_S->overlaps;
	double inv_sqrt_vol = pow(determinant(wf_R->lattice), -0.5);

	#pragma omp parallel for
	for (int w = 0; w < NUM_BANDS * NUM_KPTS; w++) {
		int ni = 0, nj = 0;
		int l1 = 0, l2 = 0;

		kpoint_t* kpt_R = wf_R->kpts[w%NUM_KPTS];
		kpoint_t* kpt_S = wf_S->kpts[w%NUM_KPTS];
		band_t* band_R = kpt_R->bands[w/NUM_KPTS];
		band_t* band_S = kpt_S->bands[BAND_NUM];

		double complex temp = 0 + 0 * I;
		int t = 0;
		for (int s = 0; s < num_M; s++) {
			ppot_t pp = pps[ref_labels[M_R[s]]];
			int s1 = M_R[s];
			int s2 = M_S[s];
			projection_t pron = band_R->projections[s1];
			projection_t ppron = band_S->projections[s2];

			for (int i = 0; i < pron.total_projs; i++) {
				for (int j = 0; j < ppron.total_projs; j++) {
					if (pron.ls[i] == ppron.ls[j]  && pron.ms[i] == ppron.ms[j]) {
						nj = pron.ns[j];
						ni = ppron.ns[i];
						temp += conj(pron.overlaps[j])
							//* (mymat[pron.num_projs*ni+nj])
							* (pp.aepw_overlap_matrix[pp.num_projs*ni+nj]
							- pp.pspw_overlap_matrix[pp.num_projs*ni+nj])
							* ppron.overlaps[i];
					}
				}
			}
		}
		overlap[2*w] = creal(temp);
		overlap[2*w+1]= cimag(temp);
		//printf("temp 1 %lf %lf\n", creal(temp), cimag(temp));

		temp = 0 + 0 * I;
		for (int s = 0; s < num_N_R; s++) {
			int site_num = N_R[s];
			projection_t pron = band_R->projections[site_num];
			projection_t ppron = band_S->wave_projections[s];
			ppot_t pp = pps[ref_labels[N_R[s]]];
			for (int i = 0; i < pp.total_projs; i++) {
				temp += ppron.overlaps[i] * conj(pron.overlaps[i]);
			}
		}
		overlap[2*w] += creal(temp);
		overlap[2*w+1]+= cimag(temp);
		//printf("temp 2 %lf %lf\n", creal(temp), cimag(temp));

		temp = 0 + 0 * I;
		for (int s = 0; s < num_N_S; s++) {
			ppot_t pp = pps[ref_labels[N_S[s]]];
			int site_num = N_S[s];
			projection_t pron = band_R->wave_projections[s];
			projection_t ppron = band_S->projections[site_num];
			for (int i = 0; i < pp.total_projs; i++) {
				temp += conj(pron.overlaps[i]) * ppron.overlaps[i];
			}
		}
		overlap[2*w] += creal(temp);
		overlap[2*w+1]+= cimag(temp);

		temp = 0 + 0 * I;
		for (int s = 0; s < num_N_RS; s++) {
			ppot_t pp = pps[ref_labels[N_RS_R[s]]];
			int site_num1 = N_RS_R[s];
			int site_num2 = N_RS_S[s];
			projection_t pron = band_R->projections[site_num1];
			projection_t ppron = band_S->projections[site_num2];
			for (int i = 0; i < pron.total_projs; i++) {
				for (int j = 0; j < ppron.total_projs; j++) {
					temp += conj(pron.overlaps[j])
						* (N_RS_overlaps[s][i*ppron.total_projs+j])
						* ppron.overlaps[i] * cexp(2*I*PI *
						dot(kpt_R->k, wf_S->dcoords + 3*s));
				}
			}
		}
		overlap[2*w] += creal(temp);
		overlap[2*w+1]+= cimag(temp);
		//printf("temp 4 %lf %lf\n", creal(temp), cimag(temp));
	}

	mkl_free_buffers();
	return overlap;
}