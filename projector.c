#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <omp.h>
#include <time.h>
#include <mkl.h>
#include <mkl_types.h>
#include "utils.h"
#include "projector.h"

#define PI 3.14159265359
#define c 0.262465831

void vc_pseudoprojection(pswf_t* wf_ref, pswf_t* wf_proj, int BAND_NUM, double* results) {

	clock_t start = clock();
	kpoint_t** kpts = wf_ref->kpts;
	kpoint_t** kptspro = wf_proj->kpts;
	int NUM_KPTS = wf_ref->nwk * wf_ref->nspin;
	int NUM_BANDS = wf_ref->nband;

	double* cband = (double*) calloc(NUM_KPTS, sizeof(double));
	double* vband = (double*) calloc(NUM_KPTS, sizeof(double));

	#pragma omp parallel for 
	for (int b = 0; b < NUM_BANDS; b++)
	{
		for (int kpt_num = 0; kpt_num < NUM_KPTS; kpt_num++)
		{
			float complex curr_overlap = 0;
			float complex* C1s = kptspro[kpt_num]->bands[0]->Cs;
			float complex* C2s = kpts[kpt_num]->bands[b]->Cs;
			int num_waves = kpts[kpt_num]->bands[b]->num_waves;
			for (int w = 0; w < num_waves; w++)
			{
				curr_overlap += C1s[w] * conj(C2s[w]);
			}
			#pragma omp critical
			{
				if (kpts[kpt_num]->bands[b]->occ > 0.5)
					vband[kpt_num] += creal((double) (curr_overlap * conj(curr_overlap)));
				else
					cband[kpt_num] += creal((double) (curr_overlap * conj(curr_overlap)));
			}
		}
	}

	double ctotal = 0.0;
	double vtotal = 0.0;
	for (int kpt_num = 0; kpt_num < NUM_KPTS; kpt_num++) {
		ctotal += cband[kpt_num] * kpts[kpt_num]->weight;
		vtotal += vband[kpt_num] * kpts[kpt_num]->weight;
	}

	printf("%lf\n", creal(kptspro[0]->bands[0]->energy));
	printf("c %lf\n", ctotal);
	printf("v %lf\n", vtotal);

	free(vband);
	free(cband);
	results[0] = vtotal;
	results[1] = ctotal;

	clock_t end = clock();
	printf("%lf seconds for band projection\n", (double)(end - start) / CLOCKS_PER_SEC);

}

float* pseudoprojection(pswf_t* wf_ref, pswf_t* wf_proj, int BAND_NUM) {

	kpoint_t** kpts = wf_ref->kpts;
	kpoint_t** kptspro = wf_proj->kpts;
	int NUM_KPTS = wf_ref->nwk * wf_ref->nspin;
	int NUM_BANDS = wf_ref->nband;

	float* projections = (float*) malloc(2*NUM_BANDS*NUM_KPTS*sizeof(float));

	#pragma omp parallel for 
	for (int b = 0; b < NUM_BANDS; b++)
	{
		for (int kpt_num = 0; kpt_num < NUM_KPTS; kpt_num++)
		{
			float complex curr_overlap = 0;
			float complex* C1s = kptspro[kpt_num]->bands[BAND_NUM]->Cs;
			float complex* C2s = kpts[kpt_num]->bands[b]->Cs;
			int num_waves = kpts[kpt_num]->bands[b]->num_waves;
			for (int w = 0; w < num_waves; w++)
			{
				curr_overlap += C1s[w] * conj(C2s[w]);
			}
			projections[2*(b*NUM_KPTS+kpt_num)] = creal(curr_overlap);
			projections[2*(b*NUM_KPTS+kpt_num)+1] = cimag(curr_overlap);
		}
	}

	return projections;
}

ppot_t* get_projector_list(int num_els, int* labels, int* ls, double* proj_grids, double* wave_grids,
	double* projectors, double* aewaves, double* pswaves) {
	ppot_t* pps = (ppot_t*) malloc(num_els * sizeof(ppot_t));
	int wt = 0;
	int pt = 0;
	int wgt = 0;
	int pgt = 0;
	for (int i = 0; i < num_els; i++) {
		pps[i].num_projs = labels[4*i+1];
		pps[i].proj_gridsize = labels[4*i+2];
		pps[i].wave_gridsize = labels[4*i+3];
		pps[i].wave_grid = (double*) malloc(pps[i].wave_gridsize*sizeof(double));
		for (int j = 0; j < pps[i].wave_gridsize; j++) {
			pps[i].wave_grid[j] = wave_grids[wgt];
			wgt++;
		}
		pps[i].proj_grid = (double*) malloc(pps[i].proj_gridsize*sizeof(double));
		for (int j = 0; j < pps[i].proj_gridsize; j++) {
			pps[i].proj_grid[j] = proj_grids[pgt];
			pgt++;
		}
		funcset_t* funcs = (funcset_t*) malloc(pps[i].num_projs*sizeof(funcset_t));
		for (int k = 0; k < pps[i].num_projs; k++) {
			funcs[k].proj = (double*) malloc(sizeof(double)*pps[i].proj_gridsize);
			funcs[k].aewave = (double*) malloc(sizeof(double)*pps[i].wave_gridsize);
			funcs[k].pswave = (double*) malloc(sizeof(double)*pps[i].wave_gridsize);
			for (int j = 0; j < pps[i].wave_gridsize; j++) {
				funcs[k].aewave[j] = aewaves[wt];
				funcs[k].pswave[j] = pswaves[wt];
				wt++;
			}
			for (int j = 0; j < pps[i].proj_gridsize; j++) {
				funcs[k].proj[j] = projectors[pt];
				pt++;
			}
		}
	}
	return pps;
}

double complex* onto_projector(int* labels, double* coords, int* G_bounds, double* lattice,
	double* kpt, int* Gs, float complex* Cs, int num_waves, int num_M, int* M, ppot_t* pps, int* fftg) {
	int dg1 = G_bounds[1] - G_bounds[0];
	int dg2 = G_bounds[3] - G_bounds[2];
	int dg3 = G_bounds[5] - G_bounds[4];
	int dg[3];
	dg[0]=dg1;
	dg[1]=dg2;
	dg[2]=dg3;
	MKL_LONG status = 0;
	DFTI_DESCRIPTOR_HANDLE handle = 0;
	MKL_Complex16* x = (MKL_Complex16*) mkl_calloc(fftg[0]*fftg[1]*fftg[2], sizeof(MKL_Complex16), 64);

	for (int w = 0; w < num_waves; w++) {
		int g1 = Gs[3*w]-G_bounds[0], g2 = Gs[3*w+1]-G_bounds[2], g3 = Gs[3*w+2]-G_bounds[4];
		x[g1*fftg[1]*fftg[2] + g2*fftg[2] + g3].real = creal(Cs[w]);
		x[g1*fftg[1]*fftg[2] + g2*fftg[2] + g3].imag = cimag(Cs[w]);
	}

	DftiCreateDescriptor(&handle, DFTI_DOUBLE, DFTI_COMPLEX, 3, dg);
	DftiCommitDescriptor(handle);
	DftiComputeBackward(handle, x);

	double kmins[3] = {G_bounds[0] + kpt[0], G_bounds[2] + kpt[1], G_bounds[4] + kpt[2]};

	int t_projs = 0;
	for (int i = 0; i < num_M; i++) {
		t_projs += pps[labels[M[i]]].num_projs * (2 * pps[labels[M[i]]].l + 1);
	}

	double complex* overlap = (double complex*) calloc(t_projs, sizeof(double complex));

	for (int i = 0; i < fftg[0]; i++) {
		for (int j = 0; j < fftg[1]; j++) {
			for (int k = 0; k  < fftg[2]; k++) {
				double frac[3] = {(double)i/fftg[0], (double)j/fftg[1], (double)k/fftg[2]};
				x[i*fftg[1]*fftg[2] + j*fftg[2] + k] *= cexp(I*2*PI*(dot(kmins, frac)));
				int t = 0;
				for (int q = 0; q < num_M; q++) {
					int p = M[q];
					ppot_t pp = pps[labels[p]];
					if (dist_from_frac(coords+3*p, frac, lattice) < pp.rmax) {
						for (int n = 0; n < pp.num_projs; n++) {
							for (int m = -pp.funcs[n].l; m <= pp.funcs[n].l; m++) {
								overlap[t] += proj_value(pp.funcs[n], m, pp.rmax, coords[3*p], frac, lattice)
											* x[i*fftg[1]*fftg[2] + j*fftg[2] + k]
											/ (fftg[0]*fftg[1]*fftg[2]);
								t++;
							}
						}
						break;
					}
				}
			}
		}
	}

	return overlap;
}

void make_pwave_overlap_matrices(ppot_t pp) {
	int size = pp.num_projs * pp.num_projs;

	double* psov = (double*) calloc(size, sizeof(double));
	double* aeov = (double*) calloc(size, sizeof(double));
	double* diov = (double*) calloc(size, sizeof(double));

	for (int i = 0; i < pp.num_projs; i++) {
		for (int j = i; j < pp.num_projs; j++) {
			if (pp.funcs[i].l == pp.funcs[j].l) {
				double* ps1 = pp.funcs[i].pswave;
				double* ps2 = pp.funcs[j].pswave;
				double* ae1 = pp.funcs[i].aewave;
				double* ae2 = pp.funcs[j].aewave;
				double dr = pp.wave_grid[0];
				double r = pp.wave_grid[0];
				for (int k = 0; k < pp.wave_gridsize; k++) {
					psov[4*i+j] += r * r * ps1[k] * ps2[k] * dr;
					aeov[4*i+j] += r * r * ae1[k] * ae2[k] * dr;
					diov[4*i+j] += r * r * (ae1[k]-ps1[k]) * (ae2[k]-ps2[k]) * dr;
					dr = pp.wave_grid[k+1] - pp.wave_grid[k];
				}
			}
		}
	}
	for (int i = 1; i < pp.num_projs; i++) {
		for (int j = 0; j < i; j++) {
			psov[4*i+j] = conj(psov[4*j+1]);
			aeov[4*i+j] = conj(aeov[4*j+1]);
			diov[4*i+j] = conj(diov[4*j+1]);
		}
	}

	pp.pspw_overlap_matrix = psov;
	pp.aepw_overlap_matrix = aeov;
	pp.diff_overlap_matrix = diov;
}

double complex* compensation_terms(int BAND_NUM, pswf_t* wf_proj, pswf_t* wf_ref, ppot_t* pps,
	int num_elems, int num_M, int num_N_R, int num_N_S, int num_N_RS,
	int* M, int* N_R, int* N_S, int* N_RS,
	int* proj_labels, double* proj_coords, int* ref_labels, double* ref_coords,
	int* fft_grid) {

	int NUM_KPTS = wf_proj->nwk * wf_proj->nspin;
	int NUM_BANDS = wf_proj->nband;
	#pragma omp parallel for 
	for (int p = 0; p < num_elems; p++) {
		make_pwave_overlap_matrices(pps[p]);
	}

	double complex* overlap = (double complex*) calloc(NUM_BANDS * NUM_KPTS, sizeof(double complex));

	double complex** lst_proj_projs = (double complex**) malloc(NUM_KPTS*sizeof(double complex*));
	#pragma omp parallel for 
	for (int kpt_num = 0; kpt_num < NUM_KPTS; kpt_num++) {
		lst_proj_projs[kpt_num] = onto_projector(proj_labels, proj_coords,
			wf_proj->G_bounds, wf_proj->lattice, wf_proj->kpts[kpt_num]->k,
			wf_proj->kpts[kpt_num]->Gs, wf_proj->kpts[kpt_num]->bands[BAND_NUM]->Cs,
			wf_proj->kpts[kpt_num]->bands[BAND_NUM]->num_waves, num_M, M, pps, fft_grid);
	}

	#pragma omp parallel for
	for (int w = 0; w < NUM_BANDS * NUM_KPTS; w++) {

		double complex* proj_projs = lst_proj_projs[w%NUM_KPTS];
		double complex* ref_projs = onto_projector(ref_labels, ref_coords,
			wf_ref->G_bounds, wf_ref->lattice, wf_ref->kpts[w%NUM_KPTS]->k,
			wf_ref->kpts[w%NUM_KPTS]->Gs, wf_ref->kpts[w%NUM_KPTS]->bands[w/NUM_KPTS]->Cs,
			wf_ref->kpts[w%NUM_KPTS]->bands[w/NUM_KPTS]->num_waves, num_M, M, pps, fft_grid);

		int t = 0;
		for (int s = 0; s < num_M; s++) {
			ppot_t pp = pps[ref_labels[M[s]]];
			int ti = 0;
			for (int i = 0; i < pp.num_projs; i++) {
				for (int temp1 = 0; temp1 < 2*pp.funcs[i].l+1; temp1++) {
					int tj = 0;
					for (int j = 0; j < pp.num_projs; j++) {
						for (int temp2 = 0; temp2 < 2*pp.funcs[j].l+1; temp2++) {
							overlap[w] += conj(ref_projs[t+tj])
								* (pp.aepw_overlap_matrix[4*i+j] - pp.pspw_overlap_matrix[4*i+j])
								* proj_projs[t+i];
							tj++;
						}
					}
					ti++;
				}
			}
			t += pp.num_projs;
		}

		free(ref_projs);
	}

	for (int kpt_num = 0; kpt_num < NUM_KPTS; kpt_num++) {
		free(lst_proj_projs[kpt_num]);
	}
	free(lst_proj_projs);

	return overlap;
}
