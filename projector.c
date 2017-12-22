#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <omp.h>
#include <time.h>
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
		funcset_t* funcs = (functset_t*) malloc(pps[i].num_projs*sizeof(funcset_t));
		for (int k = 0; k < pps[i].num_projs; k++) {
			funcs[k].proj = (double*) malloc(sizeof(double)*pps[i].proj_gridsize);
			funcs[k].aewave = (double*) malloc(sizeof(double)*pps[i].wave_gridsize);
			funcs[k].pswave = (double*) malloc(sizeof(double)*pps[i].wave_gridsize);
			for (int j = 0; j < funcs.wave_gridsize; j++) {
				funcs[k].aewave[j] = aewaves[wt];
				funcs[k].pswave[j] = pswaves[wt];
				wt++;
			}
			for (int j = 0; j < funcs.proj_gridsize; j++) {
				funcs[k].proj[j] = projectors[pt];
				pt++;
			}
		}
	}
	return pps;
}

double* onto_projector(int* labels, double* coords, int* G_bounds, int)

double* compensation_terms(pswf* wf_proj, pswf* wf_ref, ppot_t* pps,
	int num_M, int num_N_R, int num_N_S, int num_N_RS,
	int* M, int* N_R, int* N_S, int* N_RS,
	int* proj_labels, double* proj_coords, int* ref_labels, double* ref_coords) {

	NUM_KPTS = wf_proj->nwk * wf_proj->nspin;
	NUM_BANDS = wf_proj->nband;
	#pragma omp parallel for
	for (int w = 0; w < NUM_BANDS * NUM_KPTS; w++) {
		onto_projector(proj_labels, proj_coords,
			wf_proj->G_bounds, wf_proj->lattice, wf_proj->kpts[w%NUM_KPTS]->Gs,
			wf_proj->kpts[w%NUM_KPTS]->bands[w/NUM_KPTS]->Cs, num_M, M, pps);
	}
}
/*
double* read_and_project(int BAND_NUM, double* kpt_weights, char* bulkfile, char* defectfile) {
	printf("%lf\n", kpt_weights[0]);
	printf("%lf\n", kpt_weights[1]);
	printf("%lf\n", kpt_weights[5]);
	int* G_bounds = (int*) malloc(6*sizeof(double));
	double* results = (double*) malloc(2*sizeof(double));
	int NUM_SPINS, NUM_KPTS, NUM_BANDS;
	kpoint_t** kptspro = read_one_band(G_bounds, kpt_weights, &NUM_SPINS, &NUM_KPTS, &NUM_BANDS, BAND_NUM, defectfile);
	kpoint_t** kptsref = read_wavefunctions(G_bounds, kpt_weights, &NUM_SPINS, &NUM_KPTS, &NUM_BANDS, bulkfile);
	get_band_projection(BAND_NUM, NUM_KPTS, NUM_BANDS, kptsref, kptspro, G_bounds, results);
	for (int kpt_num = 0; kpt_num < NUM_KPTS; kpt_num++) {
		free_kpoint(kptsref[kpt_num]);
		free_kpoint(kptspro[kpt_num]);
	}
	free(kptsref);
	free(kptspro);
	free(G_bounds);
	return results;
}*/


