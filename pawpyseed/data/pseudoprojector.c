#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <omp.h>
#include <time.h>
#include "pseudoprojector.h"
#include "utils.h"

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

double* pseudoprojection(pswf_t* wf_ref, pswf_t* wf_proj, int BAND_NUM) {

	kpoint_t** kpts = wf_ref->kpts;
	kpoint_t** kptspro = wf_proj->kpts;
	int NUM_KPTS = wf_ref->nwk * wf_ref->nspin;
	int NUM_BANDS = wf_ref->nband;

	double* projections = (double*) malloc(2*NUM_BANDS*NUM_KPTS*sizeof(double));

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