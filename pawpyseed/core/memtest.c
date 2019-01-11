#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <string.h>
#include "utils.h"
#include "pseudoprojector.h"
#include "projector.h"
#include "reader.h"
#include "tests.h"

#define NUM_KPTS 2

void test_read_wavefunctions() {
	double kpt_wts[NUM_KPTS] = { 0.0 };
	pswf_t* wf = read_wavefunctions("WAVECAR", kpt_wts);
	free_pswf(wf);
}

void test_pseudoprojection() {
	double kpt_wts[NUM_KPTS] = { 0.0 };
	pswf_t* wf1 = read_wavefunctions("WAVECAR", kpt_wts);
	pswf_t* wf2 = read_wavefunctions("WAVECAR", kpt_wts);
	double* proj = pseudoprojection(wf2, wf1, 0);
	free(proj);
	free_pswf(wf1);
	free_pswf(wf2);
}

void test_get_projector_list() {
	int labels[4] = {0, 1, 10, 10};
	int ls[1] = {0};
	double rm[1] = {1.0};
	double pg[10] = {0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9};
	double wg[10] = {0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9};
	double ae[10] = {0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9};
	double ps[10] = {0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9};
	double p[10] = {0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9};
	ppot_t* lst = get_projector_list(1, labels, ls, pg, wg, p, ae, ps, rm, 7000);
	
	make_pwave_overlap_matrices(lst);
	free_ppot_list(lst, 1);
}
/*
void test_fft() {
	fft_check();
}
*/
void test_compensation_terms() {
	FILE* fp = fopen("potholder.txt", "r");
	int num_els = 0;
	fscanf(fp, "%d ", &num_els);
	int* labels = (int*) malloc(num_els*4 * sizeof(int));
	for (int i = 0; i < num_els; i++) {
		fscanf(fp, "%d %d %d %d ", labels+4*i, labels+4*i+1, labels+4*i+2, labels+4*i+3);
	}
	printf("%d %d %d %d\n", labels[0], labels[1], labels[2], labels[3]);

	printf("ls\n");
	int lenls = 0;
	fscanf(fp, "%d ", &lenls);
	int* ls = (int*) malloc(lenls * sizeof(double));
	for (int i = 0; i < lenls; i++) {
		fscanf(fp, "%d ", ls+i);
	}

	int lenwgrids = 0;
	fscanf(fp, "%d ", &lenwgrids);
	printf("wgrids %d\n", lenwgrids);
	double* wgrids = (double*) malloc(lenwgrids * sizeof(double));
	for (int i = 0; i < lenwgrids; i++) {
		fscanf(fp, "%lf ", wgrids + i);
	}

	int lenprojectors = 0;
	fscanf(fp, "%d ", &lenprojectors);
	printf("projs %d\n", lenprojectors);
	double* projectors = (double*) malloc(lenprojectors * sizeof(double));
	for (int i = 0; i < lenprojectors; i++) {
		fscanf(fp, "%lf ", projectors + i);
	}

	int lenaewaves = 0;
	fscanf(fp, "%d ", &lenaewaves);
	printf("aewaves %d\n", lenaewaves);
	double* aewaves = (double*) malloc(lenaewaves * sizeof(double));
	for (int i = 0; i < lenaewaves; i++) {
		fscanf(fp, "%lf ", aewaves + i);
	}

	int lenpswaves = 0;
	fscanf(fp, "%d ", &lenpswaves);
	printf("pswaves %d\n", lenpswaves);
	double* pswaves = (double*) malloc(lenpswaves * sizeof(double));
	for (int i = 0; i < lenpswaves; i++) {
		fscanf(fp, "%lf ", pswaves + i);
	}

	double rmaxs[1] = {0};
	fscanf(fp, "%lf ", rmaxs);

	int lensites = 0;
	fscanf(fp, "%d ", &lensites);
	int* selfnums = (int*) malloc(lensites * sizeof(int));
	for (int i = 0; i < lensites; i++) {
		fscanf(fp, "%d ", selfnums + i);
	}

	double* selfcoords = (double*) malloc(lensites * 3 * sizeof(double));
	for (int i = 0; i < lensites * 3; i++) {
		fscanf(fp, "%lf ", selfcoords + i);
	}

	fclose(fp);

	char* rm[1] = {"2.1042"};
	double rmm = 0;
	sscanf(rm[0], "%lf", &rmm);
	printf("pps\n");
	ppot_t* pps = get_projector_list(num_els, labels, ls, wgrids,
		projectors, aewaves, pswaves, &rmm, 7000);

	free(labels);
	free(wgrids);
	free(projectors);
	free(aewaves);
	free(pswaves);

	int M[4] = {0,1,2,3};
	int N_S[1] = {0};
	int fftg[3] = {30,30,30};
	double kws[NUM_KPTS] = {0};

	pswf_t* wf_ref = read_wavefunctions("WAVECAR", kws);
	//pswf_t* wf_proj = read_wavefunctions("WAVECAR", kws);
	double ops[18] = {1,0,0,
					0,1,0,
					0,0,1,
					1,0,0,
					0,1,0,
					0,0,1};
	int maps[2] = {0,0};
	double drs[6] = {0,0,0,0,0,0};
	int trs[2] = {0,0};
	pswf_t* wf_proj = expand_symm_wf(wf_ref, 2, maps, ops, drs, kws, trs);

	printf("terms\n");
	setup_projections(wf_proj, pps, num_els, 4, fftg, selfnums, selfcoords);
	setup_projections(wf_ref, pps, num_els, 4, fftg, selfnums, selfcoords);
	overlap_setup_real(wf_ref, wf_proj, pps, selfnums, selfnums, selfcoords, selfcoords, M, M, M, M, 4, 4, 4);
	double* terms = compensation_terms(0, wf_proj, wf_ref, pps,
		4, 0, 0, 0, M, M, N_S, N_S, N_S, N_S, selfnums, selfcoords, selfnums, selfcoords, fftg);
	double* terms2 = compensation_terms(0, wf_proj, wf_ref, pps,
		0, 4, 4, 4, N_S, N_S, M, M, M, M, selfnums, selfcoords, selfnums, selfcoords, fftg);

	free(ls);
	free_pswf(wf_ref);
	free_pswf(wf_proj);
	free_ppot_list(pps, 1);
	free(terms);
	free(terms2);
	free(selfnums);
	free(selfcoords);

}

int main(int argc, char **argv) {
	freopen("mtst.out", "w", stdout);
        setbuf(stdout, NULL);

	//test_read_wavefunctions();
	//test_pseudoprojection();
	//test_get_projector_list();
	test_compensation_terms();
}
