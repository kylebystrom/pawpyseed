#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include "util.h"
#include "projector.h"
#include "reader.h"

#define NUM_KPTS 150

int main(int argc, char **argv) {
	if (strcmp(argv[0], "read") == 0) {
		double kpt_wts[NUM_KPTS] = { 0.0 };
		pswf_t* wf = read_wavefunctions("WAVECAR", kpt_wts);
		free_pswf(wf);
	} else if (strcmp(argv[0], "pseudo") == 0) {
		double kpt_wts[NUM_KPTS] = { 0.0 };
		pswf_t* wf1 = read_wavefunctions("WAVECAR", kpt_wts);
		pswf_t* wf2 = read_wavefunctions("WAVECAR", kpt_wts);
		double* proj = pseudoprojection(wf2, wf1, 0);
		free(proj);
		free_pswf(wf);
	} else if (strcmp(argv[0], "list") == 0) {
		int labels[4] = {0, 1, 5, 5};
		int ls[1] = {0};
		double rm[1] = {1.0};
		double pg[10] = {0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9};
		double wg[10] = {0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9};
		double ae[10] = {0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9};
		double ps[10] = {0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9};
		int ls[]
		ppot_t* lst = get_projector_list(1, labels, ls, pg, wg, p, aw, ps, rm);
		free_ppot_list(lst);
	} else {
		printf("ERROR: unrecognized test");
	}
}