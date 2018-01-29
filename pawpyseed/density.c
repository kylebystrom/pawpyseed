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
#include "fft.h"
#include "quadrature.h"

#define PI 3.14159265359

double* ae_chg_density(pswf_t* wf, int num_sites, double* coords, double* labels, ppot_t* pps, int* fftg) {
	
	real_proj_site_t* sites = projector_values(num_sites, labels, coords,
		wf->lattice, wf->reclattice, pps, fftg);
	int NUM_BANDS = wf->nband;
	int NUM_KPTS = wf->nwk;

	int t_projs = 0;
	for (int i = 0; i < num_sites; i++) {
		t_projs += pps[labels[i]].total_projs;
	}

	for (int k = 0; k < NUM_KPTS; k++) {
		for (int b = 0; b < NUM_BANDS; b++) {

			MKL_Complex16* x = (MKL_Complex16*) mkl_calloc(fftg[0]*fftg[1]*fftg[2], sizeof(MKL_Complex16), 64);
			fft3d(x, wf->G_bounds, wf->lattice, wf->kpts[k]->kpt, wf->kpts[i]->Gs, wf->kpts[i]->bands[b]->Cs,
				wf->kpts[i]->bands[b]->num_waves, fftg);

			double complex* projs = onto_projector_helper(sites, labels, wf->lattice, pps, fftg);
			double sq_projs = malloc(t_projs * sizeof(double));
			for (int pnum = 0; pnum < t_projs; pnum++) {
				sq_projs[pnum] = pow(creal(projs[pnum]), 2) + pow(cimag(projs[pnum]), 2);
			}

			double complex temp = 0;;
			int t = 0;
			for (int s = 0; s < num_sites; s++) {
				ppot_t pp = pps[labels[s]];
				int ti = 0;
				for (int i = 0; i < pp.num_projs; i++) {
					l1 = pp.funcs[i].l;
					for (int m1 = -l1; m1 <= l1; m1++) {
						int tj = 0;
						for (int j = 0; j < pp.num_projs; j++) {
							l2 = pp.funcs[j].l;
							for (int m2 = -l2; m2 <= l2; m2++) {
								if (l1 == l2 && m1 == m2) {
									ae1 = wave_value(pp.funcs[i].aewave, pp.wave_grid, m1,
										l1, ion_pos, pos, lattice);
									temp += conj(projs[t+ti])
										* (pp.aepw_overlap_matrix[pp.num_projs*i+j]
										- pp.pspw_overlap_matrix[pp.num_projs*i+j])
										* projs[t+tj];
								}
								tj++;
							}
						}
						ti++;
					}
				}
				t += ti;
			}
			free(ref_projs);
			overlap[2*w] = creal(temp);
			overlap[2*w+1]= cimag(temp);

			free(projs);
			mkl_free(x);
		}
	}

	free_real_proj_site_list(sites);
	mkl_free_buffers();
}