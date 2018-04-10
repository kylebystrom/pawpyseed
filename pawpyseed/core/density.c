#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <omp.h>
#include <time.h>
#include <mkl.h>
#include <mkl_types.h>
#include "utils.h"
#include "density.h"
#include "linalg.h"

#define PI 3.14159265359

/*double* ae_chg_density(pswf_t* wf, int num_sites, double* coords, double* labels, ppot_t* pps, int* fftg) {
	
	real_proj_site_t* sites = projector_values(num_sites, labels, coords,
		wf->lattice, wf->reclattice, pps, fftg);
	int NUM_BANDS = wf->nband;
	int NUM_KPTS = wf->nwk;

	int spin_mult = 2 / wf->nspin;
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
} */

double* ae_chg_density(pswf_t* wf, ppot_t* pps, int* fftg, int* labels, double* coords) {

	int gridsize = fftg[0] * fftg[1] * fftg[2];
	double* P = mkl_calloc(gridsize, sizeof(double), 64);
	int spin_mult = 2 / wf->nspin;
	for (int k = 0; k < wf->nwk * wf->nspin; k++) {
		printf("KLOOP %d\n", k);
		for (int b = 0; b < wf->nband; b++) {
			if (wf->kpts[k]->bands[b]->occ > 0.00000001) {
				double complex* x = realspace_state(b, k, wf, pps, fftg, labels, coords);
				for (int i = 0; i < gridsize; i++) {
					P[i] += creal(x[i] * conj(x[i])) * wf->kpts[k]->weight * wf->kpts[k]->bands[b]->occ * spin_mult;
				}
				mkl_free(x);
			}
		}
	}
	mkl_free_buffers();

	return P;
}

double complex* realspace_state(int BAND_NUM, int KPOINT_NUM, pswf_t* wf, ppot_t* pps, int* fftg,
		int* labels, double* coords) {

	double complex* x = mkl_calloc(fftg[0]*fftg[1]*fftg[2], sizeof(MKL_Complex16), 64);
	printf("START FT\n");
	fft3d(x, wf->G_bounds, wf->lattice, wf->kpts[KPOINT_NUM]->k,
		wf->kpts[KPOINT_NUM]->Gs, wf->kpts[KPOINT_NUM]->bands[BAND_NUM]->Cs,
		wf->kpts[KPOINT_NUM]->bands[BAND_NUM]->num_waves, fftg);
	printf("FINISH FT\n");
	double* lattice = wf->lattice;
	double vol = determinant(lattice);

	int num_sites = wf->num_sites;
	//#pragma omp parallel for
	for (int p = 0; p < num_sites; p++) {
		projection_t pros = wf->kpts[KPOINT_NUM]->bands[BAND_NUM]->projections[p];
		printf("READ PROJECTIONS\n");
		ppot_t pp = pps[labels[p]];
		double rmax = pp.wave_grid[pp.wave_gridsize-1];
		double res[3] = {0,0,0};
		double frac[3] = {0,0,0};
		double testcoord[3] = {0,0,0};
		vcross(res, lattice+3, lattice+6);
		int grid1 = (int) (mag(res) * rmax / vol * fftg[0]) + 1;
		vcross(res, lattice+0, lattice+6);
		int grid2 = (int) (mag(res) * rmax / vol * fftg[1]) + 1;
		vcross(res, lattice+0, lattice+3);
		int grid3 = (int) (mag(res) * rmax / vol * fftg[2]) + 1;
		int center1 = (int) round(coords[3*p+0] * fftg[0]);
		int center2 = (int) round(coords[3*p+1] * fftg[1]);
		int center3 = (int) round(coords[3*p+2] * fftg[2]);
		int ii=0, jj=0, kk=0;
		for (int i = -grid1 + center1; i <= grid1 + center1; i++) {
			for (int j = -grid2 + center2; j <= grid2 + center2; j++) {
				for (int k = -grid3 + center3; k <= grid3 + center3; k++) {
					testcoord[0] = (double) i / fftg[0] - coords[3*p+0];
					testcoord[1] = (double) j / fftg[1] - coords[3*p+1];
					testcoord[2] = (double) k / fftg[2] - coords[3*p+2];
					frac_to_cartesian(testcoord, lattice);
					if (mag(testcoord) < rmax) {
						ii = (i%fftg[0] + fftg[0]) % fftg[0];
						jj = (j%fftg[1] + fftg[1]) % fftg[1];
						kk = (k%fftg[2] + fftg[2]) % fftg[2];
						frac[0] = (double) ii / fftg[0];
						frac[1] = (double) jj / fftg[1];
						frac[2] = (double) kk / fftg[2];
						projection_t pros = wf->kpts[KPOINT_NUM]->bands[BAND_NUM]->projections[p];
						for (int n = 0; n < pros.total_projs; n++) {
							x[ii*fftg[1]*fftg[2] + jj*fftg[2] + kk] += wave_value(pps[labels[p]].funcs[pros.ns[n]],
								pp.wave_gridsize, pps[labels[p]].wave_grid,
								pros.ms[n], coords+3*p, frac, lattice)
								* pros.overlaps[n];
						}
					}
				}
			}
		}
	}

	return x;
}

double* realspace_state_ri(int BAND_NUM, int KPOINT_NUM, pswf_t* wf, ppot_t* pps, int* fftg,
		int* labels, double* coords) {

	double complex* x = realspace_state(BAND_NUM, KPOINT_NUM, wf, pps, fftg, labels, coords);

	int gridsize = fftg[0]*fftg[1]*fftg[2];

	double* rpip = (double*) malloc(2*gridsize);

	for (int i = 0; i < gridsize; i++) {
		rpip[i] = creal(x[i]);
		rpip[i+gridsize] = cimag(x[i]);
	}

	return rpip;
}

void write_volumetric(char* filename, double* x, int* fftg, double scale) {

	FILE* fp = fopen(filename, "w");
	int t = 1;
	for (int k = 0; k < fftg[2]; k++) {
		for (int j = 0; j < fftg[1]; j++) {
			for (int i = 0; i < fftg[0]; i++) {
				fprintf(fp, "%E   ", x[i*fftg[1]*fftg[2] + j*fftg[2] + k] * scale);
				if (t % 5 == 0) fprintf(fp, "\n");
				t++;
			}
		}
	}
	fclose(fp);
}

double* write_realspace_state_ri_return(char* filename1, char* filename2, int BAND_NUM, int KPOINT_NUM,
	pswf_t* wf, ppot_t* pps, int* fftg,
	int* labels, double* coords) {

	double* x = realspace_state_ri(BAND_NUM, KPOINT_NUM, wf, pps, fftg, labels, coords);
	write_volumetric(filename1, x, fftg, 1);
	write_volumetric(filename2, x+fftg[0]*fftg[1]*fftg[2], fftg, 1);

	return x;
}

double* write_density_return(char* filename, pswf_t* wf, ppot_t* pps,
	int* fftg, int* labels, double* coords) {

	double* x = ae_chg_density(wf, pps, fftg, labels, coords);
	double scale = determinant(wf->lattice);
	write_volumetric(filename, x, fftg, scale);

	return x;
}

void write_realspace_state_ri_noreturn(char* filename1, char* filename2, int BAND_NUM, int KPOINT_NUM,
	pswf_t* wf, ppot_t* pps, int* fftg,
	int* labels, double* coords) {

	mkl_free(write_realspace_state_ri_return(filename1, filename2,
		BAND_NUM, KPOINT_NUM, wf, pps, fftg,
		labels, coords));
}

void write_density_noreturn(char* filename, pswf_t* wf, ppot_t* pps,
	int* fftg, int* labels, double* coords) {
	setbuf(stdout, NULL);

	double* x = write_density_return(filename, wf, pps, fftg, labels, coords);
	mkl_free(x);
}
