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

double* project_realspace_state(int BAND_NUM, int numtoproj, pswf_t* wf, pswf_t* wf_R, ppot_t* pps, int* fftg,
	int* labels, double* coords, int* labels_R, double* coords_R) {

	int nband = wf->nband;
	int nwk = wf->nwk;
	int nspin = wf->nspin;
	int gridsize = fftg[0]*fftg[1]*fftg[2];
	double* projs = (double*) malloc(2*nband*nwk*nspin*sizeof(double));
	double vol = determinant(wf->lattice);

	double complex overlap = 0;
	for (int k = 0; k < nwk * nspin; k++) {
		double complex* state = realspace_state(BAND_NUM, k, wf, pps, fftg, labels, coords);
		for (int b = 0; b < nband; b++) {
			double complex* state_R = realspace_state(b, k, wf_R, pps, fftg, labels_R, coords_R);
			cblas_zdotc_sub(gridsize, state_R, 1, state, 1, &overlap);
			overlap *= vol / gridsize;
			projs[b*nwk*nspin + k] = creal(overlap);
			projs[nband*nwk*nspin + b*nwk*nspin + k] = cimag(overlap);
			mkl_free(state_R);
		}
		mkl_free(state);
	}

	return projs;
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
	#pragma omp parallel for
	for (int i = 0; i < fftg[0]; i++) {
		double frac[3] = {0,0,0};
		double kdotr = 0;
		for (int j = 0; j < fftg[1]; j++) {
			for (int k = 0; k < fftg[2]; k++) {
				frac[0] = (double) i / fftg[0];
				frac[1] = (double) j / fftg[1];
				frac[2] = (double) k / fftg[2];
				kdotr = dot(wf->k, frac);
				x[i*fftg[1]*fftg[2] + j*fftg[1] + k] *= cexp(2*PI*I*kdotr);
			}
		}
	}

	int num_sites = wf->num_sites;
	#pragma omp parallel for
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
							x[ii*fftg[1]*fftg[2] + jj*fftg[2] + kk] += wave_value(pp.funcs[pros.ns[n]],
								pp.wave_gridsize, pp.wave_grid,
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

	double* rpip = (double*) malloc(2*gridsize * sizeof(double));

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

	double* x = write_realspace_state_ri_return(filename1, filename2,
		BAND_NUM, KPOINT_NUM, wf, pps, fftg,
		labels, coords);
	mkl_free(x);
}

void write_density_noreturn(char* filename, pswf_t* wf, ppot_t* pps,
	int* fftg, int* labels, double* coords) {
	setbuf(stdout, NULL);

	double* x = write_density_return(filename, wf, pps, fftg, labels, coords);
	mkl_free(x);
}