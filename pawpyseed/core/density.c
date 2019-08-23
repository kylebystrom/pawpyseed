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

// THE FOLLOWING TWO FUNCTIONS ARE NOT YET IMPLEMENTED
/*
double* ncl_ae_state_density(int BAND_NUM, pswf_t* wf, int* fftg, int* labels, double* coords) {
	int gridsize = fftg[0] * fftg[1] * fftg[2];
    double* P = mkl_calloc(gridsize, sizeof(double), 64);
    int spin_mult = 2 / wf->nspin;
    double complex* x = realspace_state(b, k, wf, fftg, labels, coords);
    for (int i = 0; i < gridsize; i++) {
        P[i] += creal(x[i] * conj(x[i]));
    }
    mkl_free(x);
    return P;	
}
*/

void ae_state_density(double* P, int BAND_NUM, int KPOINT_NUM, pswf_t* wf,
	int* fftg, int* labels, double* coords) {

	int gridsize = fftg[0] * fftg[1] * fftg[2];
	//double* P = mkl_calloc(gridsize, sizeof(double), 64);
	double complex* x = mkl_malloc(gridsize * sizeof(double complex), 64);
	realspace_state(x, BAND_NUM, KPOINT_NUM,
		wf, fftg, labels, coords);
	for (int i = 0; i < gridsize; i++) {
		P[i] += creal(x[i] * conj(x[i]));
	}
	mkl_free(x);
}

/*
void ae_state_density(double* P, int BAND_NUM, int KPOINT_NUM, pswf_t* wf,
	int* fftg, int* labels, double* coords) {

	ppot_t* pps = wf->pps;
	int num_sites = wf->num_sites;

	double complex* x = (double complex*) mkl_malloc(fftg[0]*fftg[1]*fftg[2]*sizeof(double complex), 64);
	fft3d(x, wf->G_bounds, wf->lattice, wf->kpts[KPOINT_NUM]->k,
		wf->kpts[KPOINT_NUM]->Gs, wf->kpts[KPOINT_NUM]->bands[BAND_NUM]->Cs,
		wf->kpts[KPOINT_NUM]->bands[BAND_NUM]->num_waves, fftg);
	//printf("FINISH FT\n");
	double* lattice = wf->lattice;
	double vol = determinant(lattice);
	for (int i = 0; i < fftg[0]; i++) {
		double frac[3] = {0,0,0};
		double kdotr = 0;
		for (int j = 0; j < fftg[1]; j++) {
			for (int k = 0; k < fftg[2]; k++) {
				P[i*fftg[1]*fftg[2] + j*fftg[2] + k] = creal(x[i*fftg[1]*fftg[2] + j*fftg[2] + k]
													* conj(x[i*fftg[1]*fftg[2] + j*fftg[2] + k]));
			}
		}
	}

	for (int p = 0; p < num_sites; p++) {
		projection_t pros = wf->kpts[KPOINT_NUM]->bands[BAND_NUM]->projections[p];
		//printf("READ PROJECTIONS\n");
		ppot_t pp = pps[labels[p]];
		double rmax = pp.wave_grid[pp.wave_gridsize-1];
		double res[3] = {0,0,0};
		vcross(res, lattice+3, lattice+6);
		int grid1 = (int) (mag(res) * rmax / vol * fftg[0]) + 1;
		vcross(res, lattice+0, lattice+6);
		int grid2 = (int) (mag(res) * rmax / vol * fftg[1]) + 1;
		vcross(res, lattice+0, lattice+3);
		int grid3 = (int) (mag(res) * rmax / vol * fftg[2]) + 1;
		int center1 = (int) round(coords[3*p+0] * fftg[0]);
		int center2 = (int) round(coords[3*p+1] * fftg[1]);
		int center3 = (int) round(coords[3*p+2] * fftg[2]);
		//printf("FINISH SETUP %d\n%d %d %d\n%d %d %d\n",p, center1, center2, center3, grid1, grid2, grid3);
		for (int i = -grid1 + center1; i <= grid1 + center1; i++) {
			double frac[3] = {0,0,0};
			double testcoord[3] = {0,0,0};
			int ii=0, jj=0, kk=0;
			double phasecoord[3] = {0,0,0};
			double phase = 0;
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
						phasecoord[0] = coords[3*p+0] + ((ii-i) / fftg[0]);
						phasecoord[1] = coords[3*p+1] + ((jj-j) / fftg[1]);
						phasecoord[2] = coords[3*p+2] + ((kk-k) / fftg[2]);
						phase = dot(phasecoord, wf->kpts[KPOINT_NUM]->k);
						for (int n = 0; n < pros.total_projs; n++) {
							for (int m = 0; m < pros.total_projs; m++) {
								P[ii*fftg[1]*fftg[2] + jj*fftg[2] + kk] += creal(
									wave_value2(pp.wave_grid,
									pp.funcs[pros.ns[n]].aewave,
									pp.funcs[pros.ns[n]].aewave_spline,
									pp.wave_gridsize,
									pros.ls[n], pros.ms[n],
									testcoord)
									* conj(wave_value2(pp.wave_grid,
									pp.funcs[pros.ns[m]].aewave,
									pp.funcs[pros.ns[m]].aewave_spline,
									pp.wave_gridsize,
									pros.ls[m], pros.ms[m],
									testcoord))
									* pros.overlaps[n] * conj(pros.overlaps[m])
									);

								P[ii*fftg[1]*fftg[2] + jj*fftg[2] + kk] -= creal(
									wave_value2(pp.wave_grid,
									pp.funcs[pros.ns[n]].pswave,
									pp.funcs[pros.ns[n]].pswave_spline,
									pp.wave_gridsize,
									pros.ls[n], pros.ms[n],
									testcoord)
									* conj(wave_value2(pp.wave_grid,
									pp.funcs[pros.ns[m]].pswave,
									pp.funcs[pros.ns[m]].pswave_spline,
									pp.wave_gridsize,
									pros.ls[m], pros.ms[m],
									testcoord))
									* pros.overlaps[n] * conj(pros.overlaps[m])
									);
							}
								
							//	wave_value(pp.funcs[pros.ns[n]],
							//	pp.wave_gridsize, pp.wave_grid,
							//	pros.ms[n], coords+3*p, frac, lattice)
							//	* pros.overlaps[n] * cexp(2*PI*I*phase);
							//	* Ylm(thetaphi[0], thetaphi[1]);
						}
					}
				}
			}
		}
	}
	mkl_free(x);
}
*/

void ae_chg_density(double* P, pswf_t* wf, int* fftg, int* labels, double* coords) {

	int gridsize = fftg[0] * fftg[1] * fftg[2];
	double complex* x = mkl_malloc(gridsize * sizeof(double complex), 64);
	//double* P = mkl_calloc(gridsize, sizeof(double), 64);
	int spin_mult = 2 / wf->nspin;
	for (int k = 0; k < wf->nwk * wf->nspin; k++) {
		//printf("KLOOP %d\n", k);
		for (int b = 0; b < wf->nband; b++) {
			if (wf->kpts[k]->bands[b]->occ > 0) {
				realspace_state(x, b, k, wf, fftg, labels, coords);
				for (int i = 0; i < gridsize; i++) {
					P[i] += creal(x[i] * conj(x[i])) * wf->kpts[k]->weight
							* wf->kpts[k]->bands[b]->occ * spin_mult;
				}
			}
		}
	}
	mkl_free(x);
	mkl_free_buffers();
}

void ncl_ae_chg_density(double* P, pswf_t* wf, int* fftg, int* labels, double* coords) {

	int gridsize = fftg[0] * fftg[1] * fftg[2];
	double complex* x = mkl_malloc(2 * gridsize * sizeof(double complex), 64);
	//double* P = mkl_calloc(gridsize, sizeof(double), 64);
	int spin_mult = 1;
	for (int k = 0; k < wf->nwk * wf->nspin; k++) {
		//printf("KLOOP %d\n", k);
		for (int b = 0; b < wf->nband; b++) {
			if (wf->kpts[k]->bands[b]->occ > 0) {
				ncl_realspace_state(x, b, k, wf, fftg, labels, coords);
				for (int i = 0; i < gridsize; i++) {
					P[i] += (creal(x[i] * conj(x[i]))
								+ creal(x[i+gridsize] * conj(x[i+gridsize])))
							* wf->kpts[k]->weight
							* wf->kpts[k]->bands[b]->occ * spin_mult;
				}
			}
		}
	}
	mkl_free(x);
	mkl_free_buffers();
}

void project_realspace_state(double complex* projs, int BAND_NUM, pswf_t* wf, pswf_t* wf_R,
	int* fftg, int* labels, double* coords, int* labels_R, double* coords_R) {

	int nband = wf->nband;
	int nwk = wf->nwk;
	int nspin = wf->nspin;
	int gridsize = fftg[0]*fftg[1]*fftg[2];
	//double* projs = (double*) malloc(2*nband*nwk*nspin*sizeof(double));
	double vol = determinant(wf->lattice);
	double complex* state = mkl_malloc(gridsize * sizeof(double complex), 64);
	double complex* state_R = mkl_malloc(gridsize * sizeof(double complex), 64);

	for (int k = 0; k < nwk * nspin; k++) {
		double complex overlap = 0;
		realspace_state(state, BAND_NUM, k, wf, fftg, labels, coords);
		for (int b = 0; b < nband; b++) {
			realspace_state(state_R, b, k, wf_R, fftg, labels_R, coords_R);
			cblas_zdotc_sub(gridsize, state_R, 1, state, 1, &overlap);
			overlap *= vol / gridsize;
			projs[b*nwk*nspin + k] = overlap;
		}
	}
	mkl_free(state_R);
	mkl_free(state);
}

void realspace_state(double complex* x, int BAND_NUM, int KPOINT_NUM,
	pswf_t* wf, int* fftg, int* labels, double* coords) {

	ppot_t* pps = wf->pps;
	//double complex* x = mkl_calloc(fftg[0]*fftg[1]*fftg[2], sizeof(double complex), 64);
	//printf("START FT\n");
	fft3d(x, wf->G_bounds, wf->lattice, wf->kpts[KPOINT_NUM]->k,
		wf->kpts[KPOINT_NUM]->Gs, wf->kpts[KPOINT_NUM]->bands[BAND_NUM]->Cs,
		wf->kpts[KPOINT_NUM]->bands[BAND_NUM]->num_waves, fftg);
	//printf("FINISH FT\n");
	double* lattice = wf->lattice;
	double vol = determinant(lattice);
	for (int i = 0; i < fftg[0]; i++) {
		double frac[3] = {0,0,0};
		double kdotr = 0;
		for (int j = 0; j < fftg[1]; j++) {
			for (int k = 0; k < fftg[2]; k++) {
				frac[0] = (double) i / fftg[0];
				frac[1] = (double) j / fftg[1];
				frac[2] = (double) k / fftg[2];
				kdotr = dot(wf->kpts[KPOINT_NUM]->k, frac);
				x[i*fftg[1]*fftg[2] + j*fftg[2] + k] *= cexp(2*PI*I*kdotr);
			}
		}
	}

	int num_sites = wf->num_sites;
	#pragma omp parallel for
	for (int p = 0; p < num_sites; p++) {
		projection_t pros = wf->kpts[KPOINT_NUM]->bands[BAND_NUM]->projections[p];
		//printf("READ PROJECTIONS\n");
		ppot_t pp = pps[labels[p]];
		double rmax = pp.wave_grid[pp.wave_gridsize-1];
		double res[3] = {0,0,0};
		vcross(res, lattice+3, lattice+6);
		int grid1 = (int) (mag(res) * rmax / vol * fftg[0]) + 1;
		vcross(res, lattice+0, lattice+6);
		int grid2 = (int) (mag(res) * rmax / vol * fftg[1]) + 1;
		vcross(res, lattice+0, lattice+3);
		int grid3 = (int) (mag(res) * rmax / vol * fftg[2]) + 1;
		int center1 = (int) round(coords[3*p+0] * fftg[0]);
		int center2 = (int) round(coords[3*p+1] * fftg[1]);
		int center3 = (int) round(coords[3*p+2] * fftg[2]);
		//printf("FINISH SETUP %d\n%d %d %d\n%d %d %d\n",p, center1, center2, center3, grid1, grid2, grid3);
		for (int i = -grid1 + center1; i <= grid1 + center1; i++) {
			double frac[3] = {0,0,0};
			double testcoord[3] = {0,0,0};
			int ii=0, jj=0, kk=0;
			double phasecoord[3] = {0,0,0};
			double phase = 0;
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
						phasecoord[0] = coords[3*p+0] + ((ii-i) / fftg[0]);
						phasecoord[1] = coords[3*p+1] + ((jj-j) / fftg[1]);
						phasecoord[2] = coords[3*p+2] + ((kk-k) / fftg[2]);
						phase = dot(phasecoord, wf->kpts[KPOINT_NUM]->k);
						for (int n = 0; n < pros.total_projs; n++) {
							x[ii*fftg[1]*fftg[2] + jj*fftg[2] + kk] +=
								wave_value2(pp.wave_grid,
								pp.funcs[pros.ns[n]].diffwave,
								pp.funcs[pros.ns[n]].diffwave_spline,
								pp.wave_gridsize,
								pros.ls[n], pros.ms[n],
								testcoord)
								* pros.overlaps[n] * cexp(2*PI*I*phase);
								
							//	wave_value(pp.funcs[pros.ns[n]],
							//	pp.wave_gridsize, pp.wave_grid,
							//	pros.ms[n], coords+3*p, frac, lattice)
							//	* pros.overlaps[n] * cexp(2*PI*I*phase);
							//	* Ylm(thetaphi[0], thetaphi[1]);
						}
					}
				}
			}
		}
	}
}

void remove_phase(double complex* x, int KPOINT_NUM, pswf_t* wf, int* fftg) {
	for (int i = 0; i < fftg[0]; i++) {
		double frac[3] = {0,0,0};
		double kdotr = 0;
		for (int j = 0; j < fftg[1]; j++) {
			for (int k = 0; k < fftg[2]; k++) {
				frac[0] = (double) i / fftg[0];
				frac[1] = (double) j / fftg[1];
				frac[2] = (double) k / fftg[2];
				kdotr = dot(wf->kpts[KPOINT_NUM]->k, frac);
				x[i*fftg[1]*fftg[2] + j*fftg[2] + k] *= cexp(-2*PI*I*kdotr);
			}
		}
	}
}

void ncl_realspace_state(double complex* x, int BAND_NUM, int KPOINT_NUM,
	pswf_t* wf, int* fftg, int* labels, double* coords) {

	ppot_t* pps = wf->pps;
	double complex* xup = x;//mkl_calloc(2*fftg[0]*fftg[1]*fftg[2], sizeof(double complex), 64);
	double complex* xdown = x + fftg[0]*fftg[1]*fftg[2];
	int num_waves = wf->kpts[KPOINT_NUM]->bands[BAND_NUM]->num_waves / 2;
	fft3d(xup, wf->G_bounds, wf->lattice, wf->kpts[KPOINT_NUM]->k,
		wf->kpts[KPOINT_NUM]->Gs, wf->kpts[KPOINT_NUM]->bands[BAND_NUM]->Cs,
		num_waves, fftg);
	fft3d(xdown, wf->G_bounds, wf->lattice, wf->kpts[KPOINT_NUM]->k,
		wf->kpts[KPOINT_NUM]->Gs, wf->kpts[KPOINT_NUM]->bands[BAND_NUM]->Cs + num_waves,
		num_waves, fftg);
	double* lattice = wf->lattice;
	double vol = determinant(lattice);
	for (int i = 0; i < fftg[0]; i++) {
		double frac[3] = {0,0,0};
		double kdotr = 0;
		for (int j = 0; j < fftg[1]; j++) {
			for (int k = 0; k < fftg[2]; k++) {
				frac[0] = (double) i / fftg[0];
				frac[1] = (double) j / fftg[1];
				frac[2] = (double) k / fftg[2];
				kdotr = dot(wf->kpts[KPOINT_NUM]->k, frac);
				xup[i*fftg[1]*fftg[2] + j*fftg[2] + k] *= cexp(2*PI*I*kdotr);
				xdown[i*fftg[1]*fftg[2] + j*fftg[2] + k] *= cexp(2*PI*I*kdotr);
			}
		}
	}

	int num_sites = wf->num_sites;
	#pragma omp parallel for
	for (int p = 0; p < num_sites; p++) {
		projection_t up_pros =
			wf->kpts[KPOINT_NUM]->bands[BAND_NUM]->up_projections[p];
		projection_t down_pros =
			wf->kpts[KPOINT_NUM]->bands[BAND_NUM]->down_projections[p];
		ppot_t pp = pps[labels[p]];
		double rmax = pp.wave_grid[pp.wave_gridsize-1];
		double res[3] = {0,0,0};
		vcross(res, lattice+3, lattice+6);
		int grid1 = (int) (mag(res) * rmax / vol * fftg[0]) + 1;
		vcross(res, lattice+0, lattice+6);
		int grid2 = (int) (mag(res) * rmax / vol * fftg[1]) + 1;
		vcross(res, lattice+0, lattice+3);
		int grid3 = (int) (mag(res) * rmax / vol * fftg[2]) + 1;
		int center1 = (int) round(coords[3*p+0] * fftg[0]);
		int center2 = (int) round(coords[3*p+1] * fftg[1]);
		int center3 = (int) round(coords[3*p+2] * fftg[2]);
		for (int i = -grid1 + center1; i <= grid1 + center1; i++) {
			double frac[3] = {0,0,0};
			double testcoord[3] = {0,0,0};
			int ii=0, jj=0, kk=0;
			double phasecoord[3] = {0,0,0};
			double phase = 0;
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
						phasecoord[0] = coords[3*p+0] + ((ii-i) / fftg[0]);
						phasecoord[1] = coords[3*p+1] + ((jj-j) / fftg[1]);
						phasecoord[2] = coords[3*p+2] + ((kk-k) / fftg[2]);
						phase = dot(phasecoord, wf->kpts[KPOINT_NUM]->k);
						for (int n = 0; n < up_pros.total_projs; n++) {
							xup[ii*fftg[1]*fftg[2] + jj*fftg[2] + kk] +=
								wave_value(pp.funcs[up_pros.ns[n]],
								pp.wave_gridsize, pp.wave_grid,
								up_pros.ms[n], coords+3*p, frac, lattice)
								* up_pros.overlaps[n] * cexp(2*PI*I*phase);
							xdown[ii*fftg[1]*fftg[2] + jj*fftg[2] + kk] +=
								wave_value(pp.funcs[down_pros.ns[n]],
								pp.wave_gridsize, pp.wave_grid,
								down_pros.ms[n], coords+3*p, frac, lattice)
								* down_pros.overlaps[n] * cexp(2*PI*I*phase);
						}
					}
				}
			}
		}
	}
}

double* realspace_state_ri(int BAND_NUM, int KPOINT_NUM,
	pswf_t* wf, int* fftg, int* labels, double* coords) {

	int gridsize = fftg[0]*fftg[1]*fftg[2];

	double complex* x = mkl_malloc(gridsize * sizeof(double complex), 64);
	realspace_state(x, BAND_NUM, KPOINT_NUM, wf, fftg, labels, coords);

	double* rpip = (double*) malloc(2 * gridsize * sizeof(double));

	for (int i = 0; i < gridsize; i++) {
		rpip[i] = creal(x[i]);
		rpip[i+gridsize] = cimag(x[i]);
	}
	mkl_free(x);

	return rpip;
}

double* realspace_state_ncl_ri(int BAND_NUM, int KPOINT_NUM,
	pswf_t* wf, int* fftg, int* labels, double* coords) {

	int gridsize = 2*fftg[0]*fftg[1]*fftg[2];

	double complex* x = mkl_malloc(2 * gridsize * sizeof(double complex), 64);
	ncl_realspace_state(x, BAND_NUM, KPOINT_NUM, wf, fftg, labels, coords);

	double* rpip = (double*) malloc(2*gridsize * sizeof(double));

	for (int i = 0; i < gridsize; i++) {
		rpip[i] = creal(x[i]);
		rpip[i+gridsize] = cimag(x[i]);
	}
	mkl_free(x);

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

void write_realspace_state_ncl_ri(char* filename1, char* filename2,
	char* filename3, char* filename4,
	int BAND_NUM, int KPOINT_NUM,
	pswf_t* wf, int* fftg, int* labels, double* coords) {

	int gridsize = fftg[0]*fftg[1]*fftg[2];

	double* x = realspace_state_ncl_ri(BAND_NUM, KPOINT_NUM, wf, fftg, labels, coords);
	write_volumetric(filename1, x+0*gridsize, fftg, 1);
	write_volumetric(filename2, x+1*gridsize, fftg, 1);
	write_volumetric(filename3, x+2*gridsize, fftg, 1);
	write_volumetric(filename4, x+3*gridsize, fftg, 1);

	free(x);
}

double* write_realspace_state_ri_return(char* filename1, char* filename2,
	int BAND_NUM, int KPOINT_NUM,
	pswf_t* wf, int* fftg, int* labels, double* coords) {

	double* x = realspace_state_ri(BAND_NUM, KPOINT_NUM, wf, fftg, labels, coords);
	write_volumetric(filename1, x, fftg, 1);
	write_volumetric(filename2, x+fftg[0]*fftg[1]*fftg[2], fftg, 1);

	return x;
}

double* write_density_return(char* filename, pswf_t* wf,
	int* fftg, int* labels, double* coords) {

	int gridsize = fftg[0] * fftg[1] * fftg[2];
	double* x = mkl_calloc(gridsize, sizeof(double), 64);
	ae_chg_density(x, wf, fftg, labels, coords);
	double scale = determinant(wf->lattice);
	write_volumetric(filename, x, fftg, scale);

	return x;
}

void write_realspace_state_ri_noreturn(char* filename1, char* filename2, int BAND_NUM, int KPOINT_NUM,
	pswf_t* wf, int* fftg, int* labels, double* coords) {

	double* x = write_realspace_state_ri_return(filename1, filename2,
		BAND_NUM, KPOINT_NUM, wf, fftg, labels, coords);
	free(x);
}

void write_density_noreturn(char* filename, pswf_t* wf,
	int* fftg, int* labels, double* coords) {
	setbuf(stdout, NULL);

	double* x = write_density_return(filename, wf, fftg, labels, coords);
	mkl_free(x);
}
