#ifndef PROJECTOR_H
#define PROJECTOR_H
#include <mkl_types.h>

ppot_t* get_projector_list(int num_els, int* labels, int* ls, double* proj_grids, double* wave_grids,
	double* projectors, double* aewaves, double* pswaves, char** rmaxs);

real_proj_site_t* projector_values(int num_sites, int* labels, double* coords,
	double* lattice, double* reclattice, ppot_t* pps, int* fftg);

void onto_projector_helper(band_t* band, MKL_Complex16* x, real_proj_site_t* sites,
	int num_sites, int* labels, double* lattice, double* reclattice, double* kpt, ppot_t* pps, int* fftg);

void onto_projector(kpoint_t* kpt, int band_num, real_proj_site_t* sites, int num_sites, int* labels,
	int* G_bounds, double* lattice, double* reclattice, ppot_t* pps, int* fftg);

void add_num_cart_gridpts(ppot_t* pp_ptr, double* lattice, int* fftg);

void make_pwave_overlap_matrices(ppot_t* pp_ptr);

void setup_projections(pswf_t* wf, ppot_t* pps, int num_elems,
	int num_sites, int* fftg, int* labels, double* coords);

void overlap_setup(pswf_t* wf_R, pswf_t* wf_S, ppot_t* pps,
	int* labels_R, int* labels_S, double* coords_R, double* coords_S,
	int* N_RS_R, int* N_RS_S, int num_N_RS);

double* compensation_terms(int BAND_NUM, pswf_t* wf_proj, pswf_t* wf_ref, ppot_t* pps,
	int num_elems, int num_M, int num_N_R, int num_N_S, int num_N_RS,
	int* M_R, int* M_S, int* N_R, int* N_S, int* N_RS_R, int* N_RS_S,
	int* proj_labels, double* proj_coords, int* ref_labels, double* ref_coords,
	int* fft_grid);

double* besselt(double* r, double* k, double* f, double encut, int N, int l);

#endif

