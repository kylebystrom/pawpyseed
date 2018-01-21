#ifndef PROJECTOR_H
#define PROJECTOR_H

void vc_pseudoprojection(pswf_t* wf_ref, pswf_t* wf_proj, int BAND_NUM, double* results);

double* pseudoprojection(pswf_t* wf_ref, pswf_t* wf_proj, int BAND_NUM);

ppot_t* get_projector_list(int num_els, int* labels, int* ls, double* proj_grids, double* wave_grids,
	double* projectors, double* aewaves, double* pswaves, double* rmaxs);

real_proj_site_t* projector_values(int num_sites, int* labels, double* coords,
	double* lattice, double* reclattice, ppot_t* pps, int* fftg);

double complex* onto_projector(real_proj_site_t* sites, int* labels, int* G_bounds, double* lattice,
	double* k, int* Gs, float complex* Cs, int num_waves, int num_M, int* M, ppot_t* pps, int* fftg);

void add_num_cart_gridpts(ppot_t* pp_ptr, double* lattice, int* fftg);

void make_pwave_overlap_matrices(ppot_t* pp_ptr);

double* compensation_terms(int BAND_NUM, pswf_t* wf_proj, pswf_t* wf_ref, ppot_t* pps,
	int num_elems, int num_M, int num_N_R, int num_N_S, int num_N_RS,
	int* M, int* N_R, int* N_S, int* N_RS,
	int* proj_labels, double* proj_coords, int* ref_labels, double* ref_coords, int* fft_grid);

#endif

