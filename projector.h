#ifndef PROJECTOR_H
#define PROJECTOR_H

void vc_pseudoprojection(pswf_t* wf_ref, pswf_t* wf_proj, int BAND_NUM, double* results);

float* pseudoprojection(pswf_t* wf_ref, pswf_t* wf_proj, int BAND_NUM);

ppot_t* get_projector_list(int num_els, int* labels, int* ls, double* proj_grids, double* wave_grids,
	double* projectors, double* aewaves, double* pswaves);

double* onto_projector(int* labels, double* coords, int* G_bounds, double* lattice,
	double* Gs, double* Cs, int num_waves, int num_M, int* M, ppot_t* pps, int* fftg);

double* make_pwave_overlap_matrices(ppot_t pp);

double* compensation_terms(pswf* wf_proj, pswf* wf_ref, ppot_t* pps,
	int num_M, int num_N_R, int num_N_S, int num_N_RS,
	int* M, int* N_R, int* N_S, int* N_RS,
	int* proj_labels, double* proj_coords, int* ref_labels, double* ref_coords, int* fft_grid);

#endif

