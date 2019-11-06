/** \file
This file contains the functions needed to perform overlap operator evaulations
between bands of different wavefunctions. It uses a similar algorithm to
VASP for the onto_projector and projector_values evaluation.
*/

#ifndef PROJECTOR_H
#define PROJECTOR_H
#include "linalg.h"

/**
Returns a point to a list of ppot_t objects, one for each element in a POTCAR
file. Called as a helper function by Wavefunction.make_c_projectors
*/
ppot_t* get_projector_list(int num_els, int* labels, int* ls, double* wave_grids,
	double* projectors, double* aewaves, double* pswaves, double* rmaxs, double grid_encut);

/**
Finds the coordinates on the FFT grid that fall within each projection sphere
and stores the values of the projectors at that those points, as well
as the coordinates of those points, in real_proj_site_t objects,
one for each site in the structure.
*/
real_proj_site_t* projector_values(int num_sites, int* labels, double* coords,
	double* lattice, double* reclattice, ppot_t* pps, int* fftg);

/**
Finds the coordinates on the FFT grid that fall with each augmentation sphere
(NOT projector sphere) and stores the values of the difference between the partial
waves at those points, as well as the coordinates of those points, in real_proj_site_t
objects, one for each site in the structure.
*/
real_proj_site_t* smooth_pw_values(int num_N, int* Nlst, int* labels, double* coords,
	double* lattice, double* reclattice, ppot_t* pps, int* fftg);

/**
Helper function for onto_projector, which performs the FFT of the wavefunction
into real space and calculates from <p_i|psit_nk> from the grid points found
in projector_values.
*/
void onto_projector_helper(band_t* band, double complex* x, real_proj_site_t* sites,
    int num_sites, double* lattice, double* reclattice, double* kpt, int num_cart_gridpts,
    int* fftg, projection_t* projections);

void get_aug_freqs_helper(band_t* band, double complex* x, real_proj_site_t* sites,
	int num_sites, double* lattice, double* reclattice, double* kpt, int num_cart_gridpts,
	int* fftg, projection_t* projections);

/**
Calculates <p_i|psit_nk> for all i={R,epsilon,l,m} in the structure for one band
*/
void onto_projector(kpoint_t* kpt, int band_num, real_proj_site_t* sites, int num_sites,
	int* G_bounds, double* lattice, double* reclattice, int num_cart_gridpts, int* fftg);

void onto_projector_ncl(kpoint_t* kpt, int band_num, real_proj_site_t* sites, int num_sites,
	int* G_bounds, double* lattice, double* reclattice, int num_cart_gridpts, int* fftg);

/**
Calculates <(phi_i-phit_i)|psit_nk> for all i={R,epsilon,l,m} in the structure for one band
*/
void onto_smoothpw(kpoint_t* kpt, int band_num, real_proj_site_t* sites, int num_sites,
	int* G_bounds, double* lattice, double* reclattice, int num_cart_gridpts, int* fftg);

void get_aug_freqs(kpoint_t* kpt, int band_num, real_proj_site_t* sites, int num_sites,
	int* G_bounds, double* lattice, double* reclattice, int num_cart_gridpts, int* fftg);

/**
Calculates the maximum number of grid points that can be contained within the projector
sphere of pp_ptr given the lattice and fftg (FFT grid dimensions, 3D vector of int),
and stores this value to pp_ptr.
*/
void add_num_cart_gridpts(ppot_t* pp_ptr, double* lattice, int* fftg);

/**
Calculates <phi_i|phi_j>, <phit_i|phit_j>, and <(phi_i-phit_i)|(phi_j-phit_j)>
for all onsite i and j for the element represented by pp_ptr.
*/
void make_pwave_overlap_matrices(ppot_t* pp_ptr);

/**
Evaluates <p_i|psit_nk> for all bands and kpoints of wf.
*/
void setup_projections(pswf_t* wf, ppot_t* pps, int num_elems,
	int num_sites, int* fftg, int* labels, double* coords);

/**
Much more efficient version of overlap_setup.
Calculates three overlap terms for when bands have different
structures:
<(phi1_i-phit1_i)|psit2_n2k>
<(phi2_i-phit2_i)|psit1_n1k>
<(phi1_i-phit1_i)|(phi2_i-phit2_i)>
*/
void overlap_setup_real(pswf_t* wf_R, pswf_t* wf_S,
	int* labels_R, int* labels_S, double* coords_R, double* coords_S,
	int* N_R, int* N_S, int* N_RS_R, int* N_RS_S, int num_N_R, int num_N_S, int num_N_RS);

void overlap_setup_recip(pswf_t* wf_R, pswf_t* wf_S,
	int* labels_R, int* labels_S, double* coords_R, double* coords_S,
	int* N_R, int* N_S, int* N_RS_R, int* N_RS_S, int num_N_R, int num_N_S, int num_N_RS);

/**
Calculates the components of the overlap operator in the augmentation
regions of each ion in the lattice, using the 'aug_real' method.
*/
void compensation_terms(double complex* overlap, int BAND_NUM, pswf_t* wf_S, pswf_t* wf_R,
	int num_M, int num_N_R, int num_N_S, int num_N_RS,
	int* M_R, int* M_S, int* N_R, int* N_S, int* N_RS_R, int* N_RS_S,
	int* proj_labels, double* proj_coords, int* ref_labels, double* ref_coords,
	int* fft_grid, int spin_flip);

/**
Calculates the components of the overlap operator in the augmentation
regions of each ion in the lattice, using the 'aug_recip' method.
*/
void compensation_terms_recip(double complex* overlap, int BAND_NUM, pswf_t* wf_S, pswf_t* wf_R,
	int num_M, int num_N_R, int num_N_S, int num_N_RS,
	int* M_R, int* M_S, int* N_R, int* N_S, int* N_RS_R, int* N_RS_S,
	int* proj_labels, double* proj_coords, int* ref_labels, double* ref_coords,
	int* fft_grid, int spin_flip);

/**
###DEPRECATED###
O(N^2) spherical Bessel transform
*/
double* besselt(double* r, double* k, double* f, double encut, int N, int l);

#endif

