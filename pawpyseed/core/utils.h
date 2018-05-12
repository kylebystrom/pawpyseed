/** \file
Contains several sets of general purpose utilities:

1) structs used by pawpyseed
as well as utility functions for freeing allocated memory
being used by those structs.

2) General vector math (i.e. dot and cross product,
determinant, etc.)

3) Other necessary functionality, such
as spline interpolation, spherical harmonics, and
plane-wave/radial function overlap, that is used
in multiple locations throughout the code.
*/

#ifndef UTILS_H
#define UTILS_H
#include <complex.h>
#include <math.h>

/**
One set of augmentation region functions, i.e.
p_i, phi_i, and phit_i (projector, all electron
partial wave, and partial wave), as well
as the difference of the partial waves, the difference
expanded in spherical Bessel functions, and spline
coefficients for each grid.
*/
typedef struct funcset {
	int l; ///< l quantum number
	double* proj; ///< projector function
	double** proj_spline; ///< projector function spline
	double* aewave; ///< all electron partial wave
	double** aewave_spline; ///< ae partial wave spline coefficients
	double* pswave; ///< pseudo partial wave
	double** pswave_spline; ///< ps partial wave spline coefficients
	double* diffwave; ///< aewave-pswave
	double** diffwave_spline; ///< spline coefficients for diffwave
	double* kwave; ///< Expansion of diffwave in spherical Bessel functions
	double** kwave_spline; ///< spline coefficients for kwave
	double* smooth_diffwave; ///< diffwave on linear grid with high-frequency components removed
	double** smooth_diffwave_spline; ///< spline coefficients for smooth_diffwave
	double* dense_kwave;
	double** dense_kwave_spline;
} funcset_t;

typedef struct ppot {
	int num_projs; ///< number of radial projector functions
	int total_projs; ///< number of projector functions
	int lmax; ///< maximum l-value of any projector
	funcset_t* funcs; ///< funcset for each projector, see funcset
	double rmax; //< maximum radius of the projector functions
	double wave_rmax; ///< maximum radius of the partial waves
	double* pspw_overlap_matrix; ///< overlap matrix for pseudo partial waves
	double* aepw_overlap_matrix; ///< overlap matrix for all electron partial waves
	double* diff_overlap_matrix; ///< overlap matrix of difference between all electron and partial waves
	int proj_gridsize; ///< number of points on projector radial grid
	int wave_gridsize; ///< number of points on partial wave radial grid
	int num_cart_gridpts; ///< number of real space grid points that can fit in the projector sphere
	double* wave_grid; ///< real radial grid for partial waves
	double* kwave_grid; ///< reciprocal radial grid for partial waves
	double* proj_grid; ///< real radial grid for projector functions
	double* smooth_grid;
	double* dense_kgrid;
} ppot_t;

typedef struct projection {
	int num_projs; ///< number of radial projector functions
	int total_projs; ///< number of projector functions
	int* ns; ///< radial projector index
	int* ls; ///< l values of projectors
	int* ms; ///< m values of projectors
	double complex* overlaps; ///< list of <p_i|psi>
} projection_t;

/**
Stores the data for a single
band, or Kohn Sham single particle
state, for a structure
*/
typedef struct band {
	int n; ///< band number
	int num_waves; ///< number of plane waves
	double occ; ///< occupancy of the band
	double N;
	double complex energy; ///< energy of the band
	float complex* Cs; ///< plane wave coefficients (normalized to 1)
	double complex* CRs; ///< wavefunction in real space
	projection_t* projections; ///< length==number of sites in structure
	projection_t* wave_projections; ///< length==number of sites in structure
} band_t;

typedef struct rayleigh_set {
	int l; ///< angular momentum quantum number
	double complex* terms; ///< rayleigh epansion terms
} rayleigh_set_t;

typedef struct kpoint {
	short int up; ///< spin
	int num_waves; ///< number of plane waves in a band
	int* Gs; ///< plane wave coefficients, in sets of three
	double* k; ///< k-point vector
	double weight; ///< k-point weight
	int num_bands; ///< number of bands
	band_t** bands; ///< bands with this k-point
	rayleigh_set_t** expansion;
} kpoint_t;

typedef struct pswf {
	int num_elems;
	int num_sites;
	ppot_t* pps;
	int* G_bounds;
	kpoint_t** kpts;
	int nspin;
	int nband;
	int nwk;
	double* lattice;
	double* reclattice;
	int* fftg;
	int num_aug_overlap_sites;
	double* dcoords;
	double complex** overlaps;
} pswf_t;

typedef struct projgrid {
	double complex* values;
} projgrid_t;

typedef struct real_proj {
	int l;
	int m;
	int func_num;
	double* paths;
	double complex* values;
} real_proj_t;

typedef struct real_proj_site {
	int index;
	int elem;
	int num_projs;
	int total_projs;
	int num_indices;
	int gridsize;
	double rmax;
	double* coord;
	int* indices;
	real_proj_t* projs;
} real_proj_site_t;

/**
Minimum of a and b
*/
int min(int a, int b);

/**
Maximum of a and b
*/
int max(int a, int b);

/**
Calculates the cross product top X bottom and places
it in res, assuming that top and bottom are both
length three
*/
void vcross(double* res, double* top, double* bottom);

/**
Returns the dot product of x1 and x2 assuming
that x1 and x2 are both of length three
*/
double dot(double* x1, double* x2);

/**
Magnitude of a three-element vector
*/
double mag(double* x1);

/**
Calculates the determinant of the matrix m with the
precondition m is 3x3 (length 9) and in row major order.
*/
double determinant(double* m);

/**
Given
coords1: fractional coordinate in lattice (length 3)
coords2: fractional coordinate in lattice (length 3)
lattice: row major lattice matrix, where each row
is a lattice unit vector (length 9),
returns the distance between coords1 and coords2.
*/
double dist_from_frac(double* coords1, double* coords2, double* lattice);

void frac_to_cartesian(double* coord, double* lattice);

void cartesian_to_frac(double* coord, double* reclattice);

void min_cart_path(double* coord, double* center, double* lattice, double* path, double* r);

double complex trilinear_interpolate(double complex* c, double* frac, int* fftg);

void free_kpoint(kpoint_t* kpt, int num_elems, ppot_t* pps);

void free_ppot(ppot_t* pp);

void free_real_proj(real_proj_t* proj);

void free_real_proj_site(real_proj_site_t* site);

void free_pswf(pswf_t* wf);

void free_ptr(void* ptr);

void free_real_proj_site_list(real_proj_site_t* sites, int length);

void free_ppot_list(ppot_t* pps, int length);

/**
Returns a list with the occupation of each band
of the wavefunction at each kpoint as shown:
loop over bands
	loop over spins
		loop over kpoints
*/
double* get_occs(pswf_t* wf);

/** Return the number of bands in the wavefunction. */
int get_nband(pswf_t* wf);

/** Return the number of kpoins in the wavefunction. */
int get_nwk(pswf_t* wf);

/** Return the number of spins in the wavefunction. */
int get_nspin(pswf_t* wf);

/** Sets the numer of sites in the structure for the wavefunction. */
void set_num_sites(pswf_t* wf, int nsites);

/** Associated legendre polynomial P_lm(x) */
double legendre(int l, int m, double x);

/** factorial */
double fac(int n);

/**
Returns the value of the complex spherical harmonic
assuming physics notation (phi is azimuthal angle)
*/
double complex Ylm(int l, int m, double theta, double phi);

/**
Same as Ylm but takes the cosine of theta instead of the angle itself.
*/
double complex Ylm2(int l, int m, double costheta, double phi);

/**
Interpolate the discretely defined projector function proj, defined on linear radial
grid x, at radius r. rmax is the maximum radius of the projector.
Uses spline interpolation, where proj_spline is the set of spline coefficients
for proj set up by spline_coeff
*/
double proj_interpolate(double r, double rmax, int size, double* x, double* proj, double** proj_spline);

/**
Interpolate the discretely defined partial wave f, defined on logarithmic radial grix x,
at radius r. Uses spline interpolation, where wave_spline is the set of spline
coefficients for f set up by spline_coeff
*/
double wave_interpolate(double r, int size, double* x, double* f, double** wave_spline);

/**
Helper function for proj_value and smooth_wave_value
*/
double complex proj_value_helper(double r, double rmax, int size,
	double* pos, double* x, double* f, double** s, int l, int m);

/**
Return the value of funcs->proj, defined on linear radial grid x centered
at 3D vector ion_pos, at position pos, given the real space lattice.
*/
double complex proj_value(funcset_t funcs, double* x, int m, double rmax,
	int size, double* ion_pos, double* pos, double* lattice);

/**
Return the value of funcs->smooth_diffwave, defined on linear radial grid x
centered at 3D vector ion_pos, at position pos, given the real space lattice
*/
double complex smooth_wave_value(funcset_t funcs, double* x, int m, double rmax,
	int size, double* ion_pos, double* pos, double* lattice);

/**
Return the value of funcs->aewave-funcs->pswave, defined on logarithmic radial
grid x center at 3D vector ion_pos, at position pos, given the real space lattice.
*/
double complex wave_value(funcset_t funcs, int size, double* x, int m,
        double* ion_pos, double* pos, double* lattice);

/**
Interpolates the value of a function on a logarithmic radial grid
at location pos.
*/
double complex wave_value2(double* x, double* wave, double** spline, int size,
	int l, int m, double* pos);

/**
Convenience function for setting up real_proj_site_t* lists
*/
void setup_site(real_proj_site_t* sites, ppot_t* pps, int num_sites, int* site_nums,
    int* labels, double* coords, double* lattice, int* fftg, int pr0_pw1);

/**
Set up spline coefficients for spline interpolation.
Essentially a translation into C of the VASP SPLCOF function.
G. Kresse and J. Hafner. Ab initio molecular dynamics for liquid metals.
Phys. Rev. B, 47:558, 1993.
*/
double** spline_coeff(double* x, double* y, int N);

/**
Find the integral of a discretely defined function that is
fitted with spline coefficients at each point. x is the independent
variable grid, a is the function values, s is the spline, and size is the
number of points on the grid.
*/
double spline_integral(double* x, double* a, double** s, int size);

void frac_from_index(int index, double* coord, int* fftg);

double sph_bessel(double k, double r, int l);

/**
Return the value of the l-order spherical bessel function at x.
*/
double sbf(double x, int l);

/**
Calculates <(phi_i-phit_i)|psit_nk> for one pseudowavefunction
band psit_nk using the projections onto plane waves calculated
in generate_rayleigh_expansion_terms
*/
double complex rayexp(double* kpt, int* Gs, float complex* Cs, int l, int m,
        int num_waves, double complex* sum_terms, double* ionp);

/**
Calculates <(phi_i-phit_i)|(k+G)> for a specific index i by
using the Rayleigh expansion of a plane wave.
*/
double complex* rayexp_terms(double* kpt, int* Gs, int num_waves,
        int l, int wave_gridsize, double* grid,
        double* wave, double** spline, double* reclattice);

/**
Calculates the overlaps <(phi_i-phit_i)|(k+G)>, where phi_i are the
AE partial waves, phit_i are the PS partial waves, and k+G are the plane
waves, in reciprocal space.
*/
void generate_rayleigh_expansion_terms(pswf_t* wf, ppot_t* pps, int num_elems);

/**
Copy the Rayleigh expansion terms from wf_R to wf.
*/
void copy_rayleigh_expansion_terms(pswf_t* wf, ppot_t* pps, int num_elems, pswf_t* wf_R);

/**
Called after a malloc or calloc call to check that
the allocation was successful.
*/
void CHECK_ALLOCATION(void* ptr);

/**
Called when a memory allocation fails to exit the program
*/
void ALLOCATION_FAILED();

#endif
