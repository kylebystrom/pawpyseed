#ifndef UTILS_H
#define UTILS_H
#include <complex.h>
#include <math.h>

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
	projection_t* projections; ///< length==number of sites in structure
	projection_t* wave_projections; ///< length==number of sites in structure
} band_t;

typedef struct rayleigh_set {
	int l;
	double complex* terms;
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
	int* G_bounds;
	kpoint_t** kpts;
	int nspin;
	int nband;
	int nwk;
	double* lattice;
	double* reclattice;
	int num_aug_overlap_sites;
	double complex** overlaps;
} pswf_t;

typedef struct funcset {
	int l;
	double* proj;
	double** proj_spline;
	double* aewave;
	double** aewave_spline;
	double* pswave;
	double** pswave_spline;
	double* diffwave;
	double** diffwave_spline;
	double* kwave;
	double** kwave_spline;
} funcset_t;

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
	double rmax;
	double* coord;
	int* indices;
	real_proj_t* projs;
} real_proj_site_t;

typedef struct ppot {
	int num_projs;
	int total_projs;
	int lmax;
	funcset_t* funcs;
	double rmax;
	double* pspw_overlap_matrix;
	double* aepw_overlap_matrix;
	double* diff_overlap_matrix;
	int proj_gridsize;
	int wave_gridsize;
	int num_cart_gridpts;
	double* wave_grid;
	double* kwave_grid;
	double* proj_grid;
} ppot_t;

/**
Minimum of a and b
*/
int min(int a, int b);

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

void free_kpoint(kpoint_t* kpt);

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
double proj_interpolate(double r, double rmax, double* x, double* proj, double** proj_spline);

/**
Interpolate the discretely defined partial wave f, defined on logarithmic radial grix x,
at radius r. Uses spline interpolation, where wave_spline is the set of spline
coefficients for f set up by spline_coeff
*/
double wave_interpolate(double r, int size, double* x, double* f, double** wave_spline);

/**
Return the value of funcs->proj, defined on linear radial grid x centered
at 3D vector ion_pos, at position pos, giving the real space lattice.
*/
double complex proj_value(funcset_t funcs, double* x, int m, double rmax,
	double* ion_pos, double* pos, double* lattice);

/**
Set up spline coefficients for spline interpolation.
Essentially a translation into C of the VASP SPLCOF function.
G. Kresse and J. Hafner. Ab initio molecular dynamics for liquid metals.
Phys. Rev. B, 47:558, 1993.
*/
double** spline_coeff(double* x, double* y, int N);

void frac_from_index(int index, double* coord, int* fftg);

double sph_bessel(double k, double r, int l);

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
Called after a malloc or calloc call to check that
the allocation was successful.
*/
void CHECK_ALLOCATION(void* ptr);

/**
Called when a memory allocation fails to exit the program
*/
void ALLOCATION_FAILED();

#endif
