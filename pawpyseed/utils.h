#ifndef UTILS_H
#define UTILS_H
#include <complex.h>
#include <math.h>
typedef struct band {
	int n;
	int num_waves;
	double occ;
	double N;
	double complex energy;
	float complex* Cs;
} band_t;

typedef struct kpoint {
	short int up;
	int* Gs;
	double* k;
	double weight;
	int num_bands;
	band_t** bands;
} kpoint_t;

typedef struct pswf {
	int* G_bounds;
	kpoint_t** kpts;
	int nspin;
	int nband;
	int nwk;
	double* lattice;
	double* reclattice;
} pswf_t;

typedef struct proj_ae_ps {
	int l;
	double* proj;
	double* aewave;
	double* pswave;
} funcset_t;

typedef struct real_proj {
	int l;
	int m;
	int func_num;
	double complex* values;
} real_proj_t;

typedef struct real_proj_site {
	int index;
	int elem;
	int num_projs;
	int total_projs;
	int num_indices;
	double rmax;
	int* indices;
	real_proj_t* projs;
} real_proj_site_t;

typedef struct pseudopot {
	int num_projs;
	funcset_t* funcs;
	double rmax;
	double* pspw_overlap_matrix;
	double* aepw_overlap_matrix;
	double* diff_overlap_matrix;
	int proj_gridsize;
	int wave_gridsize;
	int num_cart_gridpts;
	double* wave_grid;
	double* proj_grid;
} ppot_t;

void vcross(double* res, double* top, double* bottom);

double dot(double* x1, double* x2);

double mag(double* x1);

double determinant(double* m);

double dist_from_frac(double* coords1, double* coords2, double* lattice);

void frac_to_cartesian(double* coord, double* lattice);

void free_kpoint(kpoint_t* kpt);

void free_pswf(pswf_t* wf);

double* get_occs(pswf_t* wf);

int get_nband(pswf_t* wf);

int get_nwk(pswf_t* wf);

int get_nspin(pswf_t* wf);

double legendre(int l, int m, double x);

double fac(int n);

double complex Ylm(int l, int m, double theta, double phi);

double complex proj_value(funcset_t funcs, int m, double rmax,
	double* ion_pos, double* pos, double* lattice);

void ALLOCATION_FAILED();

#endif
