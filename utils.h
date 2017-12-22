#ifndef UTILS_H
#define UTILS_H

void vcross(double* res, double* top, double* bottom);

double dot(double* x1, double* x2);

double mag(double* x1);

double determinant(double* m);

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
} funcset_t

typedef struct pseudopot {
	int num_projs;
	funcset_t* funcs;
	int proj_gridsize;
	int wave_gridsize;
	double* wave_grid;
	double* proj_grid;
} ppot_t;

void free_kpoint(kpoint_t* kpt);

void ALLOCATION_FAILED();

#endif
