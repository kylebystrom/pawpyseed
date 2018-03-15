#ifndef READER_H
#define READER_H
#include "utils.h"

void setup(char* filename, int* pnrecl, int* pnspin, int* pnwk, int* pnband,
	double* nb1, double* nb2, double* nb3, double* ecut,
        double* lattice, double* reclattice);

pswf_t* read_wavefunctions(char* filename, double* kpt_weights);

kpoint_t** read_one_band(int* G_bounds, double* kpt_weights, int* ns, int* nk, int* nb, int BAND_NUM, char* filename);

#endif

