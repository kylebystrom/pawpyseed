/** \file
Functions used to read pswf_t structs from WAVECAR files (VASP output).
This code is based upon the Fortran program, WaveTrans, written by
R. M. Feenstra and M. Widom from the Dept. of Physics at Carnegie
Mellon University. To see the original work, please visit:
https://www.andrew.cmu.edu/user/feenstra/wavetrans/
*/

#ifndef READER_H
#define READER_H
#include "utils.h"

/**
Sets up variables to be used to read the pseudowavefunctions from WAVECAR
*/
void setup(char* filename, int* pnrecl, int* pnspin, int* pnwk, int* pnband,
	double* nb1, double* nb2, double* nb3, double* ecut,
        double* lattice, double* reclattice);

/**
Given char* filename pointing to a WAVECAR file (VASP output),
constructs a pswf_t* containing the plane-wave coefficients
and energies for each band at each kpoint for each spin.
*/
pswf_t* read_wavefunctions(char* filename, double* kpt_weights);

/**
DEPRECATED, DO NOT USE: function to read a single band from a WAVECAR
*/
kpoint_t** read_one_band(int* G_bounds, double* kpt_weights, int* ns, int* nk, int* nb, int BAND_NUM, char* filename);

#endif

