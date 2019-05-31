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

typedef struct WAVECAR_FILE {
	int type;
	FILE* fp;
	char* start;
	char* curr;
} WAVECAR;

WAVECAR* wcopen(char* f, int type);

void wcseek(WAVECAR* wc, long loc);

void wcread(void* ptr0, long size, long nmemb, WAVECAR* wc);

void wcclose(WAVECAR* wc);

/**
Sets up variables to be used to read the pseudowavefunctions from WAVECAR
*/
void setup(int nspin, int nwk, int nband,
	double* nb1, double* nb2, double* nb3, int* np, double ecut,
	double* lattice, double* reclattice);

/**
Handles reading WAVECAR objects, called by read_wavefunctions
and read_wavefunctions_from_str
*/
pswf_t* read_wavecar(WAVECAR* wc, double* kpt_weights);

/**
Given char* filename pointing to a WAVECAR file (VASP output),
constructs a pswf_t* containing the plane-wave coefficients
and energies for each band at each kpoint for each spin.
*/
pswf_t* read_wavefunctions(char* filename, double* kpt_weights);

/**
Read wavefunctions from a string. This is useful if the binary
WAVECAR object is opened from a .gz or .bz2 format by monty
*/
pswf_t* read_wavefunctions_from_str(char* start, double* kpt_weights);

/**
DEPRECATED, DO NOT USE: function to read a single band from a WAVECAR
*/
kpoint_t** read_one_band(int* G_bounds, double* kpt_weights, int* ns, int* nk, int* nb, int BAND_NUM, char* filename);

#endif

