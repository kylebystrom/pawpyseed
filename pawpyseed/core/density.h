/** \file
Functions for calculating the wavefunction and charge density in real space.
Also handles I/O of Kohn-Sham states and all electron (AE) charge density.
*/

#ifndef DENSITY_H
#define DENSITY_H

/**
Calculates the AE Kohn Sham state of band BAND_NUM at kpoint KPOINT_NUM in real space,
on fractional coordinate real-space grid fftg. x is the slow index.
*/
void realspace_state(double complex* x, int BAND_NUM, int KPOINT_NUM,
	pswf_t* wf, int* fftg, int* labels, double* coords);

void remove_phase(double complex* x, int KPOINT_NUM, pswf_t* wf, int* fftg);

void ae_state_density(double* P, int BAND_NUM, int KPOINT_NUM, pswf_t* wf,
	int* fftg, int* labels, double* coords);

/**
Calculates the AE Kohn Sham state of band BAND_NUM at kpoint KPOINT_NUM in real space,
on fractional coordinate real-space grid fftg, for a noncollinear VASP run.
x is the slow index.
*/
void ncl_realspace_state(double complex* x, int BAND_NUM, int KPOINT_NUM,
	pswf_t* wf, int* fftg, int* labels, double* coords);

/**
Calculates the all electron charge density by adding up realspace_state for all the bands
at each kpoint. Equivalent to the grid in AECCAR of VASP except x is the slow index instead of z.
*/
void ae_chg_density(double* P, pswf_t* wf, int* fftg, int* labels, double* coords);
void ncl_ae_chg_density(double* P, pswf_t* wf, int* fftg, int* labels, double* coords);

/**
Projects one band of wf onto all the bands of wf_R in real space. Very slow for large systems,
but a good test tool.
*/
void project_realspace_state(double complex* projs, int BAND_NUM, pswf_t* wf, pswf_t* wf_R,
	int* fftg, int* labels, double* coords, int* labels_R, double* coords_R);

void write_realspace_state_ncl_ri(char* filename1, char* filename2,
	char* filename3, char* filename4, int BAND_NUM, int KPOINT_NUM,
	pswf_t* wf, int* fftg, int* labels, double* coords);

/**
Separates the result of realspace_state into a real part and an imaginary part,
then returns it. Utilizied by the Python code.
*/
double* realspace_state_ri(int BAND_NUM, int KPOINT_NUM, pswf_t* wf, int* fftg,
		int* labels, double* coords);

/**
Writes a volumetric dataset for a system stored on grid fftg to a file called filename.
Takes in a volumetric dataset where x is the slow index, and prints out a volumetric
dataset where z is the slow index, to match VASP output for compatibility with
programs like VESTA.
*/
void write_volumetric(char* filename, double* x, int* fftg, double scale);

/**
Writes the real and imaginary parts of a Kohn Sham state two files called filename1
and filename2, and returns the (x slow index) array containing the real part followed
by the imaginary part.
*/
double* write_realspace_state_ri_return(char* filename1, char* filename2, int BAND_NUM, int KPOINT_NUM,
	pswf_t* wf, int* fftg, int* labels, double* coords);

/**
Writes the all electron charge density to a file called filename, and returns the
(x slow index) array containing the charge density.
*/
double* write_density_return(char* filename, pswf_t* wf,
	int* fftg, int* labels, double* coords);

/**
Calls write_realspace_state_ri_return and then frees the array
*/
void write_realspace_state_ri_noreturn(char* filename1, char* filename2, int BAND_NUM, int KPOINT_NUM,
	pswf_t* wf, int* fftg, int* labels, double* coords);

/**
Calls write_density_return and then frees the array
*/
void write_density_noreturn(char* filename, pswf_t* wf,
	int* fftg, int* labels, double* coords);

#endif
