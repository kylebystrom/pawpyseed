#ifndef DENSITY_H
#define DENSITY_H

double complex* realspace_state(int BAND_NUM, int KPOINT_NUM, pswf_t* wf, ppot_t* pps, int* fftg,
		int* labels, double* coords);

double* ae_chg_density(pswf_t* wf, ppot_t* pps, int* fftg, int* labels, double* coords);

double* project_realspace_state(int BAND_NUM, int numtoproj, pswf_t* wf, pswf_t* wf_R, ppot_t* pps, int* fftg,
	int* labels, double* coords, int* labels_R, double* coords_R);

double* realspace_state_ri(int BAND_NUM, int KPOINT_NUM, pswf_t* wf, ppot_t* pps, int* fftg,
		int* labels, double* coords);

void write_volumetric(char* filename, double* x, int* fftg, double scale);

double* write_realspace_state_ri_return(char* filename1, char* filename2, int BAND_NUM, int KPOINT_NUM,
	pswf_t* wf, ppot_t* pps, int* fftg,
	int* labels, double* coords);

double* write_density_return(char* filename, pswf_t* wf, ppot_t* pps,
	int* fftg, int* labels, double* coords);

void write_realspace_state_ri_noreturn(char* filename1, char* filename2, int BAND_NUM, int KPOINT_NUM,
	pswf_t* wf, ppot_t* pps, int* fftg,
	int* labels, double* coords);

void write_density_noreturn(char* filename, pswf_t* wf, ppot_t* pps,
	int* fftg, int* labels, double* coords);

#endif
