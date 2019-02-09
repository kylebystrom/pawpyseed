#ifndef TESTS_H
#define TESTS_H
#include "sbt.h"
#include "utils.h"

int fft_check(char* wavecar, double* kpt_weights, int* fftg);

void proj_check(int BAND_NUM, int KPOINT_NUM,
	pswf_t* wf, int* fftg, int* labels, double* coords);

#endif
