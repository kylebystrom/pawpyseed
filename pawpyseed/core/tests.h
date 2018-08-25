#ifndef TESTS_H
#define TESTS_H
#include "sbt.h"

double Ylmr(int l, int m, double theta, double phi);
double Ylmi(int l, int m, double theta, double phi);

int fft_check(char* wavecar, double* kpt_weights, int* fftg);

#endif
