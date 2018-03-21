#ifndef RADIAL_H
#define RADIAL_H

#include <complex.h>

double complex offsite_wave_overlap(double* dcoord, double* r1, double* f1, double** spline1, int size1,
	double* r2, double* f2, double** spline2, int size2,
	double* lattice, int l1, int m1, int l2, int m2);

#endif
