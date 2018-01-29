#ifndef RADIAL_H
#define RADIAL_H

#include <complex.h>

double complex offsite_wave_overlap(double* coord1, double* r1, double* f1, double** spline1, int size1,
	double* coord2, double* r2, double* f2, double** spline2, int size2, double* lattice, int l, int m);

#endif