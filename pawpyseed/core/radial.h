/** \file
Calculates the overlap integral of two functions, one defined on each
of two overlapping spheres
*/

#ifndef RADIAL_H
#define RADIAL_H

#include <complex.h>

/**
Given dcoord: difference between site locations (R2-R1); r1, the radial grid, size size1, of
the function f1 interpolated with spline1; likewise for f2; the lattice in which the sites
are located, and the angular quantum numbers of the two functions, finds the overlap integral
of f1 and f2.
*/
double complex offsite_wave_overlap(double* dcoord, double* r1, double* f1, double** spline1, int size1,
	double* r2, double* f2, double** spline2, int size2,
	double* lattice, int l1, int m1, int l2, int m2);

/**
Given dcoord: difference between site locations (R2-R1); k1, the reciprocal radial grid,
size size1, of the function f1 interpolated with spline1;
likewise for f2; the lattice in which the sites
are located, and the angular quantum numbers of the two functions, finds the overlap integral
of f1 and f2 in reciprocal space.
Returns integral[f2*(r-R)f1(r)d^3r]
*/
double complex reciprocal_offsite_wave_overlap(double* dcoord,
	double* k1, double* f1, double** s1, int size1,
	double* k2, double* f2, double** s2, int size2,
	double* lattice, int l1, int m1, int l2, int m2);

#endif
