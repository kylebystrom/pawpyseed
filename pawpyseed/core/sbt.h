/** \file
The following routines are based on the Fortran program NumSBT written by J. Talman.
The algorithm performs a spherical Bessel transform in O(NlnN) time. If you adapt
this code for any purpose, please cite:
Talman, J. Computer Physics Communications 2009, 180, 332-338.
The code is distributed under the Standard CPC license.
*/

#ifndef SBT_H
#define SBT_H

#include <complex.h>

/**
Contains the parameters for a fast spherical Bessel transform.
*/
typedef struct sbt_descriptor {
        double kmin; ///< Minimum reciprocal space value
        double kappamin; ///< ln(kmin)
        double rmin; ///< Minimum real space value
        double rhomin; ///< ln(rmin)
        double drho; ///< linear increment of rho = ln(r), drho == dkappa
        double dt; ///< increment of the multiplication table
        int N; ///< number of values of r and k
        double complex** mult_table; ///< M_l(t) for l up to lmax
        double* ks; ///< Reciprocal space grid
        double* rs; ///< Real space grid
		int lmax;
} sbt_descriptor_t;

/**
Creates an sbt_descriptor_t object from the real space and k-space grids.
Extrapolates the radial grid to have N values lower than the initial rmin,
and the reciprocal grid to have N values higher than the initial kmin.
encut+enbuf is the maximum value of k BEFORE extrapolation.
lmax is the maximum l-value used to construct the the mult_table.
N is the size of the grid BEFORE exprapolation (i.e. length of r and ks).
*/
sbt_descriptor_t* spherical_bessel_transform_setup(double encut, double enbuf, int lmax, int N,
	double* r, double* ks);

/**
Returns the values of radially defined function f
in reciprocal space in a spherical bessel basis.
f is defined on the grid corresponding to the r
passed into spherical_bessesl_transform_setup.
*/
double* wave_spherical_bessel_transform(sbt_descriptor_t* d, double* f, int l);

/**
Performs the inverse transform of f. Calling this function on the output
of wave_spherical_bessel_transform gives the initial input to
wave_spherical_bessel_transform on the initial radial grid passed to
spherical_bessel_transform_setup
*/
double* inverse_wave_spherical_bessel_transform(sbt_descriptor_t* d, double* f, int l);

/**
Free all the memory associated with an sbt descriptor.
*/
void free_sbt_descriptor(sbt_descriptor_t* d);

#endif
