#ifndef SBT_H
#define SBT_H

/*
The following routines are based on the Fortran program NumSBT written by J. Talman.
The algorithm performs a spherical Bessel transform in O(NlnN) time. If you adapt
this code for any purpose, please cite:
Talman, J. Computer Physics Communications 2009, 180, 332 �~@~S338.
The code is distributed under the Standard CPC license.
*/

typedef struct sbt_setup {
        double kmin;
        double kappamin;
        double rmin;
        double rhomin;
        double drho; //drho == dkappa
        double dt;
        int N;
        double complex** mult_table;
        double* ks;
        double* rs;
} sbt_descriptor_t;

sbt_descriptor_t* spherical_bessel_transform_setup(double encut, double enbuf, int lmax, int N,
	double* r, double* ks);

double* wave_spherical_bessel_transform(sbt_descriptor_t* d, double* f, int l);

#endif
