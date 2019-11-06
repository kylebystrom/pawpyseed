/** \file
Routines for calculating the overlaps of pseudowavefunctions (i.e. sums of plane-wave coefficients).
*/

#ifndef PSEUDOPROJECTOR_H
#define PSEUDOPROJECTOR_H
#include "utils.h"

/**
DEPRECATED, DO NOT USE
*/
void vc_pseudoprojection(pswf_t* wf_ref, pswf_t* wf_proj, int BAND_NUM, double* results);

/**
Takes two pswf_t objects and a band number BAND_NUM. For each kpoint and spin,
the band BAND_NUM of wf_proj will be projected onto all bands of wf_ref at that
kpoint and spin, i.e. the overlap operators are calculated and returned
as a double complex* res.

The format of the returned array is as follows:
loop over bands
	loops over spins
		loop over kpoints
*/
void pseudoprojection(double complex* projections, pswf_t* wf_ref, pswf_t* wf_proj, int BAND_NUM,
						int flip_spin);

#endif
