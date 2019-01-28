#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <omp.h>
#include "quadrature.h"
#include "utils.h"
#include "radial.h"
#include "gaunt.h"

#define PI 3.14159265358979323846
#define KGRID_SIZE 500

double complex offsite_wave_overlap(double* dcoord,
	double* r1, double* f1, double** spline1, int size1,
	double* r2, double* f2, double** spline2, int size2,
	double* lattice, int l1, int m1, int l2, int m2) {

	double temp[3] = {0,0,0};
	double r1max = r1[size1-1];
	double r2max = r2[size2-1];
	double dphi = 0;
	double THETA, PHI, R1, R2, costheta, dcostheta, sintheta, phi, integral = 0;
	double complex F1, F2;
	int NUM_RADIAL_SUBSTEPS = 1, NUM_PHI = 1, NUM_THETA = 0;

	//loop over 1st coord
	double dr = 0;

	for (int rstep = 0; rstep < size1-1; rstep++) {
		dr = (r1[rstep+1] - r1[rstep]);
		NUM_RADIAL_SUBSTEPS = max(1, (int) (dr/0.01));
		dr /= NUM_RADIAL_SUBSTEPS;
		NUM_THETA = min(80, max(6, (int) (r1[rstep] * PI / 0.01)));
		double* costhetas = QUADRATURE_POINTS[NUM_THETA-3];
		double* dcosthetas = QUADRATURE_WEIGHTS[NUM_THETA-3];
		for (int substep = 0; substep < NUM_RADIAL_SUBSTEPS; substep++) {
			R1 = r1[rstep] + (substep+0.5)*dr;
			for (int thetastep = 0; thetastep < NUM_THETA; thetastep++) {
				costheta = costhetas[thetastep];
				dcostheta = dcosthetas[thetastep];
				sintheta = pow(1 - pow(costheta, 2), 0.5);
				NUM_PHI = max((int) ((NUM_THETA+1)*2*sintheta/4)*4, 12);
				dphi = 2 * PI / NUM_PHI;
				for (int phistep = 0; phistep < NUM_PHI; phistep++) {
					phi = phistep * dphi;
					temp[0] = R1 * sintheta * cos(phi) - dcoord[0];
					temp[1] = R1 * sintheta * sin(phi) - dcoord[1];
					temp[2] = R1 * costheta - dcoord[2];
					R2 = mag(temp);
					if (R2 <= r2max && R1 <= R2 && R2 != 0) {
						THETA = acos(temp[2]/R2);
						if (R2 - fabs(temp[2]) == 0) PHI = 0;
						else PHI = acos(temp[0] / pow(temp[0]*temp[0] + temp[1]*temp[1], 0.5));
						if (temp[1] < 0) PHI = 2*PI - PHI;
						F2 = wave_interpolate(R2, size2, r2, f2, spline2) * Ylm(l2, m2, THETA, PHI);
						//F1 = (f1[rstep] + f1[rstep+1]) / 2 * Ylm2(l1, m1, costheta, phi);
						F1 = wave_interpolate(R1, size1, r1, f1, spline1) * Ylm2(l1, m1, costheta, phi);
						integral += F1 * conj(F2) * R1 / R2 * dr * dcostheta * dphi;
					}
				}
			}
		}
	}

	for (int rstep = 0; rstep < size2-1; rstep++) {
		dr = (r2[rstep+1] - r2[rstep]) / NUM_RADIAL_SUBSTEPS;
		NUM_RADIAL_SUBSTEPS = max(1, (int) (dr/0.01));
		dr /= NUM_RADIAL_SUBSTEPS;
		NUM_THETA = min(100, max(6, (int) (r2[rstep] * PI / 0.01)));
		double* costhetas = QUADRATURE_POINTS[NUM_THETA-3];
		double* dcosthetas = QUADRATURE_WEIGHTS[NUM_THETA-3];
		for (int substep = 0; substep < NUM_RADIAL_SUBSTEPS; substep++) {
			R2 = r2[rstep] + (substep+0.5)*dr;
			for (int thetastep = 0; thetastep < NUM_THETA; thetastep++) {
				costheta = costhetas[thetastep];
				dcostheta = dcosthetas[thetastep];
				sintheta = pow(1 - pow(costheta, 2), 0.5);
				NUM_PHI = max((int) ((NUM_THETA+1)*2*sintheta/4)*4, 12);
				dphi = 2 * PI / NUM_PHI;
				for (int phistep = 0; phistep < NUM_PHI; phistep++) {
					phi = phistep * dphi;
					temp[0] = R2 * sintheta * cos(phi) + dcoord[0];
					temp[1] = R2 * sintheta * sin(phi) + dcoord[1];
					temp[2] = R2 * costheta + dcoord[2];
					R1 = mag(temp);
					if (R1 <= r1max && R2 < R1) {
						THETA = acos(temp[2]/R1);
						if (R1 - fabs(temp[2]) == 0) PHI = 0;
						else PHI = acos(temp[0] / pow(temp[0]*temp[0] + temp[1]*temp[1], 0.5));
						if (temp[1] < 0) PHI = 2*PI - PHI;
						F1 = wave_interpolate(R1, size1, r1, f1, spline1) * Ylm(l1, m1, THETA, PHI);
						F2 = wave_interpolate(R2, size2, r2, f2, spline2) * Ylm2(l2, m2, costheta, phi);
						integral += F1 * conj(F2) * R2 / R1 * dr * dcostheta * dphi;
					}
				}
			}
		}
	}
	return integral;
}

double complex reciprocal_offsite_wave_overlap(double* dcoord,
	double* k1, double* f1, double** s1, int size1,
	double* k2, double* f2, double** s2, int size2,
	double* lattice, int l1, int m1, int l2, int m2) {

	int lx, ly, mx, my;
	if (l1 < l2) {
		lx = l2;
		ly = l1;
		mx = m2;
		my = m1;
	}
	else {
		lx = l1;
		ly = l2;
		mx = m1;
		my = m2;
	}
	if (my < 0) {
		mx = -mx;
		my = -my;
	}
	double kmax = k1[size1-1];
	if (kmax > k2[size2-1]) {
		kmax = k2[size2-1];
	}
	double kmin = k1[0];
	if (kmin < k2[0]) {
		kmin = k2[0];
	}

	double theta = 0, phi = 0;
	double R = mag(dcoord);
	if (R < 10e-12) {
		theta = 0;
		phi = 0;
		R = 0;
	} else {
		theta = acos(dcoord[2]/R);
		if (R - fabs(dcoord[2]) < 10e-12) phi = 0;
		else phi = acos(dcoord[0] / pow(dcoord[0]*dcoord[0] + dcoord[1]*dcoord[1], 0.5));
		if (dcoord[1] < 0) phi = 2*PI - phi;
	}

	double* kgrid = (double*) malloc(KGRID_SIZE * sizeof(double));
	double* ifunc = (double*) malloc(KGRID_SIZE * sizeof(double));
	double complex total = 0;
	double kk;
	double mult_factor = pow(-1, m1) * 8;
	for (int L = abs(l1-l2); L <= l1+l2; L+=2) {
		for (int knum = 0; knum < KGRID_SIZE; knum++) {
			kgrid[knum] = kmin * pow(kmax/kmin, (double) knum / KGRID_SIZE);
			kk = kgrid[knum];
			ifunc[knum] = wave_interpolate(kk, size1, k1, f1, s1)
				* wave_interpolate(kk, size2, k2, f2, s2)
				* kk * kk * sbf(kk*R, L);
		}

		double** ispline = spline_coeff(kgrid, ifunc, KGRID_SIZE);
		if (R > 10e-10)
			total += spline_integral(kgrid, ifunc, ispline, KGRID_SIZE)
				* SBTFACS[lx][ly][(L-abs(l1-l2))/2][lx+mx][my]
				* Ylm(L, m1-m2, theta, phi) * cpow(I, l2+L-l1) * mult_factor;
		else {
			if (L == 0 && l1 == l2 && m1 == m2)
				total += spline_integral(kgrid, ifunc, ispline, KGRID_SIZE)
					* 2 / PI;
		}
		free(ispline[0]);
		free(ispline[1]);
		free(ispline[2]);
		free(ispline);
	}
	free(kgrid);
	free(ifunc);
	return total;
}

/*double complex charge_in_sphere(double* dcoord,
	double* k, double* psf1, double* aef1, int l1, int m1,
	double* psf2, double* aef2, int l2, int m2,
	int size, double rcut) {
	return 0;
}*/

//(sin(k*R)-k*R*cos(k*R))/k^3
//-(k*R*sin(k*R)+2*cos(k*R)-2)/k^3
//double complex reciprocal_offsite_wave_overlap(double* dcoord,
//	double* k1, double* f1, double** s1, int size1,
//	double* k2, double* f2, double** s2, int size2,
//	double* lattice, int l1, int m1, int l2, int m2)
