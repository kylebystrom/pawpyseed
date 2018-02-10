#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <omp.h>
#include "quadrature.h"
#include "utils.h"
#include "radial.h"

#define PI 3.14159265358979323846

double complex offsite_wave_overlap(double* coord1, double* r1, double* f1, double** spline1, int size1,
	double* coord2, double* r2, double* f2, double** spline2, int size2,
	double* lattice, int l1, int m1, int l2, int m2) {

	double dcoord[3] = {0,0,0};
	double temp[3] = {0,0,0};
	double R = 0;
	min_cart_path(coord2, coord1, lattice, dcoord, &R);
	double r1max = r1[size1-1];
	double r2max = r2[size2-1];
	double dphi = PI / 30;
	double THETA, PHI, R1, R2, costheta, dcostheta, sintheta, phi, integral = 0;
	double complex F1, F2;

	//loop over 1st coord
	double dr = r1[0];
	double* costhetas = QUADRATURE_POINTS[26];
	double* dcosthetas = QUADRATURE_WEIGHTS[26];
	for (int rstep = 0; rstep < size1; rstep++) {
		R1 = r1[rstep];
		for (int thetastep = 0; thetastep < 29; thetastep++) {
			costheta = costhetas[thetastep];
			dcostheta = dcosthetas[thetastep];
			sintheta = pow(1 - pow(costheta, 2), 0.5);
			for (int phistep = 0; phistep < 60; phistep++) {
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
					F1 = f1[rstep] * Ylm2(l1, m1, costheta, phi);
					integral += conj(F1) * F2 * R1 / R2 * dr * dcostheta * dphi;
				}
			}
		}
		if (rstep+1 != size1) dr = r1[rstep+1] - r1[rstep];
	}

	dr = r2[0];
	for (int rstep = 0; rstep < size2; rstep++) {
		R2 = r2[rstep];
		for (int thetastep = 0; thetastep < 29; thetastep++) {
			costheta = costhetas[thetastep];
			dcostheta = dcosthetas[thetastep];
			sintheta = pow(1 - pow(costheta, 2), 0.5);
			for (int phistep = 0; phistep < 60; phistep++) {
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
					F2 = f2[rstep] * Ylm2(l2, m2, costheta, phi);
					integral += conj(F1) * F2 * R2 / R1 * dr * dcostheta * dphi;
				}
			}
		}
		if (rstep+1 != size2) dr = r2[rstep+1] - r2[rstep];
	}
	return integral;
}
