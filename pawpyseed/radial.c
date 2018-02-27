#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <omp.h>
#include "quadrature.h"
#include "utils.h"
#include "radial.h"

#define PI 3.14159265358979323846
#define NUM_THETA 13
#define NUM_PHI 30

double complex offsite_wave_overlap(double* dcoord, double* r1, double* f1, double** spline1, int size1,
	double* r2, double* f2, double** spline2, int size2,
	double* lattice, int l1, int m1, int l2, int m2) {

	double temp[3] = {0,0,0};
	double r1max = r1[size1-1];
	double r2max = r2[size2-1];
	double dphi = 2 * PI / NUM_PHI;
	double THETA, PHI, R1, R2, costheta, dcostheta, sintheta, phi, integral = 0;
	double complex F1, F2;

	//loop over 1st coord
	double dr = 0;
	double* costhetas = QUADRATURE_POINTS[NUM_THETA-3];
	double* dcosthetas = QUADRATURE_WEIGHTS[NUM_THETA-3];

	for (int rstep = 0; rstep < size1-1; rstep++) {
		dr = r1[rstep+1] - r1[rstep];
		R1 = r1[rstep] + dr/2;
		for (int thetastep = 0; thetastep < NUM_THETA; thetastep++) {
			costheta = costhetas[thetastep];
			dcostheta = dcosthetas[thetastep];
			sintheta = pow(1 - pow(costheta, 2), 0.5);
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
					F1 = (f1[rstep] + f1[rstep+1]) / 2 * Ylm2(l1, m1, costheta, phi);
					integral += F1 * conj(F2) * R1 / R2 * dr * dcostheta * dphi;
				}
			}
		}
	}

	for (int rstep = 0; rstep < size2-1; rstep++) {
		dr = r2[rstep+1] - r2[rstep];
		R2 = r2[rstep] + dr/2;
		for (int thetastep = 0; thetastep < NUM_THETA; thetastep++) {
			costheta = costhetas[thetastep];
			dcostheta = dcosthetas[thetastep];
			sintheta = pow(1 - pow(costheta, 2), 0.5);
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
					F2 = (f2[rstep] + f2[rstep]) / 2 * Ylm2(l2, m2, costheta, phi);
					integral += F1 * conj(F2) * R2 / R1 * dr * dcostheta * dphi;
				}
			}
		}
	}
	return integral;
}
