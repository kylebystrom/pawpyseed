#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <omp.h>
#include "quadrature.h"
#include "utils.h"
#include "radial.h"

#define PI 3.14159265358979323846
#define RHO_PREC 0.015
#define THETA_PREC 0.015
#define R_PREC 0.015

double complex offsite_wave_overlap(double* dcoord, double* r1, double* f1, double** spline1, int size1,
	double* r2, double* f2, double** spline2, int size2,
	double* lattice, int l1, int m1, int l2, int m2) {

	double r1max = r1[size1-1];
	double r2max = r2[size2-1];
	double min_R = 0, max_R = 0;
	double sepdist = mag(dcoord);
	double basis[9] = {0,0,0,0,0,0,0,0,0};
	if (sepdist < 10e-5) {
		basis[0] = 1;
		basis[4] = 1;
		basis[8] = 1;
	} else {
		basis[0] = dcoord[0] / sepdist;
		basis[1] = dcoord[1] / sepdist;
		basis[2] = dcoord[2] / sepdist;
		double temp = pow(basis[1]*basis[1]+basis[0]+basis[0], -0.5);
		basis[3] = -basis[1] * temp;
		basis[4] = basis[0] * temp;
		basis[5] = 0;
		vcross(basis+6, basis+0, basis+3);
	}
	if (sepdist + r2max < r1max) {
		max_R = sepdist + r2max;
	} else {
		max_R = r1max;
	}
	if (sepdist - r2max < -r1max) {
		min_R = -r1max;
	} else {
		min_R = sepdist - r2max;
	}


	double R = min_R - 0.01;
	double dR = 10e-4;
	double theta = 0, dV = 0, rho = 9, drho = 0, dtheta = 0, rhomax = 0;
	int NUM_THETA = 0;
	if (r1max > r2max) {
		rhomax = r1max;
	} else {
		rhomax = r2max;
	}
	double coord1[3] = {0,0,0};
	double coord2[3] = {0,0,0};
	double complex integral = 0;
	while (R < max_R + R_PREC) {
		/*
		if (fabs(R) < 10e-3 || fabs(R-sepdist) < 10e-3)
			dR = 10e-4;
		else {
			dR = 0.1*fabs(R);
			if (dR > R_PREC)
				dR = 0.1*fabs(R-sepdist);
			if (dR > R_PREC)
				dR = R_PREC;
		}
		*/
		dR = 0.01;
		rho = dR;
		drho = dR;
		drho = 0.01;
		while (rho < rhomax + RHO_PREC) {
			NUM_THETA = max(12, (int) (rho * PI / RHO_PREC / 4) * 4);
			dtheta = 2 * PI / NUM_THETA;
			for (int ntheta = 0; ntheta < NUM_THETA; ntheta++) {
				theta = dtheta * ntheta;
				coord1[0] = R * basis[0]
					+ rho * cos(theta) * basis[3]
					+ rho * sin(theta) * basis[6];
				coord1[1] = R * basis[1]
					+ rho * cos(theta) * basis[4]
					+ rho * sin(theta) * basis[7];
				coord1[2] = R * basis[2]
					+ rho * cos(theta) * basis[5]
					+ rho * sin(theta) * basis[8];
				coord2[0] = coord1[0] - dcoord[0];
				coord2[1] = coord1[1] - dcoord[1];
				coord2[2] = coord1[2] - dcoord[2];
				dV = dR * drho * dtheta * rho;
				if (mag(coord1) < r1max && mag(coord2) < r2max) {
					integral += wave_value2(r1, f1, spline1, size1, l1, m1, coord1)
						* wave_value2(r2, f2, spline2, size2, l2, m2, coord2) * dV;
				}
			}
			rho += drho;
			//if (drho < RHO_PREC)
			//	drho = 1.1 * drho;
		}
		R += dR;
	}
	return integral;
}


