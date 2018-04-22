#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <omp.h>
#include "quadrature.h"
#include "utils.h"
#include "radial.h"

#define PI 3.14159265358979323846
#define RHO_PREC 0.008

double complex Ylm_rad(int l, int m, double* r) {
	double R = mag(r);
	if (R < 10e-12) return 0;

	double costheta = 0, phi = 0;
	costheta = temp[2]/r;
	if (r - fabs(r[2]) < 10e-12) phi = 0;
	else phi = acos(r[0] / pow(r[0]*r[0] + r[1]*r[1], 0.5));
	if (r[1] < 0) phi = 2*PI - phi;
	return Ylm2(l, m, costheta, phi);
}

double tfac(int n, int l) {
	double total = 1;
	for (int i = 2*n; i > 0; i-=2) {
		total *= i;
	}
	if (2*n-2*l-1 > 0) {
		for (int i = 2*n-2*l-1; i > 0; i-=2) {
			total *= i;
		}
	} else if (2*n-2*l-1 < -1) {
		total *= pow(-1, n);
		for (int i = 2*l-2*n-1; i > 0; i-=2) {
			total /= i;
		}
	}
	return total;
}

double Pkern1(int l1, int l2, int l3, double r1, double r2, double r3) {
	double lambda = 0.5 * (l1+l2+l3);
	double total = 0;
	for (int n1 = 0; n1 <= lambda; n1++) {
		for (int n2 = 0; n2 <= lambda - n1; n2++) {
			int n3 = lambda - n2 - n1;
			total += pow(r1, 2*n1-l1-1) * pow(r2, 2*n2-l2-1) * pow(r3, 2*n3-l3-1)
				* tfac(n1,l2) * tfac(n2,l2) * tfac(n3,l3);
		}
	}
	return total;
}

double LKern1(int l1, int l2, int l3, double r1, double r2, double R) {
	if (r1+r2-R > 0 && r1-r2+R > 0 && r2-r1+R > 0) {
		return 2 * PI * Pkern1(l1, l2, l3, r1, r2, R);
	} else {
		return 0;
	}
}

double complex offsite_wave_overlap(double* dcoord, double* r1, double* f1, double** spline1, int size1,
	double* r2, double* f2, double** spline2, int size2,
	double* lattice, int l1, int m1, int l2, int m2) {

	double R = mag(dcoord);
	double complex total = 0;

	for (int i = 0; i < size1-1; i++) {
		dx = r1[i+1] - r1[i];
		//Ii = dx * (f1[i] + dx * (spline1[0][i]/2 + dx * (spline1[1][i]/3 + spline1[2][i]*dx/4)));
		Ii = dx * dx * (f1[i]/2 + dx * (spline1[0][i]/3 + dx * (spline1[1][i]/4 + spline1[2][i]*dx/5)));
		ri = (r1[i+1] + r1[i]) / 2;
		for (int j = 0; j < size2-1; j++) {
			dy = r2[j+1] - r2[j];
			//Ij = dy * (f2[j] + dy * (spline2[0][j]/2 + dy * (spline2[1][j]/3 + spline2[2][j]*dx/4)));
			Ij = dy * dy * (f2[j]/2 + dy * (spline2[0][j]/3 + dy * (spline2[1][j]/4 + spline2[2][j]*dx/5)));
			rj = (r2[j+1] + r2[j]) / 2;
			double complex subtotal = 0;
			for (int l = abs(l1-l2); l <= l1+l2; l++) {
				for (int m = -l; m <= l; m++) {
					subtotal += GAUNT_COEFF[l1][m1][l2][m2][l][m] * Ylm_rad(dcoord) * Lkern1(ri, rj, R);
				}
			}
			total += subtotal * Ii * Ij;
		}
	}

	return total * pow(-1, l1);
}