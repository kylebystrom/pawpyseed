#include <complex.h>
#include <math.h>
#include <mkl.h>
#include <mkl_types.h>
#include "condon_shortley.h"
#include "../core/utils.h"

//double legendre(int l, int m, double x) {
//	double l0 = 
//
//}

double laguerre(double x, int n, int alpha) {
	double l0 = 1;
	if (n == 0) return l0;
	double l1 = 1 + alpha - x;
	double lcurr = l1;
	for (int k = 1; k < n; k++) {
		l1 = lcurr;
		lcurr = ((2*k+1+alpha-x) * l1 - (k+alpha) * l0) / (k+1);
		l0 = l1;
	}
	return lcurr;
}

double hradial(int n, int l, double r) {
	return pow(4.0/n/n/n * fac(n-l-1)/fac(n+l), 0.5) * exp(-r/2/n) * pow(r/n, l) * laguerre(r, n-l-1, 2*l+1);
}

double* yk(int k, int size, double* r, double* P1, double* P2) {
	double* integrand1 = (double*) malloc(size * sizeof(double));
	double* integrand2 = (double*) malloc(size * sizeof(double));

	for (int i = 0; i < size; i++) {
		integrand1[i] = P1[i] * P2[i] * pow(r[i], k); //might need to change to k-1
	}
	for (int i = 0; i < size; i++) {
		integrand2[i] = P1[i] * P2[i] * pow(r[i], -k-1); //might need to change to -k-2
	}
	double* integrals = ykint(k, size, r, integrand1, integrand2);

	free(integrand1);
	free(integrand2);
}

double* ykchg(int k, int size, double* r, double* chg) {
	double* integrand1 = (double*) malloc(size * sizeof(double));
	double* integrand2 = (double*) malloc(size * sizeof(double));

	for (int i = 0; i < size; i++) {
		integrand1[i] = chg[i] * pow(r[i], k); //might need to change to k-1
	}
	for (int i = 0; i < size; i++) {
		integrand2[i] = chg[i] * pow(r[i], -k-1); //might need to change to -k-2
	}
	double* integrals = ykint(k, size, r, integrand1, integrand2);

	free(integrand1);
	free(integrand2);
}

double* dft_ecoul(int size, double* r, double* yks, double* chg) {

	double* coulterm = malloc
	for (int j = 0; j < size; j++) {
		coulterm[j] = yks[j] * chg[j] / r[j];
	}
	double** spline = spline_coeff(r, coulterm, size);
	return spline_integral(r, coulterm, spline, size);
}

double* ykint(int k, int size, double* r, double* integrand1, double* integrand2) {
	
	double** spline1 = spline_coeff(r, integrand1, size);
	double** spline2 = spline_coeff(r, integrand2, size);
	double* integrals1 = (double*) malloc(size * sizeof(double));
	double* integrals2 = (double*) malloc(size * sizeof(double));
	integrals1[0] = r[0] * P1[0];
	integrals2[size-1] = 0;
	double* a = integrand1, *b = spline1[0], *c = spline1[1], *d = spline1[2];
	double* e = integrand2, *f = spline2[0], *g = spline2[1], *h = spline2[2];
	double dx=0;
	int j = 0;
	for (int i = 0; i < size-1; i++) {
		dx = r[i+1] - r[i];
		integrals1[i+1] = integrals1[i] + dx * (a[i] + dx * (b[i]/2 + dx * (c[i]/3 + d[i]*dx/4)));
		j = size - i - 2;
		dx = r[j] - r[j-1];
		integrals2[j] = integrals2[j+1] + dx * (e[j] + dx * (f[j]/2 + dx * (g[i]/3 + h[j]*dx/4)));
	}
	for (int i = 0; i < size-1; i++) {
		integrals1[i+1] *= pow(r[i], -k);
		j = size - i - 2;
		integrals2[j] *= pow(r[j], k+1);
	}

	for (int i = 0; i < size; i++) {
		integrals1[i] += integrals2[i];
	}

	free(integrals2);
	free(spline1[0]);
	free(spline1[1]);
	free(spline1[2]);
	free(spline2[0]);
	free(spline2[1]);
	free(spline2[2]);
	return integrals1;

}

