#ifndef HF_UTILS_H
#define HF_UTILS_H

double legendre(int l, int m, double x);

double Ylm(int l, int m, double theta, double phi);

double Ylmd(int l, int m, double costheta, double phi);

double laguerre(double x, int n, int alpha);

double hradial(int n, int l, double r);

double* yk(int k, int size, double* r, double* P1, double* P2);

#endif

