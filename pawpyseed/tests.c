#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include "tests.h"
#include "utils.h"

double Ylmr(int l, int m, double theta, double phi) { return creal(Ylm(l, m, theta, phi)); }

double Ylmi(int l, int m, double theta, double phi) { return cimag(Ylm(l, m, theta, phi)); }