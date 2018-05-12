/**
\file
Stores gaunt coefficients and offsite overlap integral factors. Access as
l1, l2, (L-|l1-l2|)/2, l1+m1, m2. Constrained by l1 >= l2, m2 >= 0,
L >= |l1-l2|, L <= l1+l2
*/

#ifndef GAUNT_H
#define GAUNT_H

extern double GAUNT_COEFF[4][4][4][7][4];
extern double SBTFACS[4][4][4][7][4];

#endif
