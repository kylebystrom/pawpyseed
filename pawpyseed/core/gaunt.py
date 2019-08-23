"""
\file
Python script for generating Gaunt coefficients and factors used for offsite
partial wave overlap integrals.
"""

from sympy.physics.wigner import gaunt, wigner_3j
import numpy as np
from sympy import N

gcs = np.zeros([4,4,4,7,4])
facs = np.zeros([4,4,4,7,4])
print(gaunt(1,0,1,1,0,0))
print(N(gaunt(1,0,1,1,0,-1)))
print(type(N(gaunt(1,0,1,1,0,-1))))

for l1 in range(4):
	for l2 in range(l1+1):
		for l3 in range(abs(l1-l2), l1+l2+1, 2):
			for m1 in range(-l1,l1+1):
				for m2 in range(0,l2+1):
					val = N(gaunt(l1,l2,l3,m1,m2,-m1-m2))
					gcs[l1][l2][(l3-abs(l1-l2))//2][l1+m1][m2] = val
					val2 = N(wigner_3j(l1,l2,l3,0,0,0)) * N(wigner_3j(l1,l2,l3,-m1,m2,m1-m2))
					val3 = np.sqrt((2*l1+1)*(2*l2+1)*(2*l3+1)/4/np.pi)
					facs[l1][l2][(l3-abs(l1-l2))//2][l1+m1][m2] = val2 * val3
					print(val, val2 * val3)

f = open('gaunt.c', 'w')
f.write('#include "quadrature.h"\n\n')
f.write('double GAUNT_COEFF[%d][%d][%d][%d][%d] = ' % (4,4,4,7,4))
f.write((str(gcs.tolist()).replace('[', '{').replace(']', '}').replace('}, ', '},\n')) + ';\n\n')
f.write('double SBTFACS[%d][%d][%d][%d][%d] = ' % (4,4,4,7,4))
f.write((str(facs.tolist()).replace('[', '{').replace(']', '}').replace('}, ', '},\n')) + ';\n\n')
f.close()

f = open('gaunt.h', 'w')
f.write('#ifndef GAUNT_H\n#define GAUNT_H\n\n')
f.write('extern double GAUNT_COEFF[%d][%d][%d][%d][%d];\n' % (4,4,4,7,4))
f.write('extern double SBTFACS[%d][%d][%d][%d][%d];\n' % (4,4,4,7,4))
f.write('\n#endif\n')
f.close()
