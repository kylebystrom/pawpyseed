import numpy as np

MAXSIZE = 100

f = open("quadrature.c", "w")

f.write('#include "quadrature.h"\n\n')

x = []

for i in range(3, MAXSIZE+1):
	x.append( np.polynomial.legendre.leggauss(i) )

f.write('double QUADRATURE_POINTS[%d][%d] = {\n' % (MAXSIZE-2, MAXSIZE+1))

for i in range(MAXSIZE-2):
	print(i)
	currstr = str(x[i][0].tolist()).replace("[", "{").replace("]", "}")
	f.write("\t" + currstr)
	if i != MAXSIZE: f.write(',\n')

f.write('}\n\n')

f.write('double QUADRATURE_WEIGHTS[%d][%d] = {\n' % (MAXSIZE-2, MAXSIZE+1))

for i in range(MAXSIZE-2):
	currstr = str(x[i][1].tolist()).replace("[", "{").replace("]", "}")
	f.write("\t" + currstr)
	if i != MAXSIZE: f.write(',\n')

f.write('}\n')

f.close()