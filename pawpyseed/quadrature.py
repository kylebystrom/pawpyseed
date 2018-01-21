import numpy as np

f = open("quadrature.c", "w")

f.write('#include "quadrature.h"\n\n')

x = []

for i in range(3, 30):
	x.append( np.polynomial.legendre.leggauss(i) )

f.write('double QUADRATURE_POINTS[27][30] = {\n')

for i in range(27):
	print(i)
	currstr = str(x[i][0].tolist()).replace("[", "{").replace("]", "}")
	f.write("\t" + currstr)
	if i != 29: f.write(',\n')

f.write('}\n\n')

f.write('double QUADRATURE_WEIGHTS[27][30] = {\n')

for i in range(27):
	currstr = str(x[i][1].tolist()).replace("[", "{").replace("]", "}")
	f.write("\t" + currstr)
	if i != 29: f.write(',\n')

f.write('}\n')

f.close()