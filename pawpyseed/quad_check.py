import numpy as np

total = 0

for rstep in range(100):
	numtheta = min( int(rstep*27/100), 26 ) + 3
	print(numtheta)
	pts, wts = np.polynomial.legendre.leggauss(numtheta)
	for thetastep in range(numtheta):
		numphi = int((numtheta*2+2) * (1-pts[thetastep]**2)**0.5)
		for phistep in range(numphi):
			total += 1

print()
print(total)
