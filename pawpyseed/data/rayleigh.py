import numpy as np 
from scipy.special import spherical_jn, sph_harm

k = np.array([0.6,0.2,0.3])*2*np.pi

def planewave(coord):
	return np.exp(1j * (k[0] * grid[0] + k[1] * grid[1] + k[2] * grid[2])) \
		* np.exp(1j * np.dot(k, [1,1,1]))

m, l, = 1, 1

def func(coord):
	r = np.sqrt(coord[0]**2+coord[1]**2+coord[2]**2)
	theta = np.nan_to_num(np.arccos(coord[2]/r))
	phi = np.nan_to_num(np.arctan2(coord[1]/r, coord[0]/r))
	return np.exp(-r) * sph_harm(m, l, phi, theta)

def rad_int(r):
	kmag = np.linalg.norm(k)
	theta = np.nan_to_num(np.arccos(k[2]/kmag))
	phi = np.nan_to_num(np.arctan2(k[1]/kmag, k[0]/kmag))
	vals = spherical_jn(l, kmag*r) * r * r * np.exp(-r)
	print (vals)
	return np.trapz(vals, r) * np.conj(sph_harm(m, l, phi, theta)) * 4 * np.pi * (1j)**l

x = np.arange(-8, 8, 0.1)
y = np.arange(-8, 8, 0.1)
z = np.arange(-8, 8, 0.1)
x, y, z = np.meshgrid(x, y, z, indexing='ij')
grid = np.array([x, y, z])
w = planewave(grid)
f = func(grid)
integrand = w * np.conj(f)
temp1 = np.trapz(integrand, x, axis = 0)
temp2 = np.trapz(temp1, y[0,:,:], axis = 0)
final = np.trapz(temp2, z[0,0,:], axis = 0)
print(final)
r = np.arange(0,8,0.1)
integral = rad_int(r)
print(integral * np.exp(1j * np.dot(k, [1,1,1])))

