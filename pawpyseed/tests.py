import unittest
import os
import time

import numpy as np
from numpy.testing import assert_almost_equal

from scipy.special import lpmn, sph_harm

from ctypes import *

MODULE_DIR = os.path.dirname(os.path.abspath(__file__))
PAWC = CDLL(os.path.join(MODULE_DIR, "pawpy.so"))
PAWC.legendre.restype = c_double
PAWC.Ylmr.restype = c_double
PAWC.Ylmi.restype = c_double

class TestC:

	def setup(self):
		pass

	def teardown(self):
		pass

	def test_fac(self):
		assert PAWC.fac(5) == 120

	def test_legendre(self):
		xs = np.linspace(0,1,10000)
		ys = np.zeros(10000*16)
		ys2 = np.zeros(10000*16)
		t1 = time.time()
		for i in range(xs.shape[0]):
			ys[16*i:16*(i+1)] = lpmn(3,3,xs[i])[0].flatten()
		t2 = time.time()
		for i in range(xs.shape[0]):
			for l in range(4):
				for m in range(l+1):
					ys2[i*16+m*4+l] = PAWC.legendre(l, m, c_double(xs[i]))
		t3 = time.time()
		#rough speed check
		assert 2*(t2-t1) > t3 - t2
		#accuracy check
		assert_almost_equal(np.linalg.norm(ys-ys2), 0.0)

	def test_Ylm(self):
		for l in range(4):
			for m in range(0,l+1):
				xs = np.linspace(0,1,100)
				ys = np.zeros(10000, np.complex128)
				ys2 = np.zeros(10000, np.complex128)
				t1 = time.time()
				for i in range(xs.shape[0]):
					for j in range(xs.shape[0]):
						ys[i*100+j] = (-1)**m * sph_harm(m, l, xs[j]*2*np.pi, xs[i]*np.pi)
				t2 = time.time()
				i,j=0,0
				for i in range(xs.shape[0]):
					for j in range(xs.shape[0]):
						ys2[i*100+j] = PAWC.Ylmr(l, m, c_double(xs[i]*np.pi), c_double(xs[j]*np.pi*2))\
										+ 1.0j*PAWC.Ylmi(l, m, c_double(xs[i]*np.pi), c_double(xs[j]*np.pi*2))
				t3 = time.time()
				#rough speed check
				assert 2*(t2-t1) > t3 - t2
				assert_almost_equal(np.linalg.norm(ys-ys2),0.0)

	


class TestPy:

	def setup():
		pass

	def teardown():
		pass

t = TestC()
t.test_legendre()
t.test_Ylm()