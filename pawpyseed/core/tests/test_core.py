# coding: utf-8

import unittest
import os, subprocess, sys
import time
import scipy

import numpy as np
from numpy.testing import assert_almost_equal, assert_equal,\
						assert_raises

from scipy.special import lpmn, sph_harm
from scipy.interpolate import CubicSpline
from scipy.special import factorial2 as fac2
from nose import SkipTest
from nose.tools import nottest
from nose.plugins.skip import Skip

from pawpyseed.core import pawpyc
from pawpyseed.core.tests import testc

class PawpyTestError(Exception):
	"""
	Class for handling errors that occur during execution
	of Python functions in pawpyseed
	"""
	def __init__(self, msg):
		self.msg = msg

MODULE_DIR = os.path.dirname(os.path.abspath(__file__))

from pymatgen.io.vasp.inputs import Poscar, Potcar
from pymatgen.io.vasp.outputs import Vasprun, Chgcar
from pymatgen.core.structure import Structure

from pawpyseed.core.utils import *
from pawpyseed.core.wavefunction import *
from pawpyseed.core.projector import Projector
from pawpyseed.core.noncollinear import NCLWavefunction

class DummyProjector(Projector):

	def make_site_lists(self):
		M_R, M_S, N_R, N_S, N_RS = super(DummyProjector, self).make_site_lists()
		return [], [], M_R, M_S, [pair for pair in zip(M_R, M_S)]

class TestC:

	def setup(self):
		self.currdir = os.getcwd()
		os.chdir(os.path.join(MODULE_DIR, '../../../test_files'))
		self.vr = Vasprun("vasprun.xml")
		self.cr = CoreRegion(Potcar.from_file("POTCAR"))

	def teardown(self):
		os.chdir(self.currdir)

	def test_legendre(self):
		xs = np.linspace(0,1,10000)
		ys = np.zeros(10000*16)
		ys2 = np.zeros(10000*16)
		for i in range(xs.shape[0]):
			if i == 0: print (lpmn(-3,3,xs[i])[0])
			ys[16*i:16*(i+1)] = lpmn(-3,3,xs[i])[0].flatten()
		for i in range(xs.shape[0]):
			for l in range(4):
				for m in range(0,-l-1,-1):
					ys2[i*16-m*4+l] = pawpyc.legendre(l, m, xs[i])
		assert_almost_equal(np.linalg.norm(ys-ys2), 0.0)

	def test_Ylm(self):
		for l in range(4):
			for m in range(-l,l):
				xs = np.linspace(0,1,100)
				ys = np.zeros(10000, np.complex128)
				ys1 = np.zeros(10000, np.complex128)
				ys2 = np.zeros(10000, np.complex128)
				for i in range(xs.shape[0]):
					for j in range(xs.shape[0]):
						ys[i*100+j] = sph_harm(m, l, xs[j]*2*np.pi, xs[i]*np.pi)
				i,j=0,0
				for i in range(xs.shape[0]):
					for j in range(xs.shape[0]):
						ys1[i*100+j] = pawpyc.Ylm(l, m, xs[i]*np.pi, xs[j]*np.pi*2)
						ys2[i*100+j] = pawpyc.Ylm2(l, m, np.cos(xs[i]*np.pi), xs[j]*np.pi*2)
				assert_almost_equal(np.linalg.norm(ys-ys1),0.0)
				assert_almost_equal(np.linalg.norm(ys1-ys2),0.0)

	def test_unit_conversion(self):
		struct = Poscar.from_file("CONTCAR").structure
		lattice = struct.lattice.matrix
		reclattice = struct.lattice.reciprocal_lattice.matrix
		for site in struct:
			coord = site.coords
			fcoord = site.frac_coords
			temp1 = np.copy(fcoord, order='C')
			pawpyc.frac_to_cartesian(temp1, lattice)
			assert_almost_equal(np.linalg.norm(temp1-coord), 0.0)
			temp2 = np.copy(coord, order='C')
			pawpyc.cartesian_to_frac(temp2, reclattice)
			assert_almost_equal(np.linalg.norm(temp2-fcoord), 0.0)

	def test_spline(self):
		vr = self.vr
		cr = self.cr
		struct = Poscar.from_file("CONTCAR").structure
		grid = cr.pps['Ga'].projgrid
		vals = cr.pps['Ga'].realprojs[0]
		rmax = cr.pps['Ga'].rmax
		tst = np.linspace(0, max(grid), 400)
		res1 = CubicSpline(grid, vals,
			extrapolate=True, bc_type='natural')(tst)
		x, y = grid[:], vals[:]
		res2 = np.zeros(tst.shape)
		pawpyc.interpolate(res2, tst, x, y, rmax, 100, 400)
		assert_almost_equal(np.linalg.norm(res1-res2), 0, 2)
		sys.stdout.flush()

	def test_fft3d(self):
		print("TEST FFT")
		vr = self.vr 
		weights = np.array(vr.actual_kpoints_weights)
		res = testc.fft_check("WAVECAR", weights, np.array([20,20,20], dtype=np.int32, order='C'))
		assert_equal(res, 0)

	def test_sbt(self):
		from scipy.special import spherical_jn as jn
		cr = CoreRegion(Potcar.from_file("POTCAR"))
		r = cr.pps['Ga'].grid
		f = cr.pps['Ga'].aewaves[0] - cr.pps['Ga'].pswaves[0];
		ks, res = pawpyc.spherical_bessel_transform(1e6, 0, r, f)
		k = ks[180]
		vals = jn(0, r * k) * f * r
		integral = np.trapz(vals, r)
		print (integral)
		print (ks[180])
		print (res[180])
		assert_almost_equal(integral, res[180], decimal=3)
		perc = ks**2 * res**2
		perc = np.cumsum(perc)
		perc /= np.max(perc)

	def test_radial(self):
		try:
			import pawpyseed.core.tests.reference as gint
		except ImportError:
			print("No McMurchie-Davidson installed, skipping radial test")
			raise SkipTest()

		h0 = ([1], [0])
		h1 = ([2], [1])
		h2 = ([4, -2], [2, 0])
		h3 = ([8, -12], [3, 1])
		H = [h0, h1, h2, h3]

		def eval_overlap(a, ijk1, A, b, ijk2, B):
			ov = 0
			for ci1, i1 in zip(*(H[ijk1[0]])):
				for cj1, j1 in zip(*(H[ijk1[1]])):
					for ck1, k1 in zip(*(H[ijk1[2]])):
						for ci2, i2 in zip(*(H[ijk2[0]])):
							for cj2, j2 in zip(*(H[ijk2[1]])):
								for ck2, k2 in zip(*(H[ijk2[2]])):
									ov += ci1 * ci2 * cj1 * cj2 * ck1 *ck2\
											* gint.overlap(a, (i1,j1,k1), A,
															b,(i2,j2,k2), B)
			return ov

		def getf(f, r, a, n, m):
			if n == 0:
				N = 0.25
				ijks = [(0,0,0)]
				coefs = [1]
			elif n == 1 and m == 0:
				N = .25
				ijks = [(0,0,1)]
				coefs = [1]
			elif n == 1:
				N = 0.25
				ijks = [(1,0,0), (0,1,0)]
				if m == 1:
					coefs = [1,1j]
				else:
					coefs = [1,-1j]
			elif n == 2 and m == 0:
				N = 3
				ijks = [(0,0,2), (2,0,0), (0,2,0)]
				coefs = [2,-1,-1]
			elif n == 2 and m == 1:
				N = 0.25
				ijks = [(1,0,1), (0,1,1)]
				coefs = [1, 1.0j]
			elif n == 2 and m == -1:
				N = 0.25
				ijks = [(1,0,1), (0,1,1)]
				coefs = [1, -1.0j]
			elif n ==2 and m == 2:
				N = 1
				ijks = [(2,0,0), (0,2,0), (1,1,0)]
				coefs=[1, -1, 1.0j*2]
			elif n ==2 and m == -2:
				N = 1
				ijks = [(2,0,0), (0,2,0), (1,1,0)]
				coefs=[1, -1, -1.0j*2]
			elif n == 3 and m == 0:
				N = 15
				ijks = [(0,0,3), (0,2,1), (2,0,1)]
				coefs = [2, -3, -3]
			else:
				raise ValueError('Do not know that n,m pair %d %d' % (n, m))
			return r * 2**(n+2) * np.sqrt(np.pi * N / fac2(2*n+1)) * (a*r)**n * f, ijks, coefs


			"""
			elif n == 2 and m == 1:
				N = 0.25
				ijks = [(1,0,1)]
				coefs = [1]
			elif n == 2 and m == -1:
				N = 0.25
				ijks = [(0,1,1)]
				coefs = [1]
			"""
			"""
			elif n == 2 and m == 2:
				N = 1
				ijks = [(2,0,0), (0,2,0)]
				coefs = [1,-1]
			elif n == 2 and m == -2:
				N = 0.25
				ijks = [(1,1,0)]
				coefs = [1]
			"""

		# test realspace radial overlap
		# test recipspace radial overlap
		pols = ([0,0,0],[0,0,1],[0,0,2],[0,0,3],[0,1,1]) 
		ls = (0,1,2,3,2, 2, 2,2,1, 1)
		ms = (0,0,0,0,1,-1,-2,2,1,-1)
		a = 1
		b = 1
		A = [0,0,0]
		Bs = ([0,0,0], [0.5,0.5,0.5], [-0.5,0.5,-0.5], [0.123,0.543,-0.96])
		r = np.exp(np.linspace(np.log(0.001),np.log(4), 600))
		init_f1 = np.exp(-a * r**2)
		init_f2 = np.exp(-b * r**2)
		for B in Bs:
			Barr = np.array(B, dtype=np.float64, order='C')
			for l1, m1 in zip(ls, ms):
				f1, ijks1, coefs1 = getf(init_f1, r, a, l1, m1)
				for l2, m2 in zip(ls, ms):
					#print("START", l1, m1, l2, m2, a, b, B)
					f2, ijks2, coefs2 = getf(init_f2, r, b, l2, m2)
					ov1 = 0
					for coef1, ijk1 in zip(coefs1, ijks1):
						for coef2, ijk2 in zip(coefs2, ijks2):
							#print(ijk1, ijk2, coef1, coef2)
							ov1 += coef1 * np.conj(coef2) * eval_overlap(a,ijk1,A,b,ijk2,B)
					if m1 != 0:
						ov1 /= np.sqrt(2)
					if m2 != 0:
						ov1 /= np.sqrt(2)
					if np.linalg.norm(B) > 0:
						# sign convention adjustment
						if m1 > 0:
							ov1 *= (-1)**(m1)
						if m2 > 0:
							ov1 *= (-1)**(m2)

					ov2 = pawpyc.reciprocal_offsite_wave_overlap(Barr,
						r, f1, r, f2,
						l1, m1, l2, m2)
					if (np.abs(ov1) > 1e-10 and (np.abs(ov1 - ov2))/np.abs(ov1) > 1e-4):
						print(l1, m1, l2, m2, B, ov1, ov2)
					if np.abs(ov1) < 1e-10:
						assert_almost_equal(np.abs(ov2), 0, 10)
					else:
						assert_almost_equal((np.abs(ov1 - ov2))/np.abs(ov1), 0, 4)

	@nottest
	def test_sphagain(self):
		# a little snippet to check the sign of teh gaussians vs spherical harmonics
		h0 = ([1], [0])
		h1 = ([2], [1])
		h2 = ([4, -2], [2, 0])
		h3 = ([8, -12], [3, 1])
		H = [h0, h1, h2, h3]
		for n, m in [(2,1), (2,-1), (2,0), (2,2), (2,-2)]:
			if n == 0:
				N = 0.25
				ijks = [(0,0,0)]
				coefs = [1]
			elif n == 1 and m == 0:
				N = .25
				ijks = [(0,0,1)]
				coefs = [1]
			elif n == 1:
				N = 0.25
				ijks = [(1,0,0), (0,1,0)]
				if m == 1:
					coefs = [1,1j]
				else:
					coefs = [1,-1j]
			elif n == 2 and m == 0:
				N = 3
				ijks = [(0,0,2), (2,0,0), (0,2,0)]
				coefs = [2,-1,-1]
			elif n == 2 and m == 1:
				N = 0.25
				ijks = [(1,0,1), (0,1,1)]
				coefs = [1, 1.0j]
			elif n == 2 and m == -1:
				N = 0.25
				ijks = [(1,0,1), (0,1,1)]
				coefs = [1, -1.0j]
			elif n ==2 and m == 2:
				N = 1
				ijks = [(2,0,0), (0,2,0), (1,1,0)]
				coefs=[1, -1, 1.0j*2]
			elif n ==2 and m == -2:
				N = 1
				ijks = [(2,0,0), (0,2,0), (1,1,0)]
				coefs=[1, -1, -1.0j*2]
			elif n == 3 and m == 0:
				N = 15
				ijks = [(0,0,3), (0,2,1), (2,0,1)]
				coefs = [2, -3, -3]
			ov = 0
			o, p = np.pi/4, np.pi/4
			for ijk1, coef in zip(ijks, coefs):
				for ci1, i1 in zip(*(H[ijk1[0]])):
					for cj1, j1 in zip(*(H[ijk1[1]])):
						for ck1, k1 in zip(*(H[ijk1[2]])):
							ov += ci1 * cj1 * ck1 * coef\
									* (np.sin(o)*np.cos(p))**i1\
									* (np.sin(o)*np.sin(p))**j1\
									* (np.cos(o))**k1
			print("SPHHARM CHECK!!!")
			print(sph_harm(m, n, o, p))
			print(ov)


@nottest
class TestMem:

	def setup(self):
		self.currdir = os.getcwd()
		os.chdir(os.path.join(MODULE_DIR, '../../../test_files'))
		structure = Poscar.from_file("CONTCAR").structure
		cr = CoreRegion(Potcar.from_file("POTCAR"))
		pps = {}
		labels = {}
		label = 0
		for e in cr.pps:
			pps[label] = cr.pps[e]
			labels[e] = label
			label += 1
		clabels = np.array([], np.int32)
		ls = np.array([], np.int32)
		projectors = np.array([], np.float64)
		aewaves = np.array([], np.float64)
		pswaves = np.array([], np.float64)
		wgrids = np.array([], np.float64)
		pgrids = np.array([], np.float64)
		num_els = 0

		for num in pps:
			pp = pps[num]
			clabels = np.append(clabels, [num, len(pp.ls), pp.ndata, len(pp.grid)])
			ls = np.append(ls, pp.ls)
			wgrids = np.append(wgrids, pp.grid)
			pgrids = np.append(pgrids, pp.projgrid)
			num_els += 1
			for i in range(len(pp.ls)):
				proj = pp.realprojs[i]
				aepw = pp.aewaves[i]
				pspw = pp.pswaves[i]
				projectors = np.append(projectors, proj)
				aewaves = np.append(aewaves, aepw)
				pswaves = np.append(pswaves, pspw)
		print ("rmax", cr.pps['Ga'].rmax * 0.529177)

		selfnums = np.array([labels[el(s)] for s in structure], dtype=np.int32)
		selfcoords = np.array([], np.float64)

		for s in structure:
			selfcoords = np.append(selfcoords, s.frac_coords)

		f = open('potholder.txt', 'w')
		f.write('%d '%num_els)
		for i in range(len(clabels)):
			f.write('%d '%clabels[i])
		f.write('%d '%len(ls))
		for i in range(len(ls)):
			f.write('%d '%ls[i])
		f.write('%d '%len(pgrids))
		for i in range(len(pgrids)):
			f.write('%f '%pgrids[i])
		f.write('%d '%len(wgrids))
		for i in range(len(wgrids)):
			f.write('%f '%wgrids[i])
		f.write('%d '%len(projectors))
		for i in range(len(projectors)):
			f.write('%f '%projectors[i])
		f.write('%d '%len(aewaves))
		for i in range(len(aewaves)):
			f.write('%f '%aewaves[i])
		f.write('%d '%len(pswaves))
		for i in range(len(pswaves)):
			f.write('%f '%pswaves[i])
		f.write('%f '%(cr.pps['Ga'].rmax*0.529177))
		f.write('%d '%len(selfnums))
		for i in range(len(selfnums)):
			f.write('%d '%selfnums[i])
		for i in range(len(selfnums)*3):
			f.write('%f '%selfcoords[i])

		f.close()

	def teardown(self):
		os.remove('potholder.txt')
		os.chdir(self.currdir)

	def test_memory(self):
		f = open('mtest.out', 'w')
		subprocess.call('valgrind ../pawpyseed/core/memtest --track-origins=yes --leak-check=full'.split(),
			stdout=f, stderr=f)
		f.close()


class TestPy:

	def setup(self):
		self.currdir = os.getcwd()
		os.chdir(os.path.join(MODULE_DIR, '../../../test_files'))

	def teardown(self):
		os.chdir(self.currdir)

	def test_init(self):
		print("TEST INIT")
		sys.stdout.flush()
		wf = Wavefunction.from_directory('.', True)
		wf = Wavefunction.from_directory('.', False)
		wf = Wavefunction.from_files('CONTCAR', 'WAVECAR',
			'POTCAR', 'vasprun.xml', True)
		wf = Wavefunction.from_files('CONTCAR', 'WAVECAR',
			'POTCAR', 'vasprun.xml', False)
		wf = Wavefunction.from_files('CONTCAR', 'WAVECAR2.gz',
			'POTCAR', 'vasprun.xml', False)
		with assert_raises(FileNotFoundError):
			wf = Wavefunction.from_files('bubbles', 'WAVECAR',
				'POTCAR', 'vasprun.xml', True)

	def test_writestate(self):
		print("TEST WRITE")
		sys.stdout.flush()
		wf = Wavefunction.from_directory('.')
		fileprefix = ''
		b, k, s = 10, 1, 0
		state1 = wf.write_state_realspace(b, k, s, fileprefix = "", 
			dim=np.array([30,30,30]))
		wf = Wavefunction.from_directory('.', False)
		state2 = wf.write_state_realspace(b, k, s, fileprefix = "", 
			dim=np.array([30,30,30]))
		assert_almost_equal(np.linalg.norm(state1-state2),0)
		state2 = wf.write_state_realspace(b, k, s, fileprefix = "", 
			dim=np.array([30,30,30]), remove_phase=True)
		assert state1.shape[0] == 30
		assert state1.shape[1] == 30
		assert state1.shape[2] == 30
		assert state1.dtype == np.complex128
		filename_base = "%sB%dK%dS%d" % (fileprefix, b, k, s)
		filename1 = "%s_REAL" % filename_base
		filename2 = "%s_IMAG" % filename_base
		#os.remove(filename1)
		#os.remove(filename2)
		with assert_raises(ValueError):
			wf.write_state_realspace(-1, 0, 0)

	def test_density(self):
		print("TEST DENSITY")
		sys.stdout.flush()
		wf = Wavefunction.from_directory('nosym')
		#wf = wf.desymmetrized_copy()
		wf.write_density_realspace(dim=np.array([40,40,40]), scale = wf.structure.lattice.volume)
		tstchg = Chgcar.from_file("AECCAR2").data['total']# / wf.structure.volume
		chg = Chgcar.from_file("PYAECCAR").data['total']
		reldiff = np.sqrt(np.mean(((chg-tstchg)/tstchg)**2))
		newchg = chg-tstchg
		Chgcar(Poscar(wf.structure), {'total': newchg}).write_file('DIFFCHGCAR.vasp')
		print(np.sum(chg)/40**3, np.sum(tstchg)/40**3)
		assert_almost_equal(reldiff, 0, decimal=2)
		wf = Wavefunction.from_directory('nosym')
		res = wf.write_density_realspace(filename="BAND4DENS", bands=4)
		#os.remove('PYAECCAR')
		print("DENS shape", res.shape)
		assert_almost_equal(np.sum(res)*wf.structure.lattice.volume/np.cumprod(res.shape)[-1], 1, 4)

	def test_state_wf_and_density(self):
		wf = Wavefunction.from_directory('.')
		chg_from_wf = np.abs(wf.get_state_realspace(0, 0, 0))**2
		chg = wf.get_state_realspace_density(0,0,0,dim=wf.dim)
		assert_equal(chg.shape, chg_from_wf.shape)
		dv = wf.structure.volume / np.cumprod(chg.shape)[-1]
		assert_almost_equal(np.sum(chg)*dv, 1, 3)
		assert_almost_equal(np.sum(chg_from_wf)*dv, 1, 3)
		reldiff = np.sqrt(np.mean(np.abs(chg-chg_from_wf)))
		assert_almost_equal(reldiff, 0, decimal=2)

	def test_pseudoprojector(self):
		print("TEST PSEUDO")
		sys.stdout.flush()
		# test ps projections
		wf = Wavefunction.from_directory('.')
		basis = Wavefunction.from_directory('.')
		pr = Projector(wf, basis, method = "pseudo")
		res = pr.single_band_projection(6)
		assert res.shape[0] == basis.nband * basis.nspin * basis.nwk
		res = pr.defect_band_analysis(4, 10, False)
		assert len(res.keys()) == 15

	def test_projector(self):
		print("TEST PROJ")
		sys.stdout.flush()
		# test ae projections
		wf1 = Wavefunction.from_directory('.', False)
		basis = Wavefunction.from_directory('.', False)
		pr = Projector(wf1, basis)
		for b in range(wf1.nband):
			v, c = pr.proportion_conduction(b)
			if b < 6:
				assert_almost_equal(v, 1, decimal=4)
				assert_almost_equal(c, 0, decimal=8)
			else:
				assert_almost_equal(v, 0, decimal=8)
				assert_almost_equal(c, 1, decimal=4)

		generator = Projector.setup_multiple_projections('.', ['.', '.'])
		for wf_dir, wf in generator:
			wf.defect_band_analysis(4, 10, spinpol=True)

		wf1 = Wavefunction.from_directory('.', False)
		basis = Wavefunction.from_directory('.', False)
		pr = Projector(wf1, basis, 'aug_recip')
		for b in range(wf1.nband):
			v, c = pr.proportion_conduction(b)
			if b < 6:
				assert_almost_equal(v, 1, decimal=4)
				assert_almost_equal(c, 0, decimal=8)
			else:
				assert_almost_equal(v, 0, decimal=8)
				assert_almost_equal(c, 1, decimal=4)

		wf1 = Wavefunction.from_directory('.', False)
		basis = Wavefunction.from_directory('.', False)
		pr = Projector(wf1, basis, 'realspace')
		for b in range(wf1.nband):
			v, c = pr.proportion_conduction(b)
			if b < 6:
				assert_almost_equal(v, 1, decimal=4)
				assert_almost_equal(c, 0, decimal=8)
			else:
				assert_almost_equal(v, 0, decimal=8)
				assert_almost_equal(c, 1, decimal=4)

		with assert_raises(ValueError):
			pr.proportion_conduction(100)
			pr.proportion_conduction(-1)

	def test_projector_gz(self):
		print("TEST PROJGZ")
		sys.stdout.flush()
		# test ae projections
		wf1 = Wavefunction.from_files('CONTCAR', 'WAVECAR2.gz',
			'POTCAR', 'vasprun.xml', False)
		basis = Wavefunction.from_directory('.', False)
		pr = Projector(wf1, basis)
		for b in range(wf1.nband):
			v, c = pr.proportion_conduction(b)
			if b < 6:
				assert_almost_equal(v, 1, decimal=4)
				assert_almost_equal(c, 0, decimal=8)
			else:
				assert_almost_equal(v, 0, decimal=8)
				assert_almost_equal(c, 1, decimal=4)

		generator = Projector.setup_multiple_projections('.', ['.', '.'])
		for wf_dir, wf in generator:
			wf.defect_band_analysis(4, 10, spinpol=True)

	def test_offsite(self):
		Projector = DummyProjector
		print("TEST OFFSITE")
		sys.stdout.flush()
		wf1 = Wavefunction.from_directory('nosym', False)
		basis = Wavefunction.from_directory('nosym', False)
		print("ENCUT", wf1.encut, basis.encut)
		pr = Projector(wf1, basis)
		test_vals = {}
		for b in range(wf1.nband):
			v, c = pr.proportion_conduction(b)
			print("CHECK_VALS", v,c)
			test_vals[b] = (v,c)
		for b in range(wf1.nband//2):
			if b < 6:
				assert_almost_equal(test_vals[b][0], 1, decimal=2)
				assert_almost_equal(test_vals[b][1], 0, decimal=4)
			else:
				assert_almost_equal(test_vals[b][0], 0, decimal=4)
				assert_almost_equal(test_vals[b][1], 1, decimal=2)

		generator = Projector.setup_multiple_projections('.', ['.', '.'])
		for wf_dir, wf in generator:
			wf.defect_band_analysis(4, 10, spinpol=True)

		wf1 = Wavefunction.from_directory('nosym', False)
		basis = Wavefunction.from_directory('nosym', False)
		print("ENCUT", wf1.encut, basis.encut)
		pr = Projector(wf1, basis, 'aug_recip')
		test_vals = {}
		for b in range(wf1.nband):
			v, c = pr.proportion_conduction(b)
			print("CHECK_VALS", v,c)
			test_vals[b] = (v,c)
		for b in range(wf1.nband//2):
			if b < 6:
				assert_almost_equal(test_vals[b][0], 1, decimal=2)
				assert_almost_equal(test_vals[b][1], 0, decimal=4)
			else:
				assert_almost_equal(test_vals[b][0], 0, decimal=4)
				assert_almost_equal(test_vals[b][1], 1, decimal=2)

	def test_desymmetrization(self):
		print("TEST DESYM")
		sys.stdout.flush()
		# test ae projections
		wf1 = Wavefunction.from_directory('.', False)
		basis = Wavefunction.from_directory('nosym', False)
		pr = Projector(wf1, basis, unsym_wf=True, unsym_basis=True)
		test_vals = {}
		for b in range(wf1.nband):
			v, c = pr.proportion_conduction(b)
			print("CHECK_VALS", v,c)
			test_vals[b] = (v,c)
		for b in range(wf1.nband//2):
			if b < 6:
				assert_almost_equal(test_vals[b][0], 1, decimal=3)
				assert_almost_equal(test_vals[b][1], 0, decimal=7)
			else:
				assert_almost_equal(test_vals[b][0], 0, decimal=7)
				assert_almost_equal(test_vals[b][1], 1, decimal=3)

	def test_norm(self):
		wf = Wavefunction.from_directory('.', setup_projectors=True)
		wf.check_c_projectors()
		testc.proj_check(wf)

	def test_writestate_ncl(self):
		print("TEST WRITE NCL")
		sys.stdout.flush()
		wf = NCLWavefunction.from_directory('noncollinear')
		fileprefix = ''
		b, k, s = 10, 1, 0
		state1, state2 = wf.write_state_realspace(b, k, s, fileprefix = "", 
			dim=np.array([30,30,30]))
		assert state1.shape[0] == 30
		assert state1.shape[1] == 30
		assert state1.shape[2] == 30
		assert state1.dtype == np.complex128
		filename_base = "%sB%dK%dS%d" % (fileprefix, b, k, s)
		filename1 = "%s_UP_REAL" % filename_base
		filename2 = "%s_UP_IMAG" % filename_base
		filename3 = "%s_DOWN_REAL" % filename_base
		filename4 = "%s_DOWN_IMAG" % filename_base
		os.remove(filename1)
		os.remove(filename2)
		os.remove(filename3)
		os.remove(filename4)

	def test_density_ncl(self):
		print("TEST DENSITY NCL")
		sys.stdout.flush()
		print("LOAD WAVEFUNCTION")
		sys.stdout.flush()
		wf = NCLWavefunction.from_directory('noncollinear')
		print("FINISHED LOAD WAVEFUNCTION")
		sys.stdout.flush()
		#wf = wf.desymmetrized_copy()
		wf.write_density_realspace(dim=np.array([40,40,40]), scale = wf.structure.lattice.volume)
		tstchg = Chgcar.from_file("AECCAR2").data['total']# / wf.structure.volume
		chg = Chgcar.from_file("PYAECCAR").data['total']
		reldiff = np.sqrt(np.mean(((chg-tstchg)/tstchg)**2))
		newchg = chg-tstchg
		Chgcar(Poscar(wf.structure), {'total': newchg}).write_file('DIFFCHGCAR.vasp')
		print(np.sum(chg)/40**3, np.sum(tstchg)/40**3)
		assert_almost_equal(reldiff, 0, decimal=3)

	def test_flip_spin(self):

		wf = Wavefunction.from_directory('.', False)
		basis = Wavefunction.from_directory('.', False)
		pr = Projector(wf, basis)
		for b in range(wf.nband):
			res = pr.single_band_projection(b, flip_spin=True)
			res = np.abs(res)**2
			print(b,res)
			for br in range(basis.nband):
				expected = 1 if b == br else 0
				decimal = 4 if b == br else 8
				for k in range(basis.nspin * basis.nwk):
					assert_almost_equal(np.sum(res[br*basis.nwk*basis.nspin + k]),
										expected, decimal=decimal)