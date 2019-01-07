# coding: utf-8

import unittest
import os, subprocess, sys
import time
import scipy

import numpy as np
from numpy.testing import assert_almost_equal

from scipy.special import lpmn, sph_harm
from nose import SkipTest
from nose.tools import nottest

COMPILE = True

class PawpyTestError(Exception):
	"""
	Class for handling errors that occur during execution
	of Python functions in pawpyseed
	"""
	def __init__(self, msg):
		self.msg = msg

if COMPILE:
	currdir = os.getcwd()
	MODULE_DIR = os.path.dirname(os.path.abspath(__file__))
	os.chdir(MODULE_DIR+'/..')
	if not "PAWPYCC" in os.environ:
		if subprocess.call("which icc".split()) == 0:
			os.environ["PAWPYCC"] = "icc"
		elif subprocess.call("which gcc".split()) == 0:
			os.environ["PAWPYCC"] = "gcc"
		else:
			raise PawpyTestError("Can't find icc or gcc compiler!")

	status = subprocess.call('make testsuite'.split())
	if status != 0:
		raise PawpyTestError("Can't compile test pawpy.so! Check the C error output for details.")
	#status = subprocess.call('make mem'.split())
	#if status != 0:
	#	raise PawpyTestError("Can't compile memtest! Check the C error output for details.")
	os.chdir(currdir)

from pymatgen.io.vasp.inputs import Poscar, Potcar
from pymatgen.io.vasp.outputs import Vasprun, Chgcar
from pymatgen.core.structure import Structure

from pawpyseed.core.utils import *
from pawpyseed.core.wavefunction import *
from pawpyseed.core.projector import Projector

from ctypes import *

MODULE_DIR = os.path.dirname(os.path.abspath(__file__))
#PAWC = CDLL(os.path.join(MODULE_DIR, "../pawpy.so"))
PAWC.fac.restype = c_double
PAWC.legendre.restype = c_double
PAWC.Ylmr.restype = c_double
PAWC.Ylmi.restype = c_double
PAWC.spline_coeff.restype = POINTER(POINTER(c_double))
PAWC.proj_interpolate.restype = c_double
PAWC.spherical_bessel_transform_setup.restype = POINTER(None)
PAWC.wave_spherical_bessel_transform.argtypes = [c_void_p, POINTER(c_double), c_int]
PAWC.wave_spherical_bessel_transform.restype = POINTER(c_double)
PAWC.free_sbt_descriptor.argtypes = [c_void_p]
PAWC.free_sbt_descriptor.restype = None
PAWC.fft_check.argtypes = [c_char_p, c_double_p, c_int_p]

def cdouble_to_numpy(arr, length):
	arr = cast(arr, POINTER(c_double))
	newarr = np.zeros(length)
	for i in range(length):
		newarr[i] = arr[i]
	#PAWC.free_ptr(arr)
	return newarr

def cfloat_to_numpy(arr, length):
	arr = cast(arr, POINTER(c_float))
	newarr = np.zeros(length)
	for i in range(length):
		newarr[i] = arr[i]
	#PAWC.free_ptr(arr)
	return newarr

def cfloat_to_numpy(arr, length):
	arr = cast(arr, POINTER(c_int))
	newarr = np.zeros(length)
	for i in range(length):
		newarr[i] = arr[i]
	#PAWC.free_ptr(arr)
	return newarr

def numpy_to_cdouble(arr):
	newarr = (c_double * len(arr))()
	for i in range(len(arr)):
		newarr[i] = arr[i]
	return newarr

def numpy_to_cfloat(arr):
	newarr = (c_float * len(arr))()
	for i in range(len(arr)):
		newarr[i] = arr[i]
	return newarr

def numpy_to_cint(arr):
	newarr = (c_int * len(arr))()
	for i in range(len(arr)):
		newarr[i] = int(arr[i])
	return newarr


class DummyProjector(Projector):

	def setup_overlap(self):
		"""
		DEBUG VERSION OF THE setup_overlap METHOD IN THE
		Projector CLASS
		"""
		basis = self.basis
		self.check_c_projectors()
		basis.check_c_projectors()
		projector_list = self.projector_list
		basisnums = basis.nums
		basiscoords = basis.coords
		selfnums = self.nums
		selfcoords = self.coords
		
		M_R, M_S, N_R, N_S, N_RS = self.make_site_lists()
		num_N_RS = len(N_RS)
		if num_N_RS > 0:
			N_RS_R, N_RS_S = zip(*N_RS)
		else:
			N_RS_R, N_RS_S = [], []
		self.site_cat = [M_R, M_S, N_R, N_S, N_RS_R, N_RS_S]
		cfunc_call(PAWC.overlap_setup_real, None, basis.pwf.wf_ptr, self.pwf.wf_ptr,
					basisnums, selfnums, basiscoords, selfcoords,
					M_R, M_S, M_R, M_S, len(M_R), len(M_R), len(M_R))

	def single_band_projection(self, band_num):
		"""
		DEBUG VERSION OF THE single_band_projection METHOD IN THE
		Projector CLASS
		"""

		if self.pseudo:
			return self.pwf.pseudoprojection(band_num, basis.pwf)

		basis = self.basis
		nband = basis.nband
		nwk = basis.nwk
		nspin = basis.nspin
		res = cfunc_call(PAWC.pseudoprojection, 2*nband*nwk*nspin,
						basis.pwf.wf_ptr, self.pwf.wf_ptr, band_num)
		print("datsa", nband, nwk, nspin)
		sys.stdout.flush()
		projector_list = self.projector_list
		basisnums = basis.nums
		basiscoords = basis.coords
		selfnums = self.nums
		selfcoords = self.coords

		M_R, M_S, N_R, N_S, N_RS_R, N_RS_S = self.site_cat
		
		ct = cfunc_call(PAWC.compensation_terms, 2*nband*nwk*nspin,
						band_num, self.pwf.wf_ptr, basis.pwf.wf_ptr,
						len(self.cr.pps),
						0, len(M_R), len(M_S), len(M_S),
						np.array([]), np.array([]), M_R, M_S, M_R, M_S,
						selfnums, selfcoords, basisnums, basiscoords,
						self.dim)
		res += ct
		return res[::2] + 1j * res[1::2]

class TestC:

	def setup(self):
		self.currdir = os.getcwd()
		os.chdir(os.path.join(MODULE_DIR, '../../../test_files'))
		self.vr = Vasprun("vasprun.xml")
		self.cr = CoreRegion(Potcar.from_file("POTCAR"))

	def teardown(self):
		os.chdir(self.currdir)

	def test_fac(self):
		print(int(PAWC.fac(5)))
		assert int(PAWC.fac(5)) == 120

	def test_legendre(self):
		xs = np.linspace(0,1,10000)
		ys = np.zeros(10000*16)
		ys2 = np.zeros(10000*16)
		t1 = time.time()
		for i in range(xs.shape[0]):
			if i == 0: print (lpmn(-3,3,xs[i])[0])
			ys[16*i:16*(i+1)] = lpmn(-3,3,xs[i])[0].flatten()
		t2 = time.time()
		for i in range(xs.shape[0]):
			for l in range(4):
				for m in range(0,-l-1,-1):
					if i == 0: print(PAWC.legendre(l, m, c_double(xs[i])))
					ys2[i*16-m*4+l] = PAWC.legendre(l, m, c_double(xs[i]))
		t3 = time.time()
		#rough speed check
		#assert 2*(t2-t1) > t3 - t2
		#accuracy check
		assert_almost_equal(np.linalg.norm(ys-ys2), 0.0)

	def test_Ylm(self):
		for l in range(4):
			for m in range(-l,1):
				xs = np.linspace(0,1,100)
				ys = np.zeros(10000, np.complex128)
				ys2 = np.zeros(10000, np.complex128)
				t1 = time.time()
				for i in range(xs.shape[0]):
					for j in range(xs.shape[0]):
						ys[i*100+j] = sph_harm(m, l, xs[j]*2*np.pi, xs[i]*np.pi)
				t2 = time.time()
				i,j=0,0
				for i in range(xs.shape[0]):
					for j in range(xs.shape[0]):
						ys2[i*100+j] = PAWC.Ylmr(l, m, c_double(xs[i]*np.pi), c_double(xs[j]*np.pi*2))\
										+ 1.0j*PAWC.Ylmi(l, m, c_double(xs[i]*np.pi), c_double(xs[j]*np.pi*2))
				t3 = time.time()
				#rough speed check
				#assert 2*(t2-t1) > t3 - t2
				assert_almost_equal(np.linalg.norm(ys-ys2),0.0)

	def test_unit_conversion(self):
		struct = Poscar.from_file("CONTCAR").structure
		lattice = numpy_to_cdouble(struct.lattice.matrix.flatten())
		reclattice = numpy_to_cdouble(struct.lattice.reciprocal_lattice.matrix.flatten())
		for site in struct:
			coord = site.coords
			fcoord = site.frac_coords
			ccoord = numpy_to_cdouble(coord)
			cfcoord = numpy_to_cdouble(fcoord)
			PAWC.frac_to_cartesian(cfcoord, lattice)
			temp1 = cdouble_to_numpy(cfcoord, 3)
			assert_almost_equal(np.linalg.norm(temp1-coord), 0.0)
			PAWC.cartesian_to_frac(ccoord, reclattice)
			temp2 = cdouble_to_numpy(ccoord, 3)
			assert_almost_equal(np.linalg.norm(temp2-fcoord), 0.0)

	def test_spline(self):
		vr = self.vr
		cr = self.cr
		struct = Poscar.from_file("CONTCAR").structure
		grid = cr.pps['Ga'].projgrid
		vals = cr.pps['Ga'].realprojs[0]
		rmax = cr.pps['Ga'].rmax / 1.88973
		tst = np.linspace(0, max(grid), 400)
		res1 = scipy.interpolate.CubicSpline(grid, vals, extrapolate=True)(tst)
		x, y = numpy_to_cdouble(grid), numpy_to_cdouble(vals)
		cof = PAWC.spline_coeff(x, y, 100)
		res2 = (c_double * tst.shape[0])()
		for i in range(tst.shape[0]):
			res2[i] = PAWC.proj_interpolate(c_double(tst[i]), c_double(rmax),
						100, x, y, cof)
		res2 = cdouble_to_numpy(res2, tst.shape[0])
		print ('Completed spline test')
		print (res1)
		print (res2)
		print (res1-res2)
		sys.stdout.flush()

	def test_fft3d(self):
		vr = self.vr 
		weights = vr.actual_kpoints_weights
		kws = (c_double * len(weights))()
		for i in range(len(weights)):
			kws[i] = weights[i]
		PAWC.fft_check("WAVECAR".encode('utf-8'), kws,
			numpy_to_cint(np.array([40,40,40])))

	def test_sbt(self):
		from scipy.special import spherical_jn as jn
		cr = CoreRegion(Potcar.from_file("POTCAR"))
		r = cr.pps['Ga'].grid
		ks = 0 * r
		f = cr.pps['Ga'].aewaves[0] - cr.pps['Ga'].pswaves[0];
		print ('yo')
		print (len(r), len(numpy_to_cdouble(r)))
		sys.stdout.flush()
		sbtd = PAWC.spherical_bessel_transform_setup(c_double(520), c_double(5000),
					2, len(r), numpy_to_cdouble(r), numpy_to_cdouble(ks))
		print ('yo')
		sys.stdout.flush()
		res = PAWC.wave_spherical_bessel_transform(sbtd, numpy_to_cdouble(f), 0)
		k = ks[180]
		print ('hi')
		print (r, f)
		vals = jn(0, r * k) * f * r
		integral = np.trapz(vals, r)
		#print (cdouble_to_numpy(ks, 388)[334])
		print (integral)
		print (ks[180])
		print (res[180])
		assert_almost_equal(integral, res[180], decimal=3)

	@nottest
	def test_radial(self):
		# test realspace radial overlap
		# test recipspace radial overlap
		pass


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

"""
	def test_read(self):
		f = open('readtest.out', 'w')
		subprocess.call('valgrind ./memtest read'.split(), stdout=f, stderr=f)
		f.close()

	def test_pseudo(self):
		f = open('pseudotest.out', 'w')
		subprocess.call('valgrind ./memtest pseudo'.split(), stdout=f, stderr=f)
		f.close()

	def test_list(self):
		f = open('listtest.out', 'w')
		subprocess.call('valgrind ./memtest list'.split(), stdout=f, stderr=f)
		f.close()
"""

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
		wf.free_all()
		wf = Wavefunction.from_directory('.', False)
		wf.free_all()
		wf = Wavefunction.from_files('CONTCAR', 'WAVECAR',
			'POTCAR', 'vasprun.xml', 'OUTCAR', True)
		wf.free_all()
		wf = Wavefunction.from_files('CONTCAR', 'WAVECAR',
			'POTCAR', 'vasprun.xml', 'OUTCAR', False)
		wf.free_all()
		wf = Wavefunction.from_files('CONTCAR', 'WAVECAR2.gz',
			'POTCAR', 'vasprun.xml', 'OUTCAR', False)
		wf.free_all()

	def test_writestate(self):
		print("TEST WRITE")
		sys.stdout.flush()
		wf = Wavefunction.from_directory('.')
		fileprefix = ''
		b, k, s = 10, 1, 0
		print(PAWC.write_realspace_state_ri_return.argtypes, "LOOK HERE")
		sys.stdout.flush()
		state1 = wf.write_state_realspace(b, k, s, fileprefix = "", 
			dim=np.array([30,30,30]), return_wf = True)
		wf.free_all()
		wf = Wavefunction.from_directory('.', False)
		state2 = wf.write_state_realspace(b, k, s, fileprefix = "", 
			dim=np.array([30,30,30]), return_wf = True)
		wf.free_all()
		assert_almost_equal(np.linalg.norm(state1-state2),0)
		assert state1.shape[0] == 2*30*30*30
		filename_base = "%sB%dK%dS%d" % (fileprefix, b, k, s)
		filename1 = "%s_REAL" % filename_base
		filename2 = "%s_IMAG" % filename_base
		os.remove(filename1)
		os.remove(filename2)

	def test_density(self):
		print("TEST DENSITY")
		sys.stdout.flush()
		wf = Wavefunction.from_directory('.')
		wf.write_density_realspace(dim=np.array([40,40,40]))
		tstchg = Chgcar.from_file("AECCAR2").data['total']# / wf.structure.volume
		chg = Chgcar.from_file("PYAECCAR").data['total']
		reldiff = np.sqrt(np.mean(((chg-tstchg)/tstchg)**2))
		newchg = chg-tstchg
		Chgcar(Poscar(wf.structure), {'total': newchg}).write_file('DIFFCHGCAR.vasp')
		print(np.sum(chg)/40**3, np.sum(tstchg)/40**3)
		assert_almost_equal(reldiff, 0, decimal=3)
		#os.remove('PYAECCAR')

	def test_pseudoprojector(self):
		print("TEST PSEUDO")
		sys.stdout.flush()
		# test ps projections
		wf = Wavefunction.from_directory('.')
		basis = Wavefunction.from_directory('.')
		pr = Projector(wf, basis, pseudo = True)
		res = pr.single_band_projection(6)
		assert res.shape[0] == basis.nband * basis.nspin * basis.nwk
		res = pr.defect_band_analysis(4, 10, False)
		assert len(res.keys()) == 15
		wf.free_all()
		basis.free_all()

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
		basis.free_all()
		wf1.free_all()
		pr.free_all()

		generator = Projector.setup_multiple_projections('.', ['.', '.'])
		for wf_dir, wf in generator:
			wf.defect_band_analysis(4, 10, spinpol=True)

	def test_projector_gz(self):
		print("TEST PROJGZ")
		sys.stdout.flush()
		# test ae projections
		wf1 = Wavefunction.from_files('CONTCAR', 'WAVECAR2.gz',
			'POTCAR', 'vasprun.xml', 'OUTCAR', False)
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
		basis.free_all()
		wf1.free_all()
		pr.free_all()

		generator = Projector.setup_multiple_projections('.', ['.', '.'])
		for wf_dir, wf in generator:
			wf.defect_band_analysis(4, 10, spinpol=True)

	def test_offsite(self):
		Projector = DummyProjector
		print("TEST OFFSITE")
		sys.stdout.flush()
		wf1 = Wavefunction.from_directory('.', False)
		basis = Wavefunction.from_directory('.', False)
		pr = Projector(wf1, basis)
		for b in range(wf1.nband):
			v, c = pr.proportion_conduction(b)
			print("CHECK_VALS", v,c)
			if b < 6:
				assert_almost_equal(v, 1, decimal=2)
				assert_almost_equal(c, 0, decimal=4)
			else:
				assert_almost_equal(v, 0, decimal=4)
				assert_almost_equal(c, 1, decimal=2)
		basis.free_all()
		wf1.free_all()
		pr.free_all()

		generator = Projector.setup_multiple_projections('.', ['.', '.'])
		for wf_dir, wf in generator:
			wf.defect_band_analysis(4, 10, spinpol=True)
