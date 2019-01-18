from pawpyseed.core cimport pawpyc
from cpython cimport array
from libc.stdlib cimport malloc, free
from libc.stdio cimport FILE
from pymatgen.core.structure import Structure
from pymatgen.io.vasp.outputs import Vasprun
#from pawpyseed.core import projector
#from pawpyseed.core import wavefunction
from monty.io import zopen
import numpy as np
cimport numpy as np
import time

class Timer:
	@staticmethod
	def setup_time(t):
		Timer.ALL_SETUP_TIME += t
	@staticmethod
	def overlap_time(t):
		Timer.ALL_OVERLAP_TIME += t
	@staticmethod
	def augmentation_time(t):
		Timer.ALL_AUGMENTATION_TIME += t
	ALL_SETUP_TIME = 0
	ALL_OVERLAP_TIME = 0
	ALL_AUGMENTATION_TIME = 0

def el(site):
	"""
	Return the element symbol of a pymatgen
	site object
	"""
	return site.specie.symbol

cdef class Pointer:

	cdef void* ptr

	cdef set_pointer(self, void* ptr):
		self.ptr = ptr

	@staticmethod
	cdef from_pointer(void* ptr):
		self = Pointer(0)
		self.set_pointer(ptr)
		return self

	def __dealloc__(self):
		free(self.ptr)


cdef class PseudoWavefunction:
	"""
	THIS CLASS IS NOT USEFUL ON ITS OWN. IF YOU WANT TO WORK WITH
	PSEUDOWAVEFUNCTIONS, INITIALIZE A Wavefunction OBJECT WITH
	setup_projectors=False.

	Class for storing pseudowavefunction from WAVECAR file. Most important attribute
	is wf_ptr, a C pointer used in the C portion of the program for storing
	plane wave coefficients.

	Attributes:
		kpts (np.array): nx3 array of fractional kpoint vectors,
			where n is the number of kpoints
		kws (np.array): weight of each kpoint
		wf_ptr (ctypes POINTER): c pointer to pswf_t object
		ncl (bool): Whether the pseudowavefunction is from a noncollinear
			VASP calculation
		band_props (list): [band gap, conduction band minimum,
			valence band maximum, whether the band gap is direct]
	"""
	cdef pawpyc.pswf_t* wf_ptr
	cdef readonly int nband
	cdef readonly int nwk
	cdef readonly int nspin
	cdef readonly int ncl
	cdef readonly np.ndarray kws
	cdef readonly np.ndarray kpts
	cdef readonly np.ndarray band_props

	def __init__(self, filename="WAVECAR", vr="vasprun.xml"):
		if type(vr) == str:
			vr = Vasprun(vr)
		weights = vr.actual_kpoints_weights
		kpts = vr.actual_kpoints
		self.kws = np.array(weights)
		self.kpts = np.array(kpts)
		cdef double[::1] kws = self.kws
		if '.gz' in filename or '.bz2' in filename:
			f = zopen(filename, 'rb')
			contents = f.read()
			f.close() 
			self.wf_ptr = pawpyc.read_wavefunctions_from_str(
				contents, &kws[0])
		else:
			self.wf_ptr = pawpyc.read_wavefunctions(filename.encode('utf-8'), &kws[0])
		self.ncl = pawpyc.is_ncl(self.wf_ptr) > 0
		self.band_props = np.array(vr.eigenvalue_band_properties)
		self.nband = pawpyc.get_nband(self.wf_ptr)
		self.nwk = pawpyc.get_nwk(self.wf_ptr)
		self.nspin = pawpyc.get_nspin(self.wf_ptr)
		self.encut = pawpyc.get_encut(self.wf_ptr)

	def __dealloc__(self):
		pawpyc.free_pswf(self.wf_ptr)

	def pseudoprojection(self, band_num, PseudoWavefunction basis):
		"""
		Computes <psibt_n1k|psit_n2k> for all n1 and k
		and a given n2, where psibt are basis structures
		pseudowavefunctions and psit are self pseudowavefunctions

		Arguments:
			band_num (int): n2 (see description)
			basis (Pseudowavefunction): pseudowavefunctions onto whose bands
				the band of self is projected
		"""
		res = np.zeros(self.nband * self.nwk * self.nspin, dtype = np.complex128)
		cdef double complex[::1] resv = res
		pawpyc.pseudoprojection(&resv[0], basis.wf_ptr, self.wf_ptr, band_num)
		return res

cdef class CWavefunction(PseudoWavefunction):

	cdef int[::1] dimv
	cdef int gridsize

	cdef int[::1] nums
	cdef double[::1] coords
	cdef int number_projector_elements
	cdef readonly int projector_owner
	
	cdef pawpyc.ppot_t* projector_list

	def __init__(self, filename="WAVECAR", vr="vasprun.xml"):
		self.projector_owner = 0
		super(CWavefunction, self).__init__(filename, vr)

	def __dealloc__(self):
		if self.projector_owner:
			pawpyc.free_ppot_list(self.projector_list, self.number_projector_elements)

	def _c_projector_setup(self, int num_elems, int num_sites,
							double grid_encut, nums, coords, dim, pps):
		"""
		Returns a point to a list of ppot_t objects in C,
		to be used for high performance parts of the code

		Args:
			pps (dict of Pseudopotential objects): keys are integers,
				values of Pseudopotential objects

		Returns:
			c_void_p object pointing to ppot_t list with each Pseudopotential,
			ordered in the list by their numerical keys
		"""

		start = time.monotonic()
		clabels = np.array([], np.intc)
		ls = np.array([], np.intc)
		projectors = np.array([], np.double)
		aewaves = np.array([], np.double)
		pswaves = np.array([], np.double)
		wgrids = np.array([], np.double)
		augs = np.array([], np.double)
		rmaxs = np.array([], np.double)

		for num in sorted(pps.keys()):
			pp = pps[num]
			clabels = np.append(clabels, [num, len(pp.ls), pp.ndata, len(pp.grid)])
			rmaxs = np.append(rmaxs, pp.rmax)
			ls = np.append(ls, pp.ls)
			wgrids = np.append(wgrids, pp.grid)
			augs = np.append(augs, pp.augs)
			for i in range(len(pp.ls)):
				proj = pp.realprojs[i]
				aepw = pp.aewaves[i]
				pspw = pp.pswaves[i]
				projectors = np.append(projectors, proj)
				aewaves = np.append(aewaves, aepw)
				pswaves = np.append(pswaves, pspw)

		cdef int[::1] clabels_v = clabels.astype(np.intc)
		cdef int[::1] ls_v = ls.astype(np.intc)
		cdef double[::1] wgrids_v = wgrids.astype(np.double)
		cdef double[::1] projectors_v = projectors.astype(np.double)
		cdef double[::1] aewaves_v = aewaves.astype(np.double)
		cdef double[::1] pswaves_v = pswaves.astype(np.double)
		cdef double[::1] rmaxs_v = rmaxs.astype(np.double)

		print ("GRID ENCUT", grid_encut)
		self.projector_list = pawpyc.get_projector_list(
							num_elems, &clabels_v[0], &ls_v[0], &wgrids_v[0],
							&projectors_v[0], &aewaves_v[0], &pswaves_v[0],
							&rmaxs_v[0], grid_encut)
		end = time.monotonic()
		print('--------------\nran get_projector_list in %f seconds\n---------------' % (end-start))

		self.projector_owner = 1
		self.number_projector_elements = num_elems
		self.nums = np.array(nums, dtype = np.int32, copy = True)
		self.coords = np.array(coords, dtype = np.float64, copy = True)
		self.update_dimv(dim)

		pawpyc.setup_projections(
			self.wf_ptr, self.projector_list,
			num_elems, num_sites, &self.dimv[0],
			&self.nums[0], &self.coords[0]
			)

	def update_dimv(self, dim):
		dim = np.array(dim, dtype = np.int32, order = 'C', copy = False)
		self.dimv = dim
		self.gridsize = np.cumprod(dim)[-1]

	def _get_realspace_state(self, int b, int k, int s):
		res = np.zeros(self.gridsize, dtype = np.complex128, order='C')
		cdef double complex[::1] resv = res
		pawpyc.realspace_state(&resv[0], b, k+s*self.nwk,
			self.wf_ptr, self.projector_list,
			&self.dimv[0], &self.nums[0], &self.coords[0])
		res.shape = self.dimv
		return res

	def _get_realspace_density(self):
		res = np.zeros(self.gridsize, dtype = np.float64, order='C')
		cdef double[::1] resv = res
		pawpyc.ae_chg_density(&resv[0], self.wf_ptr, self.projector_list,
			&self.dimv[0], &self.nums[0], &self.coords[0])
		res.shape = self.dimv
		return res

	def _write_realspace_state(self, filename1, filename2, double scale, int b, int k, int s):
		filename1 = bytes(filename1.encode('utf-8'))
		filename2 = bytes(filename2.encode('utf-8'))
		res = self._get_realspace_state(b, k, s)
		res2 = res.view()
		res2.shape = self.gridsize
		resr = np.ascontiguousarray(np.real(res2))
		resi = np.ascontiguousarray(np.imag(res2))
		cdef double[::1] resv = resr
		pawpyc.write_volumetric(filename1, &resv[0], &self.dimv[0], scale)
		resv = resi
		pawpyc.write_volumetric(filename2, &resv[0], &self.dimv[0], scale)
		return res

	def _write_realspace_density(self, filename, double scale):
		filename = bytes(filename.encode('utf-8'))
		res = self._get_realspace_density()
		res2 = res.view()
		res2.shape = self.gridsize
		cdef double[::1] resv = res2
		pawpyc.write_volumetric(filename, &resv[0], &self.dimv[0], scale);
		return res

"""
if not arr.flags['C_CONTIGUOUS']:
	arr = np.ascontiguousarray(arr)
"""


def py_compensation_terms(b, wf, basis, M_R, M_S, N_R, N_S, N_RS_R, N_RS_S, dim):
	cdef array.array ref_labels = array.array(np.int32, basis.nums)
	cdef array.array proj_labels = array.array(np.int32, wf.nums)
	cdef array.array ref_coords = array.array(np.float64, basis.coords)
	cdef array.array proj_coords = array.array(np.float64, wf.coords)
	cdef array.array fftg = array.array('i', dim)
	return


	