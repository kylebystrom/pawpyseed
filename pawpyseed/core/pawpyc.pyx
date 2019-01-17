from pawpyseed.core cimport pawpyc
from cpython cimport array
from libc.stdlib cimport malloc, free
from libc.stdio cimport FILE
from pymatgen.io.vasp.outputs import Vasprun
#from pawpyseed.core import projector
#from pawpyseed.core import wavefunction
from monty.io import zopen
import numpy as np
cimport numpy as np

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


	