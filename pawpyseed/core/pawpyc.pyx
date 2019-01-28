# cython : profile=True

from pawpyseed.core cimport pawpyc
from cpython cimport array
from libc.stdlib cimport malloc, free
from libc.stdio cimport FILE
from pymatgen.core.structure import Structure
from pymatgen.io.vasp.outputs import Vasprun
from monty.io import zopen
import numpy as np
from numpy.testing import assert_almost_equal
cimport numpy as np
import time
import sys
from libc.stdint cimport uintptr_t
from pawpyseed.core.symmetry import *

###################
#  TIMER SECTION  #
###################

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

################################
#  MISCELLANEOUS PYTHON UTILS  #
################################  

def el(site):
	"""
	Return the element symbol of a pymatgen
	site object
	"""
	return site.specie.symbol

OPSIZE = 9

def make_c_ops(op_nums, symmops):
	ops = np.zeros(OPSIZE*len(op_nums), dtype = np.float64)
	for i in range(len(op_nums)):
		ops[OPSIZE*i:OPSIZE*(i+1)] = symmops[op_nums[i]].rotation_matrix.flatten()
	drs = np.zeros(3*len(op_nums), dtype = np.float64)
	for i in range(len(op_nums)):
		drs[3*i:3*(i+1)] = symmops[op_nums[i]].translation_vector
	return ops, drs

#################################
#  C UTILS INTERFACE FUNCTIONS  #
#################################

cpdef double legendre(int l, int m, double x):
	return pawpyc.legendre(l, m, x)

cpdef double complex Ylm(int l, int m, double theta, double phi):
	return pawpyc.Ylm(l, m, theta, phi)

cpdef double complex Ylm2(int l, int m, double costheta, double phi):
	return pawpyc.Ylm2(l, m, costheta, phi)

cpdef frac_to_cartesian(np.ndarray[double, ndim=1] coord,
	np.ndarray[double, ndim=2] lattice):

	if not coord.flags['C_CONTIGUOUS']:
		coord = np.ascontiguousarray(coord)
	latticef = lattice.flatten()
	if not latticef.flags['C_CONTIGUOUS']:
		latticef = np.ascontiguousarray(latticef)
	cdef double[::1] coordv = coord
	cdef double[::1] latticev = latticef
	pawpyc.frac_to_cartesian(&coordv[0], &latticev[0])

cpdef cartesian_to_frac(np.ndarray[double, ndim=1] coord,
	np.ndarray[double, ndim=2] reclattice):

	if not coord.flags['C_CONTIGUOUS']:
		coord = np.ascontiguousarray(coord)
	reclatticef = reclattice.flatten()
	if not reclatticef.flags['C_CONTIGUOUS']:
		reclatticef = np.ascontiguousarray(reclatticef)
	cdef double[::1] coordv = coord
	cdef double[::1] latticev = reclatticef
	pawpyc.cartesian_to_frac(&coordv[0], &latticev[0])

cpdef interpolate(np.ndarray[double, ndim=1] res,
	np.ndarray[double, ndim=1] tst,
	np.ndarray[double, ndim=1] x,
	np.ndarray[double, ndim=1] y,
	double rmax,
	int size, int tstsize):

	cdef double[::1] xv = x
	cdef double[::1] yv = y
	cdef double** coef = pawpyc.spline_coeff(&xv[0], &yv[0], size)

	for i in range(tstsize):
		res[i] = pawpyc.proj_interpolate(tst[i], rmax,
			size, &xv[0], &yv[0], coef)

	free(coef[0])
	free(coef[1])
	free(coef[2])
	free(coef)

cpdef fft_check(str wavecar, np.ndarray[double, ndim=1] kpt_weights,
	np.ndarray[int, ndim=1] fftgrid):

	cdef int[::1] fftg = fftgrid
	cdef double[::1] kws = kpt_weights
	cdef pawpyc.pswf_t* wf = pawpyc.read_wavefunctions(wavecar.encode('utf-8'), &kws[0])
	cdef double complex* x = <double complex*> malloc(fftg[0]*fftg[1]*fftg[2]*sizeof(double complex))
	pawpyc.fft3d(x, wf.G_bounds, wf.lattice, wf.kpts[0].k, wf.kpts[0].Gs,
		wf.kpts[0].bands[0].Cs, wf.kpts[0].bands[0].num_waves, &fftg[0])
	cdef int* Gs = wf.kpts[0].Gs
	cdef float complex* Cs = wf.kpts[0].bands[0].Cs
	cdef double inv_sqrt_vol = np.power(pawpyc.determinant(wf.lattice), -0.5)
	cdef double dv = pawpyc.determinant(wf.lattice) / fftg[0] / fftg[1] / fftg[2]
	cdef double* kpt = wf.kpts[0].k;
	cdef double f1 = 0
	cdef double f2 = 0
	cdef double f3 = 0
	cdef int ind
	cdef double complex temp
	il = np.arange(fftg[0], dtype=np.float64) / fftg[0]
	jl = np.arange(fftg[1], dtype=np.float64) / fftg[1]
	kl = np.arange(fftg[2], dtype=np.float64) / fftg[2]
	for i in np.arange(fftg[0]):
		for j in np.arange(fftg[1]):
			for k in np.arange(fftg[2]):
				f1 = il[i]
				f2 = jl[j]
				f3 = kl[k]
				temp = 0;
				for w in range(wf.kpts[0].bands[0].num_waves):
					temp += Cs[w] * np.exp((f1 * (Gs[3*w]) +
							f2 * (Gs[3*w+1]) +
							f3 * (Gs[3*w+2])) * 2.0j * np.pi)
				temp *= inv_sqrt_vol
				ind = i*fftg[1]*fftg[2]+j*fftg[2]+k
				assert_almost_equal(np.abs(x[ind] - temp), 0)
	free(x)
	pawpyc.free_pswf(wf)

cpdef spherical_bessel_transform(double encut, int l,
	np.ndarray[double, ndim=1] r, np.ndarray[double, ndim=1] f):

	size = r.shape[0]
	if not r.flags['C_CONTIGUOUS']:
		r = np.ascontiguousarray(r)
	if not f.flags['C_CONTIGUOUS']:
		f = np.ascontiguousarray(f)
	cdef np.ndarray k = np.zeros(size, dtype=np.float64, order='C')
	cdef np.ndarray fk = np.zeros(size, dtype=np.float64, order='C')
	cdef double[::1] rv = r
	cdef double[::1] kv = k
	cdef double[::1] fv = f
	cdef pawpyc.sbt_descriptor_t* sbtd = pawpyc.spherical_bessel_transform_setup(encut, 0,
					l, size, &rv[0], &kv[0])
	cdef double* fkv = pawpyc.wave_spherical_bessel_transform(sbtd, &fv[0], l)
	for i in range(size):
		fk[i] = fkv[i]
	free(fkv)
	return k, fk

cpdef reciprocal_offsite_wave_overlap(np.ndarray[double, ndim=1] dcoord,
	np.ndarray[double, ndim=1] r1, np.ndarray[double, ndim=1] f1,
	np.ndarray[double, ndim=1] r2, np.ndarray[double, ndim=1] f2,
	int l1, int m1, int l2, int m2):

	cdef int size1 = r1.shape[0]
	cdef int size2 = r2.shape[0]

	cdef double[::1] r1v = r1
	cdef double[::1] r2v = r2
	cdef double[::1] f1v = f1
	cdef double[::1] f2v = f2

	cdef double encut = 1e5
	k1, fk1 = spherical_bessel_transform(encut, l1, r1, f1)
	k2, fk2 = spherical_bessel_transform(encut, l2, r2, f2)

	cdef double[::1] dcoordv = dcoord
	cdef double[::1] k1v = k1
	cdef double[::1] fk1v = fk1
	cdef double[::1] k2v = k2
	cdef double[::1] fk2v = fk2

	cdef double** s1 = pawpyc.spline_coeff(&k1v[0], &fk1v[0], size1)
	cdef double** s2 = pawpyc.spline_coeff(&k2v[0], &fk2v[0], size2)

	res = pawpyc.reciprocal_offsite_wave_overlap(&dcoordv[0],
		&k1v[0], &fk1v[0], s1, size1,
		&k2v[0], &fk2v[0], s2, size2,
		NULL, l1, m1, l2, m2)

	free(s1[0])
	free(s1[1])
	free(s1[2])
	free(s1)
	free(s2[0])
	free(s2[1])
	free(s2[2])
	free(s2)

	return res

############################
#  PAWPYSEED BASE CLASSES  #
############################

cdef class PWFPointer:

	cdef pawpyc.pswf_t* ptr
	cdef readonly np.ndarray kpts
	cdef readonly np.ndarray weights
	cdef readonly np.ndarray band_props

	def __init__(self, filename = None, vr = None):
		cdef double[::1] kws
		if filename == None or vr == None:
			self.ptr = NULL
		else:
			filename = str(filename)
			if type(vr) == str:
				vr = Vasprun(vr) 
			self.weights = np.array(vr.actual_kpoints_weights, dtype=np.float64)
			self.kpts = np.array(vr.actual_kpoints, dtype=np.float64)
			self.band_props = np.array(vr.eigenvalue_band_properties)
			kws = self.weights
			if '.gz' in filename or '.bz2' in filename:
				f = zopen(filename, 'rb')
				contents = f.read()
				f.close()
				self.ptr = pawpyc.read_wavefunctions_from_str(
					contents, &kws[0])
			else:
				self.ptr = pawpyc.read_wavefunctions(filename.encode('utf-8'), &kws[0])

	@staticmethod
	cdef PWFPointer from_pointer_and_kpts(pawpyc.pswf_t* ptr,
		structure, kpts, band_props, allkpts = None, weights = None):

		return_kpts_and_weights = False
		if (allkpts is None) or (weights is None):
			return_kpts_and_weights = True
			allkpts, orig_kptnums, op_nums, symmops, trs = get_nosym_kpoints(kpts, structure)
			weights = np.ones(allkpts.shape[0], dtype=np.float64, order='C')
			# need to change this if spin orbit coupling is added in later
			for i in range(allkpts.shape[0]):
				if np.linalg.norm(allkpts[i]) < 1e-10:
					weights[i] *= 0.5
			weights /= np.sum(weights)
		else:
			orig_kptnums, op_nums, symmops, trs = get_kpt_mapping(allkpts, kpts, structure)

		ops, drs = make_c_ops(op_nums, symmops)

		kpts = np.array(allkpts, np.float64, order='C', copy=False)
		weights = np.array(weights, np.float64, order='C', copy=False)
		print(kpts, weights, "hihi")
		sys.stdout.flush()
		cdef double[::1] weights_v = weights
		cdef int[::1] orig_kptnums_v = np.array(orig_kptnums, np.int32, order='C', copy=False)
		cdef int[::1] op_nums_v = np.array(op_nums, np.int32, order='C', copy=False)
		cdef double[::1] ops_v = np.array(ops, np.float64, order='C', copy=False)
		cdef double[::1] drs_v = np.array(drs, np.float64, order='C', copy=False)
		cdef int[::1] trs_v = np.array(trs, np.int32, order='C', copy=False)

		cdef pawpyc.pswf_t* new_ptr = pawpyc.expand_symm_wf(ptr, len(orig_kptnums),
				&orig_kptnums_v[0], &ops_v[0], &drs_v[0], &weights_v[0], &trs_v[0])

		cdef PWFPointer pwfp = PWFPointer()
		pwfp.ptr = new_ptr
		pwfp.kpts = kpts
		pwfp.weights = weights
		pwfp.band_props = np.array(band_props)
		return pwfp


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

	def __init__(self, PWFPointer pwf):
		if pwf.ptr is NULL:
			raise Exception("NULL PWFPointer ptr!")
		self.wf_ptr = pwf.ptr
		self.kpts = pwf.kpts.copy(order='C')
		self.kws = pwf.weights.copy(order='C')
		self.ncl = pawpyc.is_ncl(self.wf_ptr) > 0
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
		res = np.zeros(basis.nband * basis.nwk * basis.nspin, dtype = np.complex128)
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

	def __init__(self, PWFPointer pwf):
		self.projector_owner = 0
		super(CWavefunction, self).__init__(pwf)

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

	def _desymmetrized_pwf(self, structure, band_props, allkpts = None, weights = None):
		return PWFPointer.from_pointer_and_kpts(<pawpyc.pswf_t*> self.wf_ptr, structure,
							self.kpts, band_props, allkpts, weights)

	def _get_occs(self):
		nk = self.nwk * self.nspin
		res = np.zeros(self.nband * nk, dtype=np.float64, order='C')
		for k in range(nk):
			for b in range(self.nband):
				res[b * nk + k] = self.wf_ptr.kpts[k].bands[b].occ
		return res

	def _get_energy_list(self, bands):

		energy_list = {}
		for b in bands:
			energy_list[b] = []
			for k in range(self.nwk):
				for s in range(self.nspin):
					energy_list[b].append([pawpyc.get_energy(self.wf_ptr, b, k, s),\
										pawpyc.get_occ(self.wf_ptr, b, k, s)])


cdef class CProjector:

	cdef public CWavefunction wf
	cdef public CWavefunction basis

	cdef int[::1] M_R
	cdef int[::1] M_S
	cdef int[::1] N_R
	cdef int[::1] N_S
	cdef int[::1] N_RS_R
	cdef int[::1] N_RS_S

	cdef int num_M_R
	cdef int num_M_S 
	cdef int num_N_R
	cdef int num_N_S
	cdef int num_N_RS_R
	cdef int num_N_RS_S

	def __init__(self, wf, basis):
		self.wf = wf
		self.basis = basis

	def _setup_overlap(self, site_cat):
		self.M_R = np.array(site_cat[0], dtype=np.int32, order = 'C')
		self.M_S = np.array(site_cat[1], dtype=np.int32, order = 'C')
		self.N_R = np.array(site_cat[2], dtype=np.int32, order = 'C')
		self.N_S = np.array(site_cat[3], dtype=np.int32, order = 'C')
		self.N_RS_R = np.array(site_cat[4], dtype=np.int32, order = 'C')
		self.N_RS_S = np.array(site_cat[5], dtype=np.int32, order = 'C')

		self.num_M_R, self.num_M_S = len(site_cat[0]), len(site_cat[1])
		self.num_N_R, self.num_N_S = len(site_cat[2]), len(site_cat[3])
		self.num_N_RS_R, self.num_N_RS_S = len(site_cat[4]), len(site_cat[5])

		cdef int* M_R = NULL if self.num_M_R == 0 else &self.M_R[0]
		cdef int* M_S = NULL if self.num_M_S == 0 else &self.M_S[0]
		cdef int* N_R = NULL if self.num_N_R == 0 else &self.N_R[0]
		cdef int* N_S = NULL if self.num_N_S == 0 else &self.N_S[0]
		cdef int* N_RS_R = NULL if self.num_N_RS_R == 0 else &self.N_RS_R[0]
		cdef int* N_RS_S = NULL if self.num_N_RS_S == 0 else &self.N_RS_S[0]
		
		pawpyc.overlap_setup_real(self.basis.wf_ptr, self.wf.wf_ptr,
			&self.basis.nums[0], &self.wf.nums[0], &self.basis.coords[0], &self.wf.coords[0],
			N_R, N_S, N_RS_R, N_RS_S,
			#self.N_R, self.N_S, self.N_RS_R, self.N_RS_S,
			self.num_N_R, self.num_N_S, self.num_N_RS_R)

	def _add_augmentation_terms(self, np.ndarray res, band_num):
		# declare res more specifically
		cdef double complex[::1] resv = res

		cdef int* M_R = NULL if self.num_M_R == 0 else &self.M_R[0]
		cdef int* M_S = NULL if self.num_M_S == 0 else &self.M_S[0]
		cdef int* N_R = NULL if self.num_N_R == 0 else &self.N_R[0]
		cdef int* N_S = NULL if self.num_N_S == 0 else &self.N_S[0]
		cdef int* N_RS_R = NULL if self.num_N_RS_R == 0 else &self.N_RS_R[0]
		cdef int* N_RS_S = NULL if self.num_N_RS_S == 0 else &self.N_RS_S[0]

		pawpyc.compensation_terms(&resv[0], band_num, self.wf.wf_ptr, self.basis.wf_ptr,
			self.num_M_R, self.num_N_R, self.num_N_S, self.num_N_RS_R,
			M_R, M_S, N_R, N_S, N_RS_R, N_RS_S,
			#self.M_R, self.M_S, self.N_R, self.N_S, self.N_RS_R, self.N_RS_S,
			&self.wf.nums[0], &self.wf.coords[0], &self.basis.nums[0], &self.basis.coords[0],
			&self.wf.dimv[0])

	def _realspace_projection(self, int band_num, np.ndarray dim):
		res = np.zeros(self.basis.nband * self.basis.nwk * self.basis.nspin,
			dtype=np.complex128, order='C')
		cdef double complex[::1] resv = res
		cdef int[::1] dimv
		if dim == None:
			dimv = self.wf.dimv
		else:
			dimv = np.array(dim, dtype=np.float64, order='C', copy=False)
		pawpyc.project_realspace_state(&resv[0], 
			band_num, self.wf.wf_ptr, self.basis.wf_ptr,
			self.wf.projector_list, self.basis.projector_list,
			&dimv[0], &self.wf.nums[0], &self.wf.coords[0],
			&self.basis.nums[0], &self.basis.coords[0])
		return res

