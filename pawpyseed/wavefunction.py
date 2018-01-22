from pymatgen.io.vasp.inputs import Potcar, Poscar
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.core.structure import Structure
import numpy as np
from ctypes import *
import os
import numpy as np
import json

import sys
sys.stdout.flush()

MODULE_DIR = os.path.dirname(os.path.abspath(__file__))
PAWC = CDLL(os.path.join(MODULE_DIR, "pawpy.so"))

PAWC.read_wavefunctions.restype = POINTER(None)
PAWC.get_projector_list.restype = POINTER(None)
PAWC.read_wavefunctions.restype = POINTER(None)
PAWC.compensation_terms.restype = POINTER(c_double)
PAWC.get_occs.restype = POINTER(c_double)
PAWC.get_nband.restype = c_int
PAWC.get_nwk.restype = c_int
PAWC.get_nspin.restype = c_int

PAWC.free_ptr.restype = None
PAWC.free_ppot_list.restype = None
PAWC.free_pswf.restype = None

def cdouble_to_numpy(arr, length):
	arr = cast(arr, POINTER(c_double))
	newarr = np.zeros(length)
	for i in range(length):
		newarr[i] = arr[i]
	PAWC.free_ptr(arr)
	return newarr

def cfloat_to_numpy(arr, length):
	arr = cast(arr, POINTER(c_float))
	newarr = np.zeros(length)
	for i in range(length):
		newarr[i] = arr[i]
	PAWC.free_ptr(arr)
	return newarr

def cfloat_to_numpy(arr, length):
	arr = cast(arr, POINTER(c_int))
	newarr = np.zeros(length)
	for i in range(length):
		newarr[i] = arr[i]
	PAWC.free_ptr(arr)
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

def el(site):
	return site.specie.symbol

class Pseudopotential:
	"""
	Contains important attributes from a VASP pseudopotential files. POTCAR
	"settings" can be read from the pymatgen POTCAR object
	"""

	def __init__(self, data, rmax):
		nonradial, radial = data.split("PAW radial sets", 1)
		partial_waves = radial.split("pseudo wavefunction")
		gridstr, partial_waves = partial_waves[0], partial_waves[1:]
		self.rmax = rmax
		self.pswaves = []
		self.aewaves = []
		self.recipprojs = []
		self.realprojs = []
		self.nonlocalprojs = []
		self.ls = []

		auguccstr, gridstr = gridstr.split("grid", 1)
		gridstr, aepotstr = gridstr.split("aepotential", 1)
		aepotstr, corechgstr = aepotstr.split("core charge-density", 1)
		corechgstr, kenstr = corechgstr.split("kinetic energy-density", 1)
		kenstr, pspotstr = kenstr.split("pspotential", 1)
		pspotstr, pscorechgstr = pspotstr.split("core charge-density (pseudized)", 1)
		self.grid = self.make_nums(gridstr+' 0')
		self.aepotential = self.make_nums(aepotstr)
		self.aecorecharge = self.make_nums(corechgstr)
		self.kinetic = self.make_nums(kenstr)
		self.pspotential = self.make_nums(pspotstr)
		self.pscorecharge = self.make_nums(pscorechgstr)

		for pwave in partial_waves:
			lst = pwave.split("ae wavefunction", 1)
			self.pswaves.append(self.make_nums(lst[0]+' 0'))
			self.aewaves.append(self.make_nums(lst[1]+' 0'))

		projstrs = nonradial.split("Non local Part")
		topstr, projstrs = projstrs[0], projstrs[1:]
		self.T = float(topstr[-22:-4])
		topstr, atpschgstr = topstr[:-22].split("atomic pseudo charge-density", 1)
		topstr, corechgstr = topstr.split("core charge-density (partial)", 1)
		settingstr, localstr = topstr.split("local part", 1)
		localstr, self.gradxc = localstr.split("gradient corrections used for XC", 1)
		self.gradxc = int(self.gradxc)
		self.localpart = self.make_nums(localstr)
		self.localnum = self.localpart[0]
		self.localpart = self.localpart[1:]
		self.coredensity = self.make_nums(corechgstr)
		self.atomicdensity = self.make_nums(atpschgstr)

		for projstr in projstrs:
			lst = projstr.split("Reciprocal Space Part")
			nonlocalvals, projs = lst[0], lst[1:]
			nonlocalvals = self.make_nums(nonlocalvals)
			l = nonlocalvals[0]
			count = nonlocalvals[1]
			self.nonlocalprojs.append(nonlocalvals[2:])
			for proj in projs:
				recipproj, realproj = proj.split("Real Space Part")
				self.recipprojs.append(self.make_nums(recipproj))
				self.realprojs.append(self.make_nums(realproj))
				self.ls.append(l)

		settingstr, projgridstr = settingstr.split("STEP   =")
		self.ndata = int(settingstr.split()[-1])
		projgridstr = projgridstr.split("END")[0]
		self.projgrid = self.make_nums(projgridstr)
		self.step = (self.projgrid[0], self.projgrid[1])

		self.projgrid = np.linspace(0,rmax/1.88973,self.ndata,False,dtype=np.float64)

	def make_nums(self, numstring):
		return np.array([float(num) for num in numstring.split()])

class CoreRegion:
	"""
	List of Pseudopotential objects to describe the core region of a structure.
	"""

	def __init__(self, potcar):
		self.pps = {}
		for potsingle in potcar:
			self.pps[potsingle.element] = Pseudopotential(potsingle.data[:-15], potsingle.rmax)
		

class PseudoWavefunction:
	"""
	Class for storing pseudowavefunction from WAVECAR file. Most important attribute
	is wf_ptr, a C pointer used in the C portion of the program for storing
	plane wave coefficients
	"""

	def __init__(self, filename="WAVECAR", vr="vasprun.xml"):
		if type(vr) == str:
			vr = Vasprun(vr)
		weights = vr.actual_kpoints_weights
		kws = (c_double * len(weights))()
		for i in range(len(weights)):
			kws[i] = weights[i]
		self.kws = weights
		self.wf_ptr = PAWC.read_wavefunctions(filename.encode('utf-8'), byref(kws))

class Wavefunction:
	"""
	Class for storing and manipulating all electron wave functions in the PAW
	formalism.

	Attributes:
		structure (pymatgen.core.structure.Structure): stucture of the material
			that the wave function describes
		pwf (PseudoWavefunction): Pseudowavefunction componenet
		cr (CoreRegion): Contains the pseudopotentials, with projectors and
			partials waves, for the structure
		projector: ctypes object for interfacing with C code
	"""

	def __init__(self, struct, pwf, cr, dim):
		"""
		Arguments:
			struct (pymatgen.core.Structure): structure that the wavefunction describes
			pwf (PseudoWavefunction): Pseudowavefunction componenet
			cr (CoreRegion): Contains the pseudopotentials, with projectors and
				partials waves, for the structure
		Returns:
			Wavefunction object
		"""
		self.structure = struct
		self.pwf = pwf
		self.cr = cr
		self.projector = PAWC
		self.dim = np.array(dim);
		self.projector_list = None

	def from_files(self, struct="CONTCAR", pwf="WAVECAR", cr="POTCAR", vr="vasprun.xml"):
		"""
		Construct a Wavefunction object from file paths.
		Arguments:
			struct (str): VASP POSCAR or CONTCAR file path
			pwf (str): VASP WAVECAR file path
			vr (str): VASP vasprun file path
		Returns:
			Wavefunction object
		"""
		return Wavefunction(Poscar.from_file(struct), PseudoWavefunction(pwf, vr), Potcar.from_file(cr))

	def make_site_lists(self, basis):
		"""
		Organizes sites into sets for use in the projection scheme. M_R and M_S contain site indices
		of sites which are identical in structures R (basis) and S (self). N_R and N_S contain all other
		site indices, and N_RS contains pairs of indices in R and S with overlapping augmentation
		spheres in the PAW formalism.

		Arguments:
			basis (Wavefunction object): Wavefunction in the same lattice as self.
				The bands in self will be projected onto the bands in basis
		Returns:
			M_R (numpy array): Indices of sites in basis which have an identical site in
				S (self) (same element and position to within tolerance of 0.02 Angstroms).
			M_S (numpy array): Indices of sites in self which match sites in M_R
				(i.e. M_R[i] is an identical site to M_S[i])
			N_R (numpy array): Indices of sites in basis but not in M_R
			N_S (numpy array): Indices of sites in self but not in M_S
			N_RS (numpy array): Pairs of indices (one in basis and one in self) which
				are not identical but have overlapping augmentation regions
		"""
		ref_sites = basis.structure.sites
		sites = self.structure.sites
		M_R = []
		M_S = []
		for i in range(len(ref_sites)):
			for j in range(len(sites)):
				if ref_sites[i].distance(sites[j]) <= 0.02 and el(ref_sites[i]) == el(sites[j]):
					M_R.append(i)
					M_S.append(j)
		N_R = []
		N_S = []
		for i in range(len(ref_sites)):
			if not i in M_R:
				N_R.append(i)
			for j in range(len(sites)):
				if (not j in N_S) and (not j in M_S):
					N_S.append(j)
		N_RS = []
		for i in N_R:
			for j in N_S:
				if ref_sites[i].distance(sites[j]) < self.cr.pps[el(ref_sites[i])].rmax + self.cr.pps[el(sites[j])].rmax:
					N_RS.append((i,j))
		return M_R, M_S, N_R, N_S, N_RS


	def full_projection(self, basis):
		#self.make_site_lists(basis)
		res = (c_double * 2)()
		#self.projector.full_pseudoprojection(basis, self.pwf, res)
		return np.array((res[0], res[1]))

	def single_band_projection(self, band_num, basis):
		res = self.projector.pseudoprojection(basis.pwf.wf_ptr, self.pwf.wf_ptr, band_num)
		nband = self.projector.get_nband(basis.pwf.wf_ptr)
		nwk = self.projector.get_nwk(basis.pwf.wf_ptr)
		nspin = self.projector.get_nspin(basis.pwf.wf_ptr)
		print("datsa", nband, nwk, nspin)
		res = cdouble_to_numpy(res, 2*nband*nwk*nspin)
		M_R, M_S, N_R, N_S, N_RS = self.make_site_lists(basis)
		projector_list, selfnums, selfcoords, basisnums, basiscoords = self.make_c_projectors(basis)
		ct = self.projector.compensation_terms(band_num, self.pwf.wf_ptr, basis.pwf.wf_ptr, projector_list, 
			len(self.cr.pps), len(M_R), len(N_R), len(N_S), len(M_S), numpy_to_cint(M_R),
			numpy_to_cint(N_R), numpy_to_cint(N_S), numpy_to_cint(M_S), numpy_to_cint(selfnums),
			numpy_to_cdouble(selfcoords), numpy_to_cint(basisnums), numpy_to_cdouble(basiscoords),
			numpy_to_cint(self.dim))
		ct = cdouble_to_numpy(ct, 2*nband*nwk*nspin)
		occs = cdouble_to_numpy(self.projector.get_occs(basis.pwf.wf_ptr), nband*nwk*nspin)
		c, v = 0, 0
		for i in range(nband*nwk*nspin):
			temp = (ct[2*i] + res[2*i]) + 1j * (ct[2*i+1] + res[2*i+1])
			#temp = (ct[2*i]) + 1j * (ct[2*i+1])
			if occs[i] > 0.5:
				v += np.absolute(temp) ** 2 * self.pwf.kws[i%nwk] / nspin
			else:
				c += np.absolute(temp) ** 2 * self.pwf.kws[i%nwk] / nspin
		print (res.shape, res)
		print (ct.shape, ct)
		print ('c, v', c, v)
		self.projector_list = projector_list

	def make_c_projectors(self, basis=None):
		"""
		Uses the CoreRegion objects in self and basis (if not None)
		to construct C representations of the projectors and partial waves
		for a structure. Also assigns numerical labels for each element and
		returns a list of indices and positions which can be easily converted
		to C lists for projection functions.

		Arguments:
			basis (None or Wavefunction): an additional structure from which
				to include pseudopotentials. E.g. can be useful if a basis contains
				some different elements than self.
		Returns:
			projector_list (C pointer): describes the pseudopotential data in C
			selfnums (int32 numpy array): numerical element label for each site in
				the structure
			selfcoords (float64 numpy array): flattened list of coordinates of each site
				in self
			basisnums (if basis != None): same as selfnums, but for basis
			basiscoords (if basis != None): same as selfcoords, but for basis
		"""
		pps = {}
		labels = {}
		label = 0
		for e in self.cr.pps:
			pps[label] = self.cr.pps[e]
			labels[e] = label
			label += 1
		print (pps, labels)
		if basis != None:
			for e in basis.cr.pps:
				if not e in labels:
					pps[label] = basis.cr.pps[e]
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
		print ("rmax", self.cr.pps['Ga'].rmax * 0.529177)

		projector_list = self.projector.get_projector_list(num_els, numpy_to_cint(clabels),
			numpy_to_cint(ls), numpy_to_cdouble(pgrids), numpy_to_cdouble(wgrids),
			numpy_to_cdouble(projectors), numpy_to_cdouble(aewaves), numpy_to_cdouble(pswaves),
			numpy_to_cdouble(np.array([2.10244619055465565])))
		selfnums = np.array([labels[el(s)] for s in self.structure], dtype=np.int32)
		basisnums = np.array([labels[el(s)] for s in basis.structure], dtype=np.int32)
		selfcoords = np.array([], np.float64)
		basiscoords = np.array([], np.float64)

		for s in self.structure:
			selfcoords = np.append(selfcoords, s.frac_coords)
		if basis != None:
			for s in basis.structure:
				basiscoords = np.append(basiscoords, s.frac_coords)
			return projector_list, selfnums, selfcoords, basisnums, basiscoords
		return projector_list, selfnums, selfcoords

	def proportion_conduction(self, band_num, bulk):
		pass

	def free_all(self):
		self.projector.free_pswf(self.pwf.wf_ptr)
		if self.projector_list != None:
			self.projector.free_ppot_list(self.projector_list, len(self.cr.pps))

posb = Poscar.from_file("CONTCAR").structure
posd = Poscar.from_file("CONTCAR").structure
pot = Potcar.from_file("POTCAR")
pwf1 = PseudoWavefunction("WAVECAR", "vasprun.xml")
pwf2 = PseudoWavefunction("WAVECAR", "vasprun.xml")

wf1 = Wavefunction(posb, pwf1, CoreRegion(pot), (40,40,40))
wf2 = Wavefunction(posd, pwf2, CoreRegion(pot), (40,40,40))
for i in range(0,1):
	wf2.single_band_projection(i, wf1)

wf1.free_all()
wf2.free_all()

#For each structure
#numerical element label for each site
#Fractional coordinates of sites
#grid for each element
#Labels (element, l, ndata, length) for each projector-partial set for each element
#projector functions for each label above for each element
#AE partial waves for each label above for each element
#PS partial waves for each label above for each element
