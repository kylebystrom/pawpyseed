from pymatgen.io.vasp.inputs import Potcar
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.core.structure import Structure
import numpy as np
from ctypes import *
import os
import numpy as np
import json

MODULE_DIR = os.path.dirname(os.path.abspath(__file__))

def cdouble_to_numpy(arr, length):
	arr = cast(arr, POINTER(c_double))
	newarr = np.zeros(length)
	for i in range(length):
		newarr[i] = arr[i]
	return newarr

def cfloat_to_numpy(arr, length):
	arr = cast(arr, POINTER(c_float))
	newarr = np.zeros(length)
	for i in range(length):
		newarr[i] = arr[i]
	return newarr

def numpy_to_cdouble(arr):
	newarr = (c_double * len(weights))()
	for i in range(len(arr)):
		newarr[i] = arr[i]
	return newarr

def numpy_to_cfloat(arr):
	newarr = (c_float * len(weights))()
	for i in range(len(arr)):
		newarr[i] = arr[i]
	return newarr

def el(site):
	return site.specie.element.symbol

class Pseudopotential:

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
		self.grid = self.make_nums(gridstr)
		self.aepotential = self.make_nums(aepotstr)
		self.aecorecharge = self.make_nums(corechgstr)
		self.kinetic = self.make_nums(kenstr)
		self.pspotential = self.make_nums(pspotstr)
		self.pscorecharge = self.make_nums(pscorechgstr)

		for pwave in partial_waves:
			lst = pwave.split("ae wavefunction", 1)
			self.pswaves.append(self.make_nums(lst[0]))
			self.aewaves.append(self.make_nums(lst[1]))

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

		projgridstr = settingstr.split("STEP   =")[-1]
		projgridstr = projgridstr.split("END")[0]
		self.projgrid = self.make_nums(projgridstr)
		self.step = (self.projgrid[0], self.projgrid[1])
		self.projgrid = self.projgrid[2:]

		self.projgrid = np.linspace(0,rmax/1.88973,100,False)
		y = self.realprojs[1]

	def make_nums(self, numstring):
		return np.array([float(num) for num in numstring.split()])

class CoreRegion:

	def __init__(self, potcar):
		self.pps = {}
		for potsingle in potcar:
			self.pps[potsingle.element] = Pseudopotential(potsingle.data[:-15], potsingle.rmax)
		

class PseudoWavefunction:

	def __init__(self, filename="WAVECAR", vr="vasprun.xml"):
		self.reader = CDLL(os.path.join(MODULE_DIR, "pawpy.so"))
		if type(vr) == str:
			vr = Vasprun(vr)
		weights = vr.actual_kpoints_weights
		kws = (c_double * len(weights))()
		for i in range(len(weights)):
			kws[i] = weights[i]
		self.kws = weights
		self.wf_ptr = self.reader.read_wavefunctions(filename, byref(kws))

class Wavefunction:

	def __init__(self, struct, pwf, cr):
		self.structure = struct
		self.pwf = pwf
		self.cr = cr
		self.projector = pwf.reader

	def from_files(self, pwf="WAVECAR", cr="POTCAR", vr="vasprun.xml"):
		return Wavefunction(PseudoWavefunction(pwf, vr), Potcar.from_file(cr))

	def make_site_lists(self, basis):
		ref_sites = basis.structure.sites
		sites = self.structure.sites
		M_R = []
		M_S = []
		for i in range(len(ref_sites)):
			for j in range(len(sites)):
				if ref_sites[i].distance(sites[j]) <= 0.02:
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
				if ref_sites[i].distance(sites[j]) < self.cr.pps(el(ref_sites[i])).rmax + self.cr.pps(el(sites[j])).rmax:
					N_RS.append((i,j))
		return zip(M_R, M_S), N_R, N_S, N_RS


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
		res = cfloat_to_numpy(res, 2*nband*nwk*nspin)
		print res
		re = res[::2]
		im = res[1::2]
		mag = re**2 + im**2
		for i in range(len(mag)):
			mag[i] *= self.pwf.kws[i%len(self.pwf.kws)]
		print sum(mag[:256*32]), sum(mag[256*32:])
		#pass

	def proportion_conduction(self, band_num, bulk):
		pass


pwf1 = PseudoWavefunction("bulk/WAVECAR", "bulk/vasprun.xml")
pwf2 = PseudoWavefunction("charge_0/WAVECAR", "charge_0/vasprun.xml")

wf1 = Wavefunction(None, pwf1, None)
wf2 = Wavefunction(None, pwf2, None)
for i in range(253,257):
	wf2.single_band_projection(i, wf1)
