# coding: utf-8

## @package pawpyseed.core.wavefunction
# Base class containing Python classes for parsing files
# and storing and analyzing wavefunction data.

from pymatgen.io.vasp.inputs import Potcar, Poscar
from pymatgen.io.vasp.outputs import Vasprun, Outcar
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.core.operations import SymmOp
import numpy as np
from ctypes import *
from pawpyseed.core.utils import *
import os, time
import numpy as np
import json
from monty.io import zopen

import sys

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

class Pseudopotential:
	"""
	Contains important attributes from a VASP pseudopotential files. POTCAR
	"settings" can be read from the pymatgen POTCAR object

	Note: for the following attributes, 'index' refers to an energy
	quantum number epsilon and angular momentum quantum number l,
	which define one set consisting of a projector function, all electron
	partial waves, and pseudo partial waves.

	Attributes:
		rmaxstr (str): Maximum radius of the projection operators, as string
			of double precision float
		grid (np.array): radial grid on which partial waves are defined
		aepotential (np.array): All electron potential defined radially on grid
		aecorecharge (np.array): All electron core charge defined radially
			on grid (i.e. charge due to core, and not valence, electrons)
		kinetic (np.array): Core kinetic energy density, defined raidally on grid
		pspotential (np.array): pseudopotential defined on grid
		pscorecharge (np.array): pseudo core charge defined on grid
		ls (list): l quantum number for each index
		pswaves (list of np.array): pseudo partial waves for each index
		aewaves (list of np.array): all electron partial waves for each index
		projgrid (np.array): radial grid on which projector functions are defined
		recipprojs (list of np.array): reciprocal space projection operators
			for each index
		realprojs (list of np.array): real space projection operators
			for each index
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
		self.rmaxstrs = []

		auguccstr, gridstr = gridstr.split("grid", 1)
		gridstr, aepotstr = gridstr.split("aepotential", 1)
		aepotstr, corechgstr = aepotstr.split("core charge-density", 1)
		try:
			corechgstr, kenstr = corechgstr.split("kinetic energy-density", 1)
			kenstr, pspotstr = kenstr.split("pspotential", 1)
		except:
			kenstr = "0 0"
			corechgstr, pspotstr = corechgstr.split("pspotential", 1)
		pspotstr, pscorechgstr = pspotstr.split("core charge-density (pseudized)", 1)
		self.grid = self.make_nums(gridstr)
		self.aepotential = self.make_nums(aepotstr)
		self.aecorecharge = self.make_nums(corechgstr)
		self.kinetic = self.make_nums(kenstr)
		self.pspotential = self.make_nums(pspotstr)
		self.pscorecharge = self.make_nums(pscorechgstr)

		augstr, uccstr = auguccstr.split('uccopancies in atom', 1)
		head, augstr = augstr.split('augmentation charges (non sperical)', 1)
		self.augs = self.make_nums(augstr)

		for pwave in partial_waves:
			lst = pwave.split("ae wavefunction", 1)
			self.pswaves.append(self.make_nums(lst[0]))
			self.aewaves.append(self.make_nums(lst[1]))

		projstrs = nonradial.split("Non local Part")
		topstr, projstrs = projstrs[0], projstrs[1:]
		self.T = float(topstr[-22:-4])
		topstr, atpschgstr = topstr[:-22].split("atomic pseudo charge-density", 1)
		try:
			topstr, corechgstr = topstr.split("core charge-density (partial)", 1)
			settingstr, localstr = topstr.split("local part", 1)
		except:
			corechgstr = "0 0"
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
			self.rmaxstr = c_char_p()
			self.rmaxstr.value = nonlocalvals.split()[2].encode('utf-8')
			self.rmax = self.make_nums(nonlocalvals.split()[2])[0]
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
		return np.fromstring(numstring, dtype = np.float64, sep = ' ')

class CoreRegion:
	"""
	List of Pseudopotential objects to describe the core region of a structure.

	Attributes:
		pps (dict of Pseudopotential): keys are element symbols,
			values are Pseudopotential objects
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

	Attributes:
		kpts (np.array): nx3 array of fractional kpoint vectors,
			where n is the number of kpoints
		kws (np.array): weight of each kpoint
		wf_ptr (ctypes POINTER): c pointer to pswf_t object
	"""

	def __init__(self, filename="WAVECAR", vr="vasprun.xml"):
		if type(vr) == str:
			vr = Vasprun(vr)
		weights = vr.actual_kpoints_weights
		kws = numpy_to_cdouble(weights)
		self.kws = np.array(weights, dtype = np.float64)
		self.kpts = np.array(vr.actual_kpoints, dtype=np.float64)
		if '.gz' in filename or '.bz2' in filename:
			f = zopen(filename, 'rb')
			contents = f.read()
			f.close()
			self.wf_ptr = PAWC.read_wavefunctions_from_str(
				contents, kws)
		else:
			self.wf_ptr = PAWC.read_wavefunctions(filename.encode('utf-8'), kws)
		self.ncl = PAWC.is_ncl(self.wf_ptr) > 0

	def pseudoprojection(self, band_num, basis):
		"""
		Computes <psibt_n1k|psit_n2k> for all n1 and k
		and a given n2, where psibt are basis structures
		pseudowavefunctions and psit are self pseudowavefunctions

		Arguments:
			band_num (int): n2 (see description)
			basis (Pseudowavefunction): pseudowavefunctions onto whose bands
				the band of self is projected
		"""
		nband = PAWC.get_nband(c_void_p(basis.wf_ptr))
		nwk = PAWC.get_nwk(c_void_p(basis.wf_ptr))
		nspin = PAWC.get_nspin(c_void_p(basis.wf_ptr))

		res = PAWC.pseudoprojection(c_void_p(basis.wf_ptr), c_void_p(self.wf_ptr), band_num)
		res = cdouble_to_numpy(res, 2*nband*nwk*nspin)
		return res[::2] + 1j * res[1::2]


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
		wf_ptr (C pointer): pointer to the pswf_t C object for this wavefunction
		dim (np.ndarray, length 3): dimension of the FFT grid used by VASP
			and therefore for FFTs in this code
	"""

	def __init__(self, struct, pwf, cr, outcar, setup_projectors=False):
		"""
		Arguments:
			struct (pymatgen.core.Structure): structure that the wavefunction describes
			pwf (PseudoWavefunction): Pseudowavefunction component
			cr (CoreRegion): Contains the pseudopotentials, with projectors and
				partials waves, for the structure
			outcar (pymatgen.io.vasp.outputs.Outcar): Outcar object for reading ngf
			setup_projectors (bool, False): Whether to set up the core region
				components of the wavefunctions (leave as False if passing this
				object to Projector, which will do the setup automatically)
		Returns:
			Wavefunction object
		"""

		if pwf.ncl:
			raise PAWpyError("Pseudowavefunction is noncollinear! Call NCLWavefunction(...) instead")
		self.structure = struct
		self.pwf = pwf
		self.cr = cr
		if type(outcar) == Outcar:
			self.dim = outcar.ngf
			self.dim = np.array(self.dim).astype(np.int32) // 2
		else:
			#assume outcar is actually ngf, will fix later
			self.dim = outcar
			self.dim = np.array(self.dim).astype(np.int32)
		self.projector_owner = False
		self.projector_list = None
		self.nums = None
		self.coords = None
		self.nband = PAWC.get_nband(c_void_p(pwf.wf_ptr))
		self.nwk = PAWC.get_nwk(c_void_p(pwf.wf_ptr))
		self.nspin = PAWC.get_nspin(c_void_p(pwf.wf_ptr))
		self.encut = PAWC.get_encut(c_void_p(pwf.wf_ptr))
		if setup_projectors:
			self.check_c_projectors()
		self.num_proj_els = None
		self.freed = False

	@staticmethod
	def from_files(struct="CONTCAR", pwf="WAVECAR", cr="POTCAR",
		vr="vasprun.xml", outcar="OUTCAR", setup_projectors=False):
		"""
		Construct a Wavefunction object from file paths.

		Arguments:
			struct (str): VASP POSCAR or CONTCAR file path
			pwf (str): VASP WAVECAR file path
			cr (str): VASP POTCAR file path
			vr (str): VASP vasprun file path
			outcar (str): VASP OUTCAR file path

		Returns:
			Wavefunction object
		"""
		return Wavefunction(Poscar.from_file(struct).structure,
			PseudoWavefunction(pwf, vr),
			CoreRegion(Potcar.from_file(cr)),
			Outcar(outcar), setup_projectors)

	@staticmethod
	def from_directory(path, setup_projectors = False):
		"""
		Assumes VASP output has the default filenames and is located
		in the directory specificed by path.

		Arguments:
			path (str): VASP output directory
			setup_projectors (bool, False): Whether to set up the core region
				components of the wavefunctions (leave as False if passing this
				object to Projector, which will do the setup automatically)

		Returns:
			Wavefunction object
		"""
		filepaths = []
		for d in ["CONTCAR", "WAVECAR", "POTCAR", "vasprun.xml", "OUTCAR"]:
			filepaths.append(str(os.path.join(path, d)))
		args = filepaths + [setup_projectors]
		return Wavefunction.from_files(*args)

	@staticmethod
	def from_atomate_directory(path, setup_projectors = False):
		"""
		Assumes VASP output has the default filenames and is located
		in the directory specificed by path. Checks for
		gzipped files created by atomate

		Arguments:
			path (str): VASP output directory
			setup_projectors (bool, False): Whether to set up the core region
				components of the wavefunctions (leave as False if passing this
				object to Projector, which will do the setup automatically)

		Returns:
			Wavefunction object
		"""

		files = ["CONTCAR", "WAVECAR", "POTCAR", "vasprun.xml", "OUTCAR"]
		paths = []

		for file in files:
		    filepat = os.path.join( path, file +'.relax2.gz')
		    if not os.path.exists( filepat):
		        filepat = os.path.join( path, file +'.relax1.gz')
		    if not os.path.exists( filepat):
		        filepat = os.path.join( path, file +'.gz')
		    if not os.path.exists( filepat):
		        filepat = os.path.join( path, file)
		    if not os.path.exists( filepat):
		        print('Could not find {}! Skipping this defect...'.format(file))
		        return False

		    paths.append(filepat)

		args = paths + [setup_projectors]
		wf = Wavefunction.from_files(*args)

		return wf

	def get_c_projectors_from_pps(self, pps):
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
		clabels = np.array([], np.int32)
		ls = np.array([], np.int32)
		projectors = np.array([], np.float64)
		aewaves = np.array([], np.float64)
		pswaves = np.array([], np.float64)
		wgrids = np.array([], np.float64)
		pgrids = np.array([], np.float64)
		augs = np.array([], np.float64)
		rmaxstrs = (c_char_p * len(pps))()
		rmaxs = np.array([], np.float64)
		num_els = 0

		for num in sorted(pps.keys()):
			pp = pps[num]
			clabels = np.append(clabels, [num, len(pp.ls), pp.ndata, len(pp.grid)])
			rmaxstrs[num_els] = pp.rmaxstr
			rmaxs = np.append(rmaxs, pp.rmax)
			ls = np.append(ls, pp.ls)
			wgrids = np.append(wgrids, pp.grid)
			pgrids = np.append(pgrids, pp.projgrid)
			augs = np.append(augs, pp.augs)
			num_els += 1
			for i in range(len(pp.ls)):
				proj = pp.realprojs[i]
				aepw = pp.aewaves[i]
				pspw = pp.pswaves[i]
				projectors = np.append(projectors, proj)
				aewaves = np.append(aewaves, aepw)
				pswaves = np.append(pswaves, pspw)

		#print (num_els, clabels, ls, pgrids, wgrids, rmaxs)
		grid_encut = (2 * np.pi * self.dim / self.structure.lattice.abc)**2 / 0.262
		projector_list = cfunc_call(PAWC.get_projector_list, None,
							num_els, clabels, ls, pgrids, wgrids,
							projectors, aewaves, pswaves,
							rmaxs, max(grid_encut))
		end = time.monotonic()
		print('--------------\nran get_projector_list in %f seconds\n---------------' % (end-start))
		return projector_list
			

	def make_c_projectors(self):
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
		
		projector_list = self.get_c_projectors_from_pps(pps)

		selfnums = np.array([labels[el(s)] for s in self.structure], dtype=np.int32)
		selfcoords = np.array([], np.float64)

		self.num_proj_els = len(pps)
		for s in self.structure:
			selfcoords = np.append(selfcoords, s.frac_coords)
		return projector_list, selfnums, selfcoords

	def check_c_projectors(self):
		"""
		Check to see if the projector functions have been read in and set up.
		If not, do so.
		"""
		if not self.projector_list:
			start = time.monotonic()
			self.projector_owner = True
			self.projector_list, self.nums, self.coords = self.make_c_projectors()
			cfunc_call(PAWC.setup_projections, None,
					self.pwf.wf_ptr, self.projector_list,
					self.num_proj_els, len(self.structure), self.dim,
					self.nums, self.coords)
			end = time.monotonic()
			print('--------------\nran setup_projections in %f seconds\n---------------' % (end-start))

	def get_state_realspace(self, b, k, s, dim=None):
		"""
		Returns the real and imaginary parts of a given band.
		Args:
			b (int): band number
			k (int): kpoint number
			s (int): spin number
			dim (numpy array of 3 ints): dimensions of the FFT grid
		Returns:
			An array (x slow-indexed) where the first half of the values
				are the real part and second half of the values are the
				imaginary part
		"""

		self.check_c_projectors()
		if type(dim) == type(None):
			dim = self.dim
		return cfunc_call(PAWC.realspace_state_ri, 2*dim[0]*dim[1]*dim[2], b, k+s*self.nwk,
			self.pwf.wf_ptr, self.projector_list,
			dim, self.nums, self.coords)

	def _convert_to_vasp_volumetric(self, filename, dim):
		

		#from pymatgen VolumetricData class
		p = Poscar(self.structure)
		lines = filename + '\n'
		lines += "   1.00000000000000\n"
		latt = self.structure.lattice.matrix
		lines += " %12.6f%12.6f%12.6f\n" % tuple(latt[0, :])
		lines += " %12.6f%12.6f%12.6f\n" % tuple(latt[1, :])
		lines += " %12.6f%12.6f%12.6f\n" % tuple(latt[2, :])
		lines += "".join(["%5s" % s for s in p.site_symbols]) + "\n"
		lines += "".join(["%6d" % x for x in p.natoms]) + "\n"
		lines += "Direct\n"
		for site in self.structure:
			lines += "%10.6f%10.6f%10.6f\n" % tuple(site.frac_coords)
		lines += " \n"
	
		f = open(filename, 'r')
		nums = f.read()
		f.close()
		f = open(filename, 'w')
		dimstr = '%d %d %d\n' % (dim[0], dim[1], dim[2])
		#pos = Poscar(self.structure, velocities = None)
		#posstr = pos.get_string() + '\n'
		f.write(lines + dimstr + nums)
		f.close()

	def write_state_realspace(self, b, k, s, fileprefix = "", dim=None, return_wf = False):
		"""
		Writes the real and imaginary parts of a given band to two files,
		prefixed by fileprefix

		Args:
			b (int): band number (0-indexed!)
			k (int): kpoint number (0-indexed!)
			s (int): spin number (0-indexed!)
			dim (numpy array of 3 ints): dimensions of the FFT grid
			fileprefix (string, optional): first part of the file name
			return_wf (bool): whether to return the wavefunction
		Returns:
			(if return_wf==True) An array (x slow-indexed) where the first half of the values
				are the real part and second half of the values are the
				imaginary part
			The wavefunction is written with z the slow index.
		"""
		res = None
		print("PARAMETERS", self.nums, self.coords, dim)
		sys.stdout.flush()
		self.check_c_projectors()
		if type(dim) == type(None):
			dim = self.dim
		filename_base = "%sB%dK%dS%d" % (fileprefix, b, k, s)
		filename1 = "%s_REAL" % filename_base
		filename2 = "%s_IMAG" % filename_base
		print("PARAMETERS", self.nums, self.coords, dim)
		sys.stdout.flush()
		if return_wf:
			res = cfunc_call(PAWC.write_realspace_state_ri_return, 2*dim[0]*dim[1]*dim[2],
				filename1, filename2,
				b, k+s*self.nwk,
				self.pwf.wf_ptr, self.projector_list,
				dim, self.nums, self.coords)
		else:
			cfunc_call(PAWC.write_realspace_state_ri_noreturn, None, filename1, filename2,
				b, k+s*self.nwk,
				self.pwf.wf_ptr, self.projector_list,
				dim, self.nums, self.coords)
		self._convert_to_vasp_volumetric(filename1, dim)
		self._convert_to_vasp_volumetric(filename2, dim)
		return res

	def write_density_realspace(self, filename = "PYAECCAR", dim=None, return_wf = False):
		"""
		Writes the AE charge density to a file. Returns it if desired.

		Args:
			b (int): band number
			k (int): kpoint number
			s (int): spin number
			dim (numpy array of 3 ints): dimensions of the FFT grid
			filename (string, "PYAECCAR"): charge density filename
			return_wf (bool): whether to return the wavefunction
		Returns:
			(if return_wf==True) An array (x slow-indexed, as in VASP)
				with the charge densities
			The charge density is written with z the slow index.
		"""

		res = None
		self.check_c_projectors()
		if type(dim) == type(None):
			dim = self.dim
		if return_wf:
			res = cfunc_call(PAWC.write_density_return, dim[0]*dim[1]*dim[2], filename,
				self.pwf.wf_ptr, self.projector_list, dim, self.nums, self.coords)
		else:
			cfunc_call(PAWC.write_density_noreturn, None, filename,
				self.pwf.wf_ptr, self.projector_list, dim, self.nums, self.coords)
			res = None
		self._convert_to_vasp_volumetric(filename, dim)
		return res

	def get_symmops(self, symprec):
		sga = SpacegroupAnalyzer(self.structure, symprec)
		symmops = sga.get_symmetry_operations(cartesian = True)
		lattice = self.structure.lattice.matrix
		invlattice = self.structure.lattice.inv_matrix
		newops = []
		for op in symmops:
			newrot = np.dot(lattice, op.rotation_matrix)
			newrot = np.dot(newrot, invlattice)
			newtrans = np.dot(op.translation_vector, invlattice)
			newops.append(SymmOp.from_rotation_and_translation(
				newrot, newtrans))
		return newops

	def get_nosym_kpoints(self, init_kpts = None, symprec=1e-5,
		gen_trsym = True, fil_trsym = True):

		kpts = np.array(self.pwf.kpts)
		allkpts = [] if init_kpts == None else [kpt for kpt in init_kpts]
		orig_kptnums = []
		op_nums = []
		symmops = self.get_symmops(symprec)
		trs = []
		for i, op in enumerate(symmops):
			for k, kpt in enumerate(kpts):
				newkpt = np.dot(op.rotation_matrix, kpt)
				newkpt -= np.around(newkpt)
				newkpt[ abs(newkpt + 0.5) < 1e-5 ] = 0.5
				#if ((newkpt > 0.5+1e-6) + (newkpt < -0.5+1e-6)).any():
				#	continue
				if fil_trsym:
					if newkpt[2] < -1e-6 or \
						(abs(newkpt[2]) < 1e-6 and newkpt[1] < -1e-6) or \
						(abs(newkpt[2]) < 1e-6 and abs(newkpt[1]) < 1e-6 and newkpt[0] < -1e-6):
						continue
				unique = True
				for nkpt in allkpts:
					diff = (newkpt - nkpt) % 1
					oppdiff = 1 - diff
					tst = (np.abs(diff) < 1e-4) + (np.abs(oppdiff) < 1e-4)
					if ( tst.all() ):
						unique = False
						break
				if unique:
					allkpts.append(newkpt)
					orig_kptnums.append(k)
					op_nums.append(i)
					trs.append(0)
		if gen_trsym:
			for i, op in enumerate(symmops):
				for k, kpt in enumerate(kpts):
					newkpt = np.dot(op.rotation_matrix, kpt) * -1
					newkpt -= np.around(newkpt)
					newkpt[ abs(newkpt + 0.5) < 1e-5 ] = 0.5
					#if ((newkpt > 0.5+1e-6) + (newkpt < -0.5+1e-6)).any():
					#	continue
					if fil_trsym:
						if newkpt[2] < -1e-10 or \
							(abs(newkpt[2]) < 1e-6 and newkpt[1] < -1e-6) or \
							(abs(newkpt[2]) < 1e-6 and abs(newkpt[1]) < 1e-6 and newkpt[0] < -1e-6):
							continue
					unique = True
					for nkpt in allkpts:
						diff = (newkpt - nkpt) % 1
						oppdiff = 1 - diff
						tst = (np.abs(diff) < 1e-4) + (np.abs(oppdiff) < 1e-4)
						if ( tst.all() ):
							unique = False
							break
					if unique:
						allkpts.append(newkpt)
						orig_kptnums.append(k)
						op_nums.append(i)
						trs.append(1)
		self.nosym_kpts = allkpts
		self.orig_kptnums = orig_kptnums
		self.op_nums = op_nums
		self.symmops = symmops
		return np.array(allkpts), orig_kptnums, op_nums, symmops, trs

	def get_kpt_mapping(self, allkpts, symprec=1e-5, gen_trsym = True):
		symmops = self.get_symmops(symprec)
		kpts = np.array(self.pwf.kpts)
		orig_kptnums = []
		op_nums = []
		trs = []
		for nkpt in allkpts:
			match = False
			for i, op in enumerate(symmops):
				for k, kpt in enumerate(kpts):
					newkpt = np.dot(op.rotation_matrix, kpt)
					#if ((newkpt > 0.5+1e-6) + (newkpt < -0.5+1e-6)).any():
					#	continue
					diff = (newkpt - nkpt) % 1
					oppdiff = 1 - diff
					tst = (np.abs(diff) < 1e-4) + (np.abs(oppdiff) < 1e-4)
					if tst.all():
						match = True
						orig_kptnums.append(k)
						op_nums.append(i)
						trs.append(0)
						break
				if match:
					break
			if match:
				continue
			for i, op in enumerate(symmops):
				for k, kpt in enumerate(kpts):
					newkpt = np.dot(op.rotation_matrix, kpt) * -1
					#if ((newkpt > 0.5+1e-6) + (newkpt < -0.5+1e-6)).any():
					#	continue
					diff = (newkpt - nkpt) % 1
					oppdiff = 1 - diff
					tst = (np.abs(diff) < 1e-4) + (np.abs(oppdiff) < 1e-4)
					if tst.all():
						match = True
						orig_kptnums.append(k)
						op_nums.append(i)
						trs.append(1)
						break
				if match:
					break
			if not match:
				raise PAWpyError("Could not find kpoint mapping to %s" % str(nkpt))
		return orig_kptnums, op_nums, symmops, trs

	@property
	def kpts(self):
		return self.pwf.kpts

	@property
	def kws(self):
		return self.pwf.kws

	def free_all(self):
		"""
		Frees all of the C structures associated with the Wavefunction object.
		After being called, this object is not usable.
		"""
		PAWC.free_pswf(c_void_p(self.pwf.wf_ptr))
		if self.projector_owner:
			PAWC.free_ppot_list(c_void_p(self.projector_list), len(self.cr.pps))
		self.freed = True

