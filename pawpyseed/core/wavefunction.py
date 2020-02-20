# coding: utf-8

## @package pawpyseed.core.wavefunction
# Base class containing Python classes for parsing files
# and storing and analyzing wavefunction data.

from pymatgen.io.vasp.inputs import Potcar, Poscar
from pymatgen.io.vasp.outputs import Vasprun, Outcar
from pymatgen.core.structure import Structure
import numpy as np
from pawpyseed.core.utils import *
import pawpyseed.core.symmetry as pawpy_symm
import os, time
import numpy as np
import json

import sys

from pawpyseed.core import pawpyc

class Pseudopotential:
	"""
	Contains important attributes from a VASP pseudopotential files. POTCAR
	"settings" can be read from the pymatgen POTCAR object

	If you use pymatgen, you can think of this as correlating with
	the PotcarSingle object.

	Note: for the following attributes, 'index' refers to an energy
	quantum number epsilon and angular momentum quantum number l,
	which define one set consisting of a projector function, all electron
	partial waves, and pseudo partial waves.

	Attributes:
		rmax (np.float64): Maximum radius of the projection operators
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

	def __init__(self, data):
		"""
		Initializer for Pseudopotential.
		Should only be used by CoreRegion.

		Arguments:
			data (str): single-element pseudopotential
				(POTCAR) as a string
		"""
		nonradial, radial = data.split("PAW radial sets", 1)
		partial_waves = radial.split("pseudo wavefunction")
		gridstr, partial_waves = partial_waves[0], partial_waves[1:]
		self.pswaves = []
		self.aewaves = []
		self.recipprojs = []
		self.realprojs = []
		self.nonlocalprojs = []
		self.ls = []

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
		#self.aepotential = self.make_nums(aepotstr)
		#self.aecorecharge = self.make_nums(corechgstr)
		#self.kinetic = self.make_nums(kenstr)
		#self.pspotential = self.make_nums(pspotstr)
		#self.pscorecharge = self.make_nums(pscorechgstr)

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
		"""
		if "gradient corrections used for XC" in localstr:
			localstr, self.gradxc = localstr.split("gradient corrections used for XC", 1)
			self.gradxc = int(self.gradxc)
		else:
			self.gradxc = None
		self.localpart = self.make_nums(localstr)
		self.localnum = self.localpart[0]
		self.localpart = self.localpart[1:]
		self.coredensity = self.make_nums(corechgstr)
		self.atomicdensity = self.make_nums(atpschgstr)
		"""

		for projstr in projstrs:
			lst = projstr.split("Reciprocal Space Part")
			nonlocalvals, projs = lst[0], lst[1:]
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
		#projgridstr = projgridstr.split("END")[0]
		self.projgrid = np.arange(len(self.realprojs[0])) * self.rmax / len(self.realprojs[0])
		self.step = (self.projgrid[0], self.projgrid[1])

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
		"""
		Returns a new CoreRegion object from a
		pymatgen.io.vasp.inputs.Potcar class.

		Arguments:
			potcar (pymatgen.io.vasp.inputs.Potcar): Potcar file for
				the VASP calculation

		Returns:
			CoreRegion object based on potcar
		"""
		self.pps = {}
		for potsingle in potcar:
			self.pps[potsingle.element] = Pseudopotential(potsingle.data[:-15])


class Wavefunction(pawpyc.CWavefunction):
	"""
	Class for storing and manipulating all electron wave functions in the PAW
	formalism.

	Attributes:
		structure (pymatgen.core.structure.Structure): stucture of the material
			that the wave function describes
		cr (CoreRegion): Contains the pseudopotentials, with projectors and
			partials waves, for the structure
		dim (np.ndarray, length 3): dimension of the FFT grid used by VASP
			and therefore for FFTs in this code
		band_props (np.ndarray): 4-item array of containing the information
			(band gap, cbm, vbm, is_band_gap_direct). This object contains the same
			information as pymatgen.io.vasp.outputs.Vasprun.eigenvalue_band_properties
	"""

	def __init__(self, struct, pwf, cr, dim, symprec = 1e-4, setup_projectors=False):
		"""
		Arguments:
			struct (pymatgen.core.Structure): structure that the wavefunction describes
			pwf (pawpyc.PWFPointer): holder class for pswf_t and k-points/k-point weights
			cr (CoreRegion): Contains the pseudopotentials, with projectors and
				partials waves, for the structure
			dim (pymatgen.io.vasp.outputs.Outcar OR np.ndarry OR list of length 3):
				Outcar object for reading ngf or the dimensions NG* of the FFT grid
			symprec (float, 1e-4): precision tolerance for symmetry operations
			setup_projectors (bool, False): Whether to set up the core region
				components of the wavefunctions. Pawpyseed will set up the projectors
				automatically when they are first needed, so this generally
				can be left as False.
		Returns:
			Wavefunction object
		"""
		self.band_props = pwf.band_props.copy(order = 'C')
		super(Wavefunction, self).__init__(pwf)
		if self.ncl:
			raise PAWpyError("Pseudowavefunction is noncollinear! Call NCLWavefunction(...) instead")
		self.structure = struct
		self.symprec = symprec
		self.cr = cr
		self.dim = np.array(dim).astype(np.int32)
		if len(dim) != 3:
			raise PAWpyError("Grid dimensions must be length 3")
		if setup_projectors:
			self.check_c_projectors()

	def check_band_index(self, b):
		if b < 0 or b >= self.nband:
			raise ValueError("Invalid band {}. Should be in range [{}, {}]".format(
								b, 0, self.nband-1))

	def check_kpoint_index(self, k):
		if k < 0 or k >= self.nwk:
			raise ValueError("Invalid kpoint index {}. Should be in range [{}, {}]".format(
								k, 0, self.nwk-1))

	def check_spin_index(self, s):
		if s < 0 or s >= self.nspin:
			raise ValueError("Spin must be 0 for non-spin-polarized or 0 or 1 for spin-polarized.")

	def check_bks_spec(self, b, k, s):
		self.check_band_index(b)
		self.check_kpoint_index(k)
		self.check_spin_index(s)

	def update_dim(self, dim):
		self.dim = np.array(dim, dtype=np.int32)
		self.update_dimv(dim)

	def desymmetrized_copy(self, allkpts=None, weights=None, symprec=None,
							time_reversal_symmetry=True):
		"""
		Returns a copy of self with a k-point mesh that is not reduced
		using crystal symmetry.

		Arguments:
			allkpts (optional, None): An optional k-point mesh to map
				onto. Used by the Projector class for some cases
			weights (optional, None): If allkpts is not None, weights
				should contain the k-point weights of each k-point,
				with the sum normalized to 1.
			symprec: Symmetry precision to use when determining the space group.
				If None, the symmetry precision used to generate the
				Wavefunction will be used (the default).
			time_reversal_symmetry: Whether time reversal symmetry is used.
		"""
		if not symprec:
			symprec = self.symprec

		pwf = self._desymmetrized_pwf(self.structure, self.band_props, allkpts, weights,
										symprec, time_reversal_symmetry)
		new_wf = Wavefunction(self.structure, pwf, self.cr, self.dim, symprec=symprec)
		return new_wf

	@staticmethod
	def from_files(struct="CONTCAR", wavecar="WAVECAR", cr="POTCAR",
		vr="vasprun.xml", setup_projectors=False):
		"""
		Construct a Wavefunction object from file paths.

		Arguments:
			struct (str): VASP POSCAR or CONTCAR file path
			wavecar (str): VASP WAVECAR file path
			cr (str): VASP POTCAR file path
			vr (str): VASP vasprun file path
			outcar (str): VASP OUTCAR file path
			setup_projectors (bool, False): Whether to set up the core region
				components of the wavefunctions. Pawpyseed will set up the projectors
				automatically when they are first needed, so this generally
				can be left as False.

		Returns:
			Wavefunction object
		"""
		for fname in [struct, wavecar, cr, vr]:
			if not os.path.isfile(fname):
				raise FileNotFoundError("File {} does not exist.".format(fname))
		vr = Vasprun(vr)
		dim = np.array([vr.parameters["NGX"], vr.parameters["NGY"], vr.parameters["NGZ"]])
		symprec = vr.parameters["SYMPREC"]
		pwf = pawpyc.PWFPointer(wavecar, vr)
		return Wavefunction(Poscar.from_file(struct).structure,
			pwf, CoreRegion(Potcar.from_file(cr)),
			dim, symprec, setup_projectors)

	@staticmethod
	def from_directory(path, setup_projectors = False):
		"""
		Assumes VASP output has the default filenames and is located
		in the directory specificed by path.

		Arguments:
			path (str): VASP output directory
			setup_projectors (bool, False): Whether to set up the core region
				components of the wavefunctions. Pawpyseed will set up the projectors
				automatically when they are first needed, so this generally
				can be left as False.

		Returns:
			Wavefunction object
		"""
		filepaths = []
		for d in ["CONTCAR", "WAVECAR", "POTCAR", "vasprun.xml"]:
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
				components of the wavefunctions. Pawpyseed will set up the projectors
				automatically when they are first needed, so this generally
				can be left as False.

		Returns:
			Wavefunction object
		"""

		files = ["CONTCAR", "WAVECAR", "POTCAR", "vasprun.xml"]
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

	def _make_c_projectors(self):
		"""
		Uses the CoreRegion objects in self
		to construct C representations of the projectors and partial waves
		for a structure. Also assigns numerical labels for each element and
		setups up a list of indices and positions which can be easily converted
		to C lists for projection routines.
		"""

		pps = {}
		labels = {}
		label = 0
		for e in self.cr.pps:
			pps[label] = self.cr.pps[e]
			labels[e] = label
			label += 1

		nums = np.array([labels[el(s)] for s in self.structure], dtype=np.int32)
		coords = np.array([], dtype = np.float64)

		self.num_sites = len(self.structure)
		self.num_elems = len(pps)
		for s in self.structure:
			coords = np.append(coords, s.frac_coords)

		grid_encut = (np.pi * self.dim / self.structure.lattice.abc)**2 / 0.262

		self._c_projector_setup(self.num_elems, self.num_sites, max(grid_encut),
								nums, coords, self.dim, pps)

	def check_c_projectors(self):
		"""
		Check to see if the projector functions have been read in and set up.
		If not, do so.
		"""
		if not self.projector_owner:
			start = time.monotonic()
			self._make_c_projectors()
			end = time.monotonic()
			print('--------------\nran setup_projections in %f seconds\n---------------' % (end-start))

	def get_state_realspace(self, b, k, s, dim=None, remove_phase = False):
		"""
		Returns the real and imaginary parts of a given band.
		Args:
			b (int): band number
			k (int): kpoint number
			s (int): spin number
			dim (numpy array of 3 ints): dimensions of the FFT grid
		Returns:
			A 3D array (indexed by x,y,z where x,y,z are fractional coordinates)
				with complex double values for the realspace wavefunction
		"""

		self.check_c_projectors()
		if dim is not None:
			self.update_dim(np.array(dim))
		return self._get_realspace_state(b, k, s, remove_phase)

	def get_state_realspace_density(self, b, k, s, dim=None):
		"""
		Returns the real and imaginary parts of a given band.
		Args:
			b (int): band number
			k (int): kpoint number
			s (int): spin number
			dim (numpy array of 3 ints): dimensions of the FFT grid
		Returns:
			A 3D array (indexed by x,y,z where x,y,z are fractional coordinates)
				with complex double values for the realspace wavefunction
		"""

		self.check_c_projectors()
		if dim is not None:
			self.update_dim(np.array(dim)//2)
		return self._get_realspace_state_density(b, k, s)

	def get_realspace_density(self, dim = None, bands = None):
		"""
		Returns the all electron charge density.
		Args:
			dim (numpy array of 3 ints, None): dimensions of the FFT grid
		Returns:
			A 3D array (indexed by x,y,z where x,y,z are fractional coordinates)
				with real double values for the all electron charge density
		"""
		self.check_c_projectors()
		if dim is not None:
			self.update_dim(np.array(dim)//2)
		return self._get_realspace_density()

	def _convert_to_vasp_volumetric(self, filename, dim):
		"""
		Utility function to convert pawpyseed volumetric
		output to VASP volumetric output.
		"""

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

	def write_state_realspace(self, b, k, s, fileprefix = "", dim=None,
							  scale = 1, remove_phase=False):
		"""
		Writes the real and imaginary parts of a given band to two files,
		prefixed by fileprefix

		Args:
			b (int): band number (0-indexed!)
			k (int): kpoint number (0-indexed!)
			s (int): spin number (0-indexed!)
			fileprefix (string, ""): first part of the file name
			dim (numpy array of 3 ints, None): dimensions of the FFT grid
			scale (scalar, 1): number to multiply the realspace wavefunction by.
				For example, VASP multiplies charge density by the volume
				of the structure.
			remove_phase (False): If True, removes the e^(ikr) phase
				from the wavefunction (this does not necessarily mean
				the wavefunction is real). This is useful if you want
				to visualize the wavefunction because the e^(ikr) phase
				makes the wavefunction non-periodic
		Returns:
			A 3D array (indexed by x,y,z where x,y,z are fractional coordinates)
				with complex double values for the realspace wavefunction
			The wavefunction is written in two files with z the slow index.
		"""
		self.check_c_projectors()
		if dim is not None:
			self.update_dim(np.array(dim))
		filename_base = "%sB%dK%dS%d" % (fileprefix, b, k, s)
		filename1 = "%s_REAL" % filename_base
		filename2 = "%s_IMAG" % filename_base
		res = self._write_realspace_state(filename1, filename2, scale,
										  b, k, s, remove_phase)
		self._convert_to_vasp_volumetric(filename1, self.dim)
		self._convert_to_vasp_volumetric(filename2, self.dim)
		return res

	def write_density_realspace(self, filename = "PYAECCAR", dim=None,
								scale = 1, bands=None):
		"""
		Writes the real and imaginary parts of a given band to two files,
		prefixed by fileprefix

		Args:
			b (int): band number (0-indexed!)
			k (int): kpoint number (0-indexed!)
			s (int): spin number (0-indexed!)
			fileprefix (string, ""): first part of the file name
			dim (numpy array of 3 ints, None): dimensions of the FFT grid
			scale (scalar, 1): number to multiply the realspace wavefunction by.
				For example, VASP multiplies charge density by the volume
				of the structure.
			bands (int or [int], None): Only calculate the density for a specific
				band or set of bands
		Returns:
			A 3D array (indexed by x,y,z where x,y,z are fractional coordinates)
				with complex double values for the realspace wavefunction
			The charge density is written with z the slow index.
		"""

		self.check_c_projectors()
		if dim is not None:
			self.update_dim(np.array(dim)//2)
		res = self._write_realspace_density(filename, scale, bands)
		self._convert_to_vasp_volumetric(filename, self.dim*2)
		return res

	def get_nosym_kpoints(self, init_kpts = None, symprec=None,
		gen_trsym = True, fil_trsym = True):
		"""
		Helper function to get a non-symmetry-reduced k-point
		mesh based on the symmetry-reduced mesh of self.
		"""

		if symprec == None:
			symprec = self.symprec
		return pawpy_symm.get_nosym_kpoints(kpts, self.structure, init_kpts,
										symprec, gen_trsym, fil_trsym)

	def get_kpt_mapping(self, allkpts, symprec=None, gen_trsym = True):
		"""
		Helper function to find the mappings from self.kpts to
		allkpts using the symmetry operations of self.structure
		"""
		
		if symprec == None:
			symprec = self.symprec
		return pawpy_symm.get_kpt_mapping(allkpts, self.kpts, self.structure,
										symprec, gen_trsym)
