from pawpyseed.core.wavefunction import *

class NCLWavefunction(pawpyc.CNCLWavefunction, Wavefunction):

	def __init__(self, struct, pwf, cr, dim, symprec=1e-4, setup_projectors=False):
		"""
		Arguments:
			struct (pymatgen.core.Structure): structure that the wavefunction describes
			pwf (pawpyc.PWFPointer): holder class for pswf_t and k-points/k-point weights
			cr (CoreRegion): Contains the pseudopotentials, with projectors and
				partials waves, for the structure
			dim (pymatgen.io.vasp.outputs.Outcar OR np.ndarry OR list of length 3):
				Outcar object for reading ngf or the dimensions NG* of the FFT grid
			setup_projectors (bool, False): Whether to set up the core region
				components of the wavefunctions. Pawpyseed will set up the projectors
				automatically when they are first needed, so this generally
				can be left as False.
		Returns:
			Wavefunction object
		"""
		self.band_props = pwf.band_props.copy(order = 'C')
		super(Wavefunction, self).__init__(pwf)
		if not self.ncl:
			raise PAWpyError("Pseudowavefunction is collinear! Call Wavefunction(...) instead")
		self.structure = struct
		self.cr = cr
		self.dim = np.array(dim).astype(np.int32)
		if setup_projectors:
			self.check_c_projectors()

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
		vr = Vasprun(vr)
		dim = np.array([vr.parameters["NGX"], vr.parameters["NGY"], vr.parameters["NGZ"]])
		symprec = vr.parameters["SYMPREC"]
		pwf = pawpyc.PWFPointer(wavecar, vr)
		return NCLWavefunction(Poscar.from_file(struct).structure,
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
		return NCLWavefunction.from_files(*args)

	def desymmetrized_copy(self, allkpts = None, weights = None):
		raise NotImplementedError()

	def write_state_realspace(self, b, k, s, fileprefix = "", dim=None, scale = 1,
								remove_phase=False):
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
		filename1 = "%s_UP_REAL" % filename_base
		filename2 = "%s_DOWN_REAL" % filename_base
		filename3 = "%s_UP_IMAG" % filename_base
		filename4 = "%s_DOWN_IMAG" % filename_base
		res0, res1 = self._write_realspace_state(filename1, filename2, filename3, filename4,
											scale, b, k, s)
		self._convert_to_vasp_volumetric(filename1, dim)
		self._convert_to_vasp_volumetric(filename2, dim)
		self._convert_to_vasp_volumetric(filename3, dim)
		self._convert_to_vasp_volumetric(filename4, dim)
		return res0, res1

	def write_density_realspace(self, filename = "PYAECCAR", dim=None, scale = 1):
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
		Returns:
			A 3D array (indexed by x,y,z where x,y,z are fractional coordinates)
				with complex double values for the realspace wavefunction
			The charge density is written with z the slow index.
		"""

		self.check_c_projectors()
		if dim is not None:
			self.update_dim(np.array(dim))
		res = self._write_realspace_density(filename, scale)
		self._convert_to_vasp_volumetric(filename, dim)
		return res
