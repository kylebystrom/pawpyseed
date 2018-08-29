from pawpyseed.core.wavefunction import *

class NCLPseudoWavefunction(PseudoWavefunction):
	"""
	Class for storing noncollinear pseudowavefunction from WAVECAR file.
	Most important attribute
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
		self.kws = weights
		self.kpts = vr.actual_kpoints
		self.wf_ptr = PAWC.read_wavefunctions(filename.encode('utf-8'), kws)

class NCLWavefunction(Wavefunction):

	def __init__(self, struct, pwf, cr, outcar, setup_projectors=True):
		"""
		Arguments:
			struct (pymatgen.core.Structure): structure that the wavefunction describes
			pwf (PseudoWavefunction): Pseudowavefunction componenet
			cr (CoreRegion): Contains the pseudopotentials, with projectors and
				partials waves, for the structure
			outcar (pymatgen.io.vasp.outputs.Outcar): Outcar object for reading ngf
		Returns:
			Wavefunction object
		"""
		
		if type(pwf) != NCLPseudoWavefunction:
			raise PAWpyError("Need NCLPseudoWavefunction to initialize NCLWavefunction")
		self.structure = struct
		self.pwf = pwf
		self.cr = cr
		self.dim = outcar.ngf
		self.dim = np.array(self.dim).astype(np.int32) // 2
		self.projector_owner = False
		self.projector_list = None
		self.nums = None
		self.coords = None
		if setup_projectors:
			self.check_c_projectors()
		self.nband = PAWC.get_nband(c_void_p(pwf.wf_ptr))
		self.nwk = PAWC.get_nwk(c_void_p(pwf.wf_ptr))
		self.nspin = PAWC.get_nspin(c_void_p(pwf.wf_ptr))
		self.num_proj_els = None

	@staticmethod
	def from_files(struct="CONTCAR", pwf="WAVECAR", cr="POTCAR",
		vr="vasprun.xml", outcar="OUTCAR", setup_projectors=True):
		"""
		Construct a Wavefunction object from file paths.
		Arguments:
			struct (str): VASP POSCAR or CONTCAR file path
			pwf (str): VASP WAVECAR file path
			vr (str): VASP vasprun file path
		Returns:
			Wavefunction object
		"""
		return NCLWavefunction(Poscar.from_file(struct).structure,
			NCLPseudoWavefunction(pwf, vr),
			CoreRegion(Potcar.from_file(cr)),
			Outcar(outcar), setup_projectors)

	@staticmethod
	def from_directory(path, setup_projectors=True):
		"""
		Assumes VASP output has the default filenames and is located
		in the directory specificed by path.
		"""
		filepaths = []
		for d in ["CONTCAR", "WAVECAR", "POTCAR", "vasprun.xml", "OUTCAR"]:
			filepaths.append(str(os.path.join(path, d)))
		args = filepaths + [setup_projectors]
		return NCLWavefunction.from_files(*args)

	def write_state_realspace(self, b, k, s, fileprefix = "", dim=None):
		"""
		Writes the real and imaginary parts of a given band to two files,
		prefixed by fileprefix

		Args:
			b (int): band number (0-indexed!)
			k (int): kpoint number (0-indexed!)
			s (int): spin number (0-indexed!)
			dim (numpy array of 3 ints): dimensions of the FFT grid
			fileprefix (string, optional): first part of the file name
		Returns:
			None
			The wavefunction is written with z the slow index.
		"""
		print("PARAMETERS", self.nums, self.coords, dim)
		sys.stdout.flush()
		self.check_c_projectors()
		if type(dim) == type(None):
			dim = self.dim
		filename_base = "%sB%dK%dS%d" % (fileprefix, b, k, s)
		filename1 = "%s_UP_REAL" % filename_base
		filename2 = "%s_DOWN_REAL" % filename_base
		filename3 = "%s_UP_IMAG" % filename_base
		filename4 = "%s_DOWN_IMAG" % filename_base
		print("PARAMETERS", self.nums, self.coords, dim)
		sys.stdout.flush()
		"""
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
		"""
		cfunc_call(PAWC.write_realspace_state_ncl_ri, None, filename1, filename2,
			filename3, filename4,
			b, k+s*self.nwk,
			self.pwf.wf_ptr, self.projector_list,
			dim, self.nums, self.coords)
		self._convert_to_vasp_volumetric(filename1, dim)
		self._convert_to_vasp_volumetric(filename2, dim)
		self._convert_to_vasp_volumetric(filename3, dim)
		self._convert_to_vasp_volumetric(filename4, dim)

