from pawpyseed.core.wavefunction import *

class Projector(Wavefunction):

	def __init__(wf, basis, setup_basis=True):
		"""
		Projector is a class to projector KS states
		from wf onto the KS states of basis
		(both wavefunction objects). Projector extends
		Wavefunction, so a Projector is essentially
		a Wavefunction object that is set up to be
		projected onto basis.
		"""
		self.structure = wf.structure
		self.pwf = wf.pwf
		self.cr = wf.cr
		self.dim = wf.dim
		self.projector_owner = True
		self.projector_list = None
		self.nband = wf.nband
		self.nwk = wf.nwk
		self.nspin = wf.nspin
		self.nums = wf.nums
		self.coords = wf.coords
		self.basis = basis
		self.wf = wf

		self.setup_projection(basis, setup_basis)

	@staticmethod
	def from_files(basis, struct="CONTCAR", pwf="WAVECAR", cr="POTCAR",
		vr="vasprun.xml", outcar="OUTCAR", setup_basis=True):
		"""
		Construct a Projector object from file paths.
		Arguments:
			basis (Wavefunction): the basis Wavefunction
			struct (str): VASP POSCAR or CONTCAR file path
			pwf (str): VASP WAVECAR file path
			vr (str): VASP vasprun file path
			outcar (str): VASP OUTCAR file path
			setup_basis (bool): whether to set up the basis
		Returns:
			Projector object
		"""
		return Projector(Poscar.from_file(struct).structure,
			PseudoWavefunction(pwf, vr),
			CoreRegion(Potcar.from_file(cr)),
			Outcar(outcar), setup_basis)

	@staticmethod
	def from_directory(basis, path, setup_basis=True):
		"""
		Assumes VASP output has the default filenames and is located
		in the directory specificed by path.
		Arguments:
			basis (Wavefunction): the basis Wavefunction
			path (str): path to the VASP calculation directory
			setup_basis (bool): whether to set up the basis
		Returns:
			Projector object
		"""
		filepaths = []
		for d in ["CONTCAR", "WAVECAR", "POTCAR", "vasprun.xml", "OUTCAR"]:
			filepaths.append(str(os.path.join(path, d)))
		args = filepaths + [setup_basis]
		return Projector.from_files(*args)

	def setup_projection(self, basis, setup_basis=True):
		"""
		Evaluates projectors <p_i|psi>, as well
		as <(phi-phit)|psi> and <(phi_i-phit_i)|(phi_j-phit_j)>,
		when needed

		Arguments:
			basis (Wavefunction): wavefunction onto which bands of self
			will be projected.
		"""

		#if not basis.projection_data:
		#	basis.projection_data = self.make_c_projectors(basis)
		#projector_list, selfnums, selfcoords, basisnums, basiscoords = basis.projection_data
		if setup_basis:
			self.projector_list, self.nums, self.coords,\
				basis.nums, basis.coords = self.make_c_projectors(basis)
		projector_list = self.projector_list
		basisnums = basis.nums
		basiscoords = basis.coords
		selfnums = self.nums
		selfcoords = self.coords

		print(hex(projector_list), hex(self.pwf.wf_ptr))
		sys.stdout.flush()
		print ("TYPETHING", basis.pwf.wf_ptr, type(basis.pwf.wf_ptr))
		
		if setup_basis:
			cfunc_call(PAWC.setup_projections, None,
						basis.pwf.wf_ptr, projector_list,
						self.num_proj_els, len(basis.structure), self.dim,
						basisnums, basiscoords)
		start = time.monotonic()
		cfunc_call(PAWC.setup_projections_copy_rayleigh, None,
					self.pwf.wf_ptr, basis.pwf.wf_ptr,
					projector_list, self.num_proj_els, len(self.structure),
					self.dim, selfnums, selfcoords)
		end = time.monotonic()
		print('--------------\nran setup_projections in %f seconds\n---------------' % (end-start))
		Timer.setup_time(end-start)
		M_R, M_S, N_R, N_S, N_RS = self.make_site_lists(basis)
		num_N_RS = len(N_RS)
		if num_N_RS > 0:
			N_RS_R, N_RS_S = zip(*N_RS)
		else:
			N_RS_R, N_RS_S = [], []
		self.site_cat = [M_R, M_S, N_R, N_S, N_RS_R, N_RS_S]
		start = time.monotonic()
		cfunc_call(PAWC.overlap_setup_real, None, basis.pwf.wf_ptr, self.pwf.wf_ptr,
					projector_list, basisnums, selfnums, basiscoords, selfcoords,
					N_R, N_S, N_RS_R, N_RS_S, len(N_R), len(N_S), len(N_RS_R))
		#cfunc_call(PAWC.overlap_setup_real, None, basis.pwf.wf_ptr, self.pwf.wf_ptr,
		#			projector_list, basisnums, selfnums, basiscoords, selfcoords,
		#			M_R, M_S, M_R, M_S, len(M_R), len(M_R), len(M_R))
		end = time.monotonic()
		Timer.overlap_time(end-start)
		print('-------------\nran overlap_setup in %f seconds\n---------------' % (end-start))

	def single_band_projection(self, band_num, basis):
		"""
		All electron projection of the band_num band of self
		onto all the bands of basis. Returned as a numpy array,
		with the overlap operator matrix elements ordered as follows:
		loop over band
			loop over spin
				loop over kpoint

		Arguments:
			band_num (int): band which is projected onto basis
			basis (Wavefunction): basis Wavefunction object

		Returns:
			res (np.array): overlap operator expectation values
				as described above
		"""

		nband = PAWC.get_nband(c_void_p(basis.pwf.wf_ptr))
		nwk = PAWC.get_nwk(c_void_p(basis.pwf.wf_ptr))
		nspin = PAWC.get_nspin(c_void_p(basis.pwf.wf_ptr))
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
		
		start = time.monotonic()
		ct = cfunc_call(PAWC.compensation_terms, 2*nband*nwk*nspin,
						band_num, self.pwf.wf_ptr, basis.pwf.wf_ptr,
						projector_list, len(self.cr.pps),
						len(M_R), len(N_R), len(N_S), len(N_RS_R),
						M_R, M_S, N_R, N_S, N_RS_R, N_RS_S,
						selfnums, selfcoords, basisnums, basiscoords,
						self.dim)
		end = time.monotonic()
		Timer.augmentation_time(end-start)
		print('---------\nran compensation_terms in %f seconds\n-----------' % (end-start))
		"""
		ct = cfunc_call(PAWC.compensation_terms, 2*nband*nwk*nspin,
						band_num, self.pwf.wf_ptr, basis.pwf.wf_ptr,
						projector_list, len(self.cr.pps),
						0, len(M_R), len(M_S), len(M_S),
						np.array([]), np.array([]), M_R, M_S, M_R, M_S,
						selfnums, selfcoords, basisnums, basiscoords,
						self.dim)
		"""
		res += ct
		return res[::2] + 1j * res[1::2]

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
			print ("THIS IS THE NUM %d" % num)
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
		#return PAWC.get_projector_list(num_els, numpy_to_cint(clabels),
		#	numpy_to_cint(ls), numpy_to_cdouble(pgrids), numpy_to_cdouble(wgrids),
		#	numpy_to_cdouble(projectors), numpy_to_cdouble(aewaves), numpy_to_cdouble(pswaves),
		#	numpy_to_cdouble(rmaxs), max(grid_encut))
		return cfunc_call(PAWC.get_projector_list, None,
							num_els, clabels, ls, pgrids, wgrids,
							projectors, aewaves, pswaves,
							rmaxs, max(grid_encut))

	@staticmethod
	def setup_multiple_projections(basis_dir, wf_dirs, ignore_errors = False):
		"""
		A convenient generator function for processing the Kohn-Sham wavefunctions
		of multiple structures with respect to one structure used as the basis.
		All C memory is freed after each yield for the wavefunctions to be analyzed,
		and C memory associated with the basis wavefunction is freed when
		the generator is called after all wavefunctions have been yielded.

		Args:
			basis_dir (str): path to the VASP output to be used as the basis structure
			wf_dirs (list of str): paths to the VASP outputs to be analyzed
			ignore_errors (bool, False): whether to ignore errors in setting up
				Wavefunction objects by skipping over the directories for which
				setup fails.

		Returns:
			list -- wf_dir, basis, wf
			Each iteration of the generator function returns a directory name from
			wf_dirs (wf_dir), the basis Wavefunction object (basis), and the Wavefunction
			object associated with wf_dir (wf), fully setup to project bands of wf
			onto bands of basis.
		"""

		basis = Wavefunction.from_directory(basis_dir)
		crs = [basis.cr] + [CoreRegion(Potcar.from_file(os.path.join(wf_dir, 'POTCAR'))) \
			for wf_dir in wf_dirs]

		pps = {}
		labels = {}
		label = 0
		for cr in crs:
			for e in cr.pps:
				if not e in labels:
					pps[label] = cr.pps[e]
					labels[e] = label
					label = label + 1

		#print (pps)
		projector_list = basis.get_c_projectors_from_pps(pps)
		basisnums = np.array([labels[el(s)] for s in basis.structure], dtype=np.int32)
		basiscoords = np.array([], np.float64)
		for s in basis.structure:
			basiscoords = np.append(basiscoords, s.frac_coords)
		projector_list = basis.projector_list
		basis.nums = basisnums
		basis.coords = basiscoords
		basis.num_proj_els = len(pps)

		sys.stdout.flush()	
		#PAWC.setup_projections(c_void_p(basis.pwf.wf_ptr),
		#	c_void_p(projector_list), label,
		#	len(basis.structure), numpy_to_cint(basis.dim), numpy_to_cint(basisnums),
		#	numpy_to_cdouble(basiscoords))
		cfunc_call(PAWC.setup_projections, None,
					basis.pwf.wf_ptr, projector_list, label,
					len(basis.structure), basis.dim, basisnums, basiscoords)

		for wf_dir, cr in zip(wf_dirs, crs[1:]):

			try:
				files = {}
				for f in ['CONTCAR', 'OUTCAR', 'vasprun.xml', 'WAVECAR']:
					files[f] = os.path.join(wf_dir, f)
				struct = Poscar.from_file(files['CONTCAR']).structure
				pwf = PseudoWavefunction(files['WAVECAR'], files['vasprun.xml'])
				outcar = Outcar(files['OUTCAR'])

				wf = Wavefunction(struct, pwf, cr, outcar, False)

				selfnums = np.array([labels[el(s)] for s in wf.structure], dtype=np.int32)
				selfcoords = np.array([], np.float64)

				for s in wf.structure:
					selfcoords = np.append(selfcoords, s.frac_coords)
				wf.nums = selfnums
				wf.coords = selfcoords
				wf.num_proj_els = len(pps)

				pr = Projector(wf, basis, setup_basis=False)

				yield [wf_dir, basis, pr]
				wf.free_all()
				pr.free_all()
			except Exception as e:
				if ignore_errors:
					errcount += 1
				else:
					raise PAWpyError('Unable to setup wavefunction in directory %s' % wf_dir\
										+'\nGot the following error:\n'+str(e))

		basis.free_all()
		print("Number of errors:", errcount)
			

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
		#print (pps, labels)
		if basis != None:
			for e in basis.cr.pps:
				if not e in labels:
					pps[label] = basis.cr.pps[e]
					labels[e] = label
					label += 1
		
		projector_list = self.get_c_projectors_from_pps(pps)

		selfnums = np.array([labels[el(s)] for s in self.structure], dtype=np.int32)
		selfcoords = np.array([], np.float64)
		if basis != None:
			basisnums = np.array([labels[el(s)] for s in basis.structure], dtype=np.int32)
			basiscoords = np.array([], np.float64)

		self.num_proj_els = len(pps)
		if basis != None:
			basis.num_proj_els = len(pps)
		for s in self.structure:
			selfcoords = np.append(selfcoords, s.frac_coords)
		if basis != None:
			for s in basis.structure:
				basiscoords = np.append(basiscoords, s.frac_coords)
			return projector_list, selfnums, selfcoords, basisnums, basiscoords
		return projector_list, selfnums, selfcoords

	def proportion_conduction(self, band_num, bulk, pseudo = False, spinpol = False):
		"""
		Calculates the proportion of band band_num in self
		that projects onto the valence states and conduction
		states of bulk. Designed for analysis of point defect
		wavefunctions.

		Arguments:
			band_num (int): number of defect band in self
			bulk (Wavefunction): wavefunction of bulk crystal
				with the same lattice and basis set as self

		Returns:
			v, c (int, int): The valence (v) and conduction (c)
				proportion of band band_num
		"""

		nband = PAWC.get_nband(c_void_p(bulk.pwf.wf_ptr))
		nwk = PAWC.get_nwk(c_void_p(bulk.pwf.wf_ptr))
		nspin = PAWC.get_nspin(c_void_p(bulk.pwf.wf_ptr))
		occs = cdouble_to_numpy(PAWC.get_occs(c_void_p(bulk.pwf.wf_ptr)), nband*nwk*nspin)

		if pseudo:
			res = self.pwf.pseudoprojection(band_num, bulk.pwf)
		else:
			res = self.single_band_projection(band_num, bulk)

		if spinpol:
			c, v = np.zeros(nspin), np.zeros(nspin)
			for b in range(nband):
				for s in range(nspin):
					for k in range(nwk):
						i = b*nspin*nwk + s*nwk + k
						if occs[i] > 0.5:
							v += np.absolute(res[i]) ** 2 * self.pwf.kws[i%nwk] / nspin
						else:
							c += np.absolute(res[i]) ** 2 * self.pwf.kws[i%nwk] / nspin
		else:
			c, v = 0, 0
			for i in range(nband*nwk*nspin):
				if occs[i] > 0.5:
					v += np.absolute(res[i]) ** 2 * self.pwf.kws[i%nwk] / nspin
				else:
					c += np.absolute(res[i]) ** 2 * self.pwf.kws[i%nwk] / nspin
		if pseudo:
			t = v+c
			v /= t
			c /= t
		if spinpol:
			v = v.tolist()
			c = c.tolist()
		return v, c

	def defect_band_analysis(self, bulk, num_below_ef=20, num_above_ef=20, spinpol = False):
		"""
		Identifies a set of 'interesting' bands in a defect structure
		to analyze by choosing any band that is more than bound conduction
		and more than bound valence in the pseudoprojection scheme,
		and then fully analyzing these bands using single_band_projection

		Args:
			bulk (Wavefunction object): bulk structure wavefunction
			num_below_ef (int, 20): number of bands to analyze below the fermi level
			num_above_ef (int, 20): number of bands to analyze above the fermi level
			spinpol (bool, False): whether to return spin-polarized results (only allowed
				for spin-polarized DFT output)
		"""
		nband = PAWC.get_nband(c_void_p(bulk.pwf.wf_ptr))
		nwk = PAWC.get_nwk(c_void_p(bulk.pwf.wf_ptr))
		nspin = PAWC.get_nspin(c_void_p(bulk.pwf.wf_ptr))
		#totest = set()
		occs = cdouble_to_numpy(PAWC.get_occs(c_void_p(bulk.pwf.wf_ptr)), nband*nwk*nspin)
		vbm = 0
		for i in range(nband):
			if occs[i*nwk*nspin] > 0.5:
				vbm = i
		min_band, max_band = vbm - num_below_ef, vbm + num_above_ef
		if min_band < 0 or max_band >= nband:
			raise PAWpyError("The min or max band is too large/small with min_band=%d, max_band=%d, nband=%d" % (min_band, max_band, nband))
		"""
		for b in range(nband):
			v, c = self.proportion_conduction(b, bulk, pseudo = True, spinpol = False)
			if v > bound and c > bound:
				totest.add(b)
				totest.add(b-1)
				totest.add(b+1)
		"""
		totest = [i for i in range(min_band,max_band+1)]
		print("NUM TO TEST", len(totest))

		results = {}
		for b in totest:
			results[b] = self.proportion_conduction(b, bulk, pseudo = False, spinpol = spinpol)

		return results

	def check_c_projectors(self):
		"""
		Check to see if the projector functions have been read in and set up.
		If not, do so.
		"""
		pass

	def free_all(self):
		"""
		Frees all of the C structures associated with the Wavefunction object.
		After being called, this object is not usable.
		"""
		PAWC.free_ppot_list(c_void_p(self.projector_list), len(self.cr.pps))