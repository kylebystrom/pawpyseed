# coding: utf-8

## @package pawpyseed.core.projector
# Defines the Projector class for evaluating
# AE and PS projection operators.

from pawpyseed.core.wavefunction import *
from pawpyseed.core import pawpyc
from pawpyseed.core.pawpyc import Timer
import warnings

class Projector(pawpyc.CProjector):
	"""
	Projector is a class to project KS states
	from wf onto the KS states of basis
	(both wavefunction objects).

	Attributes:
		wf (Wavefunction): Wavefunction object
		basis (Wavefunction): Wavefunction object onto which the
			wavefunctions of wf are to be projected
		method (str): The method used to perform the projections.
			'aug_real' (default): Filter high-frequency components out
				of the partial waves and then project them onto
				the pseudo wavefunctions in real space inside
				the augmentation spheres
			'aug_recip': Filter high-frequency components out of
				the partial waves, then sum all the partial wave components
				in real space on the FFT grid. Fourier transform
				the result and project it onto the pseudowavefunctions
				in reciprocal space.
			'realspace': Project the full wavefunctions onto realspace
				grids, then integrate over real space.
			'pseudo': Perform projections
				using only plane-wave coefficient components of the
				wavefunctions (pseudo wavefunctions).
				Sacrifices orthogonalization and
				normalization for speed. Not recommended except for
				very rough, qualitative informtation.
	"""

	METHODS = ["pseudo", "realspace", "aug_recip", "aug_real"]

	def __init__(self, wf, basis,
		unsym_basis = False, unsym_wf = False, method = "aug_real"):
		"""
		Arguments:
			wf (Wavefunction): The wavefunction objects whose
				bands are to be projected onto basis
			basis (Wavefunction): The wavefunction whose bands
				serve as the basis set for projection
			unsym_basis (bool, False): If True, makes a copy
				of basis in which the k-point mesh is not symmetrically
				reduced, and then frees the original basis
			unsym_wf (bool, False): If True, makes a copy of
				wf in which the k-point mesh is not symmetrically
				reduced, and the frees the original wf
			method (str, "aug_recip"): Options: "pseudo", "realspace", "aug_recip", "aug_real";
				The method to use for the projections. See method
				options in the Attributes section.

		Returns:
			Projector object
		"""
		self.method = method
		if self.method == "pseudo":
			self._single_band_projection = self._single_band_projection_pseudo
		elif self.method == "realspace":
			self._single_band_projection = self._single_band_projection_realspace
		elif self.method == "aug_recip":
			self._single_band_projection = self._single_band_projection_aug_recip
		elif self.method == "aug_real":
			self._single_band_projection = self._single_band_projection_aug_real
		else:
			raise PAWpyError("method not recognized for Projector")

		if wf.ncl or basis.ncl:
			raise PAWpyError("Projection not supported for noncollinear case!")

		if unsym_basis and unsym_wf:
			basis = basis.desymmetrized_copy()
			wf = wf.desymmetrized_copy(basis.kpts, basis.kws)
		elif unsym_wf and not unsym_basis:
			if basis.kpts.shape[0] < wf.kpts.shape[0]:
				raise PAWpyError("Basis doesn't have enough kpoints, needs to be desymmetrized!")
			allkpts = basis.kpts
			weights = basis.kws
			wf = wf.desymmetrized_copy(allkpts, weights)
		elif unsym_basis and not unsym_wf:
			if wf.kpts.shape[0] < basis.kpts.shape[0]:
				raise PAWpyError("Defect doesn't have enough kpoints, needs to be desymmetrized!")
			allkpts = wf.kpts
			weights = wf.kws
			basis = basis.desymmetrized_copy(allkpts, weights)

		if np.linalg.norm(basis.kpts - wf.kpts) > 1e-10:
			raise PAWpyError("k-point grids for projection are not matched.")
		if np.linalg.norm(basis.kws - wf.kws) > 1e-10:
			raise PAWpyError("k-point weights for projection are not matched.")
		if wf.structure.lattice != basis.structure.lattice:
			raise PAWpyError("Need the lattice to be the same for projections, and they are not")

		if self.method != "pseudo":
			basis.check_c_projectors()
			wf.check_c_projectors()

		super(Projector, self).__init__(wf, basis)

		if "aug" in self.method:
			self.setup_overlap()

	def make_site_lists(self):
		"""
		Organizes sites into sets for use in the projection scheme. M_R and M_S contain site indices
		of sites which are identical in structures R (basis) and S (self). N_R and N_S contain all other
		site indices, and N_RS contains pairs of indices in R and S with overlapping augmentation
		spheres in the PAW formalism. R si for self.basis, S is for self.wf

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

		ref_sites = self.basis.structure.sites
		sites = self.wf.structure.sites
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
				rmax1 = self.basis.cr.pps[el(ref_sites[i])].rmax
				rmax2 = self.wf.cr.pps[el(sites[i])].rmax
				if ref_sites[i].distance(sites[j]) < rmax1 + rmax2:
					N_RS.append((i,j))
		return M_R, M_S, N_R, N_S, N_RS

	def setup_overlap(self):
		"""
		Evaluates projectors <p_i|psi>, as well
		as <(phi-phit)|psi> and <(phi_i-phit_i)|(phi_j-phit_j)>,
		when needed
		"""
		M_R, M_S, N_R, N_S, N_RS = self.make_site_lists()
		num_N_RS = len(N_RS)
		if num_N_RS > 0:
			N_RS_R, N_RS_S = zip(*N_RS)
		else:
			N_RS_R, N_RS_S = [], []
		self.site_cat = [M_R, M_S, N_R, N_S, N_RS_R, N_RS_S]
		start = time.monotonic()
		if self.method == "aug_recip":
			self._setup_overlap(self.site_cat, True)
		elif self.method == "aug_real":
			self._setup_overlap(self.site_cat, False)
		else:
			raise PAWpyError("method must be aug type for setup_overlap call")
		end = time.monotonic()
		Timer.overlap_time(end-start)
		print('-------------\nran overlap_setup in %f seconds\n---------------' % (end-start))

	def _single_band_projection_pseudo(self, band_num):
		"""
		Very rough approximation for the projection of the band_num band of self
		onto all the bands of basis. Simply sums over the plane
		waves of the pseudo wavefunctions to calculate overlaps,
		and the pseudo wavefunctions are not orthonormal. Only use if
		you just want a quick, rough look at band character.
		"""
		return self.wf.pseudoprojection(band_num, self.basis)

	def _single_band_projection_realspace(self, band_num, dim = None):
		"""
		All electron projection of the band_num band of self
		onto all the bands of basis. The wavefunctions are evaluated
		and integrated on a real space FFT grid, with the default
		dimension being the fine FFT grid from VASP.
		"""
		if dim == None:
			dim = self.wf.dim * 2
		return self._realspace_projection(band_num, dim)

	def _single_band_projection_aug_real(self, band_num, flip_spin=False):
		"""
		All electron projection of the band_num band of self
		onto all the bands of basis. High frequency components
		are filtered out from the partial waves, which are then
		projected onto the pseudo wavefunctions in real space.
		"""
		res = self.wf.pseudoprojection(band_num, self.basis, flip_spin)
		start = time.monotonic()
		self._add_augmentation_terms(res, band_num, flip_spin)
		end = time.monotonic()
		Timer.augmentation_time(end-start)
		#print('---------\nran compensation_terms in %f seconds\n-----------' % (end-start))
		return res

	def _single_band_projection_aug_recip(self, band_num, flip_spin=False):
		"""
		All electron projection of the band_num band of self
		onto all the bands of basis. High frequency components
		are filtered out from the partial waves, which are then
		projected into real space, summed, projected back into reciprocal
		space, and then projected onto the pseudo wavefunction
		in reciprocal space.
		"""
		res = self.wf.pseudoprojection(band_num, self.basis, flip_spin)
		self._projection_recip(res, band_num, flip_spin)
		return res

	def single_band_projection(self, band_num, **kwargs):
		"""
		Projection of the band_num band of self
		onto all the bands of basis. Returned as a numpy array,
		with the overlap operator matrix elements ordered as follows:
		loop over band
			loop over spin
				loop over kpoint

		Arguments:
			band_num (int): band which is projected onto basis
			kwargs:
				For method=='realspace':
					dim (tuple or numpy array): Set a custom grid dimension
				For method=='aug_*':
					flip_spin (bool): If true, swap the spins of basis
						for the projection. I.e. instead of
						<basis;b,k,s|wf;b0,k,s>, you get
						<basis;b,k,1-s|wf;b0,k,s> if and only if nspin==2

		Returns:
			(np.array): overlap operator expectation values
				as described above
		
		Example:
			# Get overlap of band b0, k-point k, spin s (0 or 1 index)
			# of wf with band b, k-point k, spin s of basis
			# <basis;b,k,s|wf;b0,k,s>
			>>> pr = Projector(wf, basis)
			>>> res = pr.single_band_projection(b0)
			>>> print(res[b*pr.nwk*pr.nspin + s*pr.nwk + k])
		"""
		if band_num >= self.wf.nband or band_num < 0:
			raise ValueError("Band index out of range (0-indexed)")
		return self._single_band_projection(band_num, **kwargs)

	@staticmethod
	def setup_bases(basis_dirs, desymmetrize = True,
		atomate_compatible = True):
		"""
		This convenience function performs the setup
		of all the bases in the basis_dirs list.

		Arguments:
			basis_dir (list of str): paths to the VASP outputs
				to be used as the basis structures
			desymmetrize (bool, False): If True, constructs
				Wavefunction objects in which the k-point mesh
				is not symmetrically reduced
			atomate_compatible (bool, True): If True, checks for the gzipped
				files created the atomate workflow tools and reads the most
				recent run based on title

		Returns:
			list of Wavefunction objects, each basis in the same
				order as the basis_dirs list
		"""

		bases = []
		crs = []
		for bdir in basis_dirs:
			if atomate_compatible:
				basis = Wavefunction.from_atomate_directory(bdir, False)
			else:
				basis = Wavefunction.from_directory(bdir, False)

			if desymmetrize:
				basis = basis.desymmetrized_copy()

			basis.check_c_projectors()
			bases.append(basis)

		return bases

	@staticmethod
	def setup_multiple_projections(basis_dir, wf_dirs, method = "aug_real",
									ignore_errors = False,
									desymmetrize = False,
									atomate_compatible = True):
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
			desymmetrize (bool, False): If True, constructs
				Wavefunction objects in which the k-point mesh
				is not symmetrically reduced
			atomate_compatible (bool, True): If True, checks for the gzipped
				files created by the atomate workflow tools and reads the most
				recent run based on title

		Returns:
			list -- wf_dir, basis, wf
			Each iteration of the generator function returns a directory name from
			wf_dirs (wf_dir), the basis Wavefunction object (basis), and the Wavefunction
			object associated with wf_dir (wf), fully setup to project bands of wf
			onto bands of basis.
		"""

		if atomate_compatible:
			basis = Wavefunction.from_atomate_directory(basis_dir, False)
		else:
			basis = Wavefunction.from_directory(basis_dir, False)
		
		if desymmetrize:
			basis = basis.desymmetrized_copy()

		errcount = 0
		pr = None
		for wf_dir in wf_dirs:

			try:
				if atomate_compatible:
					wf = Wavefunction.from_atomate_directory(wf_dir, False)
				else:
					wf = Wavefunction.from_directory(wf_dir, False)

				if desymmetrize:
					pr = Projector(wf, basis, method = method, unsym_wf = True)
				else:
					pr = Projector(wf, basis, method = method)

				yield [wf_dir, pr]
			except Exception as e:
				if ignore_errors:
					errcount += 1
				else:
					raise PAWpyError('Unable to setup wavefunction in directory %s' % wf_dir\
										+'\nGot the following error:\n'+str(e))
		if not pr:
			raise PAWpyError("Could not generate any projector setups")
		print("Number of errors:", errcount)

	def proportion_conduction(self, band_num, spinpol = False):
		"""
		Calculates the proportion of band band_num in self
		that projects onto the valence states and conduction
		states of self.basis (should be the bulk structure).
		Designed for analysis of point defect
		wavefunctions.

		Arguments:
			band_num (int): number of defect bands in self
			spinpol (bool, False): whether to return separate
				values of the projection for spin up and spin down

		Returns:
			v, c (int, int): The valence (v) and conduction (c)
				proportion of band band_num
		"""
		basis = self.basis
		nband = basis.nband
		nwk = basis.nwk
		nspin = basis.nspin
		occs = self.basis._get_occs()

		res = self.single_band_projection(band_num)

		if spinpol:
			c, v = np.zeros(nspin), np.zeros(nspin)
			for b in range(nband):
				for s in range(nspin):
					ind = b*nspin+s
					prop = np.absolute(res[ind*nwk:(ind+1)*nwk])**2
					c[s] += np.dot(prop, (1-occs[ind*nwk:(ind+1)*nwk]) * self.wf.kws)
					v[s] += np.dot(prop, occs[ind*nwk:(ind+1)*nwk] * self.wf.kws)
		else:
			c, v = 0, 0
			for i in range(nband*nwk*nspin):
				if occs[i] > 0.5:
					v += np.absolute(res[i]) ** 2 * self.wf.kws[i%nwk] / nspin
				else:
					c += np.absolute(res[i]) ** 2 * self.wf.kws[i%nwk] / nspin
		if self.method == "pseudo":
			t = v+c
			v /= t
			c /= t
		if spinpol:
			v = v.tolist()
			c = c.tolist()
		return v, c

	def defect_band_analysis(self, num_below_ef=20,
		num_above_ef=20, spinpol = False, return_energies = False,
		vbmband = None, band_list = None, analyze_all = False):
		"""
		Identifies a set of 'interesting' bands in a defect structure
		to analyze by choosing any band that is more than bound conduction
		and more than bound valence in the pseudoprojection scheme,
		and then fully analyzing these bands using single_band_projection.
		NOTE: ALL BANDS ARE ZERO-INDEXED!

		Args:
			num_below_ef (int, 20): number of bands to analyze below the fermi level
			num_above_ef (int, 20): number of bands to analyze above the fermi level
			spinpol (bool, False): whether to return spin-polarized results (only allowed
				for spin-polarized DFT output)
			return_energies (bool, False): whether to return the energy levels
				of the bands analyzed in the form
				{band : [kpoint label: {spin label: (energy, occupation)}},
				where the kpoint label and spin label are integers
			vbmband (int, None): Optionally allows a VBM band number to be specified.
				If None, the VBM of wf is determined and used.
			band_list (list of int, None): If not None, overrides num_below_ef,
				num_above_ef, and vbmband. Specifies the set of bands to analyze.
			analyze_all (bool, False): If True, overrides num_below_ef,
				num_above_ef, vbmband, and band_list. Whether to perform
				analysis on all bands in wf
		"""
		if num_below_ef < 0 or num_above_ef < 0:
			raise ValueError("num_above_ef and num_below_ef must both be nonnegative.")

		basis = self.basis
		nband = basis.nband
		nwk = basis.nwk
		nspin = basis.nspin
		occs = self.wf._get_occs()
		
		if analyze_all:
			totest = [i for i in range(nband)]
		elif band_list:
			totest = band_list[:]
		else:
			vbm = 0
			for i in range(self.wf.nband):
				if occs[i*self.wf.nwk*self.wf.nspin] > 0.5:
					vbm = i
			if vbmband != None:
				vbm = vbmband
			min_band, max_band = max(vbm - num_below_ef, 0), min(vbm + num_above_ef, self.wf.nband - 1)
			totest = [i for i in range(min_band,max_band+1)]

		results = {}
		for b in totest:
			results[b] = self.proportion_conduction(b, spinpol = spinpol)

		if return_energies:
			return results, self.wf._get_energy_list(totest)
		else:
			return results
