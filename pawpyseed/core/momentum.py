from pawpyseed.core.wavefunction import Wavefunction
from pawpyseed.core import pawpyc

class MomentumMatrix(pawpyc.CMomentumMatrix):

	def __init__(self, wf, encut = None):
		"""
		Class for scattering applications, i.e. calculating
		< psi1 | exp(iGr) | psi2 > for some set of G
		and bands psi1 and psi2. Note that because all of these
		functions are relatively computationally expensive,
		all setup is done in the Cython basis class
		CMomentumMatrix (see pawpyc.pyx), and the numerical
		routines are performed in momentum.c

		Args:
			wf (Wavefunction) : Wavefunction object for system
				of interest
			encut (number, None) : Energy cutoff. If None, defaults
				to 4 * wf.encut. Energy cutoff determines the size
				of self.momentum_grid, which is used for two purposes:
				1) When get_momentum_matrix_elems is called,
				< psi1 | exp(iGr) | psi2 > is calculated for every
				value of G in self.momentum_grid
				2) When get_reciprocal_fullfw is called, the wavefunction
				is projected onto self.momentum_grid.
				NOTE: self.momentum_grid is NOT k-point dependent, unlike
				the G-grids in VASP. That means that it includes all
				G-points with energy <= encut, but the set (k+G) will
				include or exclude a few edge cases. If this is a problem,
				you can increase encut a bit to make sure all the necessary
				plane waves are included.
		"""
		wf.check_c_projectors()
		if encut == None:
			encut = 4 * wf.encut
		super(MomentumMatrix, self).__init__(wf, encut)

	@property
	def momentum_grid(self):
		"""
		Get the G-point grid used by get_momentum_matrix_elems
		and get_reciproal_fullfw. See __init__ docs for detals.
		"""
		grid = self._get_ggrid()
		return grid.reshape((grid.shape[0]//3, 3))

	def get_momentum_matrix_elems(self, b1, k1, s1, b2, k2, s2):
		"""
		< b1, k1, s1 | exp(i * (G + k1 - k2) * r) | b2, k2, s2 >
		for each G in self.momentum_grid

		Args:
			b1, k1, s1: band number, k-point, and spin indices for
				the bra band (0-indexed).
			b2, k2, s2: band number, k-point, and spin indices for
				the ket band (0-indexed).
		Returns:
			numpy array of shape (self.momentum_grid.shape[0],),
			containing the above matrix elements in the same
			order as self.momentum_grid.
		"""
		self.wf.check_bks_spec(b1, k1, s1)
		self.wf.check_bks_spec(b2, k2, s2)
		return self._get_momentum_matrix_elems(b1, k1, s1, b2, k2, s2)

	def get_reciprocal_fullfw(self, b, k, s):
		"""
		Calculates C(b,k,s,G) such that
		| b, k, s > = 1/sqrt(V) \sum_G C(b,k,s,G) exp(i * (k + G) * r).
		If the basis set is complete, \sum_G C(b,k,s,G) = 1.

		Args:
			b, k, s: band number, k-point, and spin indices (0-indexed).
		Returns:
			numpy array of shape (self.momentum_grid.shape[0],),
			containing C(b,k,s,G) in the same
			order as self.momentum_grid.
		"""
		self.wf.check_bks_spec(b, k, s)
		return self._get_reciprocal_fullfw(b, k, s)

	def g_from_wf(self, b1, k1, s1, b2, k2, s2, G):
		"""
		This function is very inefficient and is primarily for testing
		purposes. It first projects psi1 and psi2 into reciprocal space,
		then takes the dot product of the coefficients to calculate
		< b1, k1, s1 | exp(i * (G + k1 - k2) * r) | b2, k2, s2 >.
		"""
		self.wf.check_bks_spec(b1, k1, s1)
		self.wf.check_bks_spec(b2, k2, s2)
		return self._get_g_from_fullfw(b1, k1, s1, b2, k2, s2, G)