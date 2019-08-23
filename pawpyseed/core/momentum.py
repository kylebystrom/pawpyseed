from pawpyseed.core.wavefunction import Wavefunction
from pawpyseed.core import pawpyc

class MomentumMatrix(pawpyc.CMomentumMatrix):

	def __init__(self, wf, encut = None):
		wf.check_c_projectors()
		if encut == None:
			encut = 4 * wf.encut
		super(MomentumMatrix, self).__init__(wf, encut)

	@property
	def momentum_grid(self):
		grid = self._get_ggrid()
		return grid.reshape((grid.shape[0]//3, 3))

	def get_momentum_matrix_elems(self, b1, k1, s1, b2, k2, s2):
		return self._get_momentum_matrix_elems(b1, k1, s1, b2, k2, s2)

	def get_reciprocal_fullfw(self, b, k, s):
		return self._get_reciprocal_fullfw(b, k, s)

	def g_from_wf(self, b1, k1, s1, b2, k2, s2, G):
		return self._get_g_from_fullfw(b1,k1,s1,b2,k2,s2,G)