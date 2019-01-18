## @package pawpyseed.core.utils
# Miscellaneous utilities file for the Python portion of the code..

import numpy as np
import os

class PAWpyError(Exception):
	"""
	Class for handling errors that occur during execution
	of Python functions in pawpyseed
	"""
	def __init__(self, msg):
		self.msg = msg

class PAWpyWarning(Warning):
	"""
	Class for handling warnings taht occur during execution
	of Python functions in pawpyseed
	"""
	def __init__(self, msg):
		self.msg = msg

def check_spin(spin, nspin):
	"""
	Utility to check if the spin input parameter to single_band_projection
	and similar functions is allowed given nspin of the wavefunction object
	being analyzed. Returns a new value of spin if spin must be changed,
	raises an error if spin is not allowed.
	"""
	if spin >= 0:
		if spin >= nspin:
			raise PAWpyError('spin must be less than nspin. spin is %d, nspin is %d' % (spin, nspin))
		else:
			return 1
	return nspin

def el(site):
	"""
	Return the element symbol of a pymatgen
	site object
	"""
	return site.specie.symbol
