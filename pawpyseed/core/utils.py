## @package pawpyseed.core.utils
# Main utilities file for the Python portion of the code.
# This files stores 1) the PAWC ctypes module, which
# contains all of the C functions used to read and analyze
# PAW wavefunctions, 2) converter functions than transfer
# data from numpy arrays to C pointers and vice versa,
# 3) cfunc_call, which is used to conveniently call ctypes
# functions, and 4) a few other utilties employed mainly
# by the wavefunction classes.


import numpy as np
from ctypes import *

import os

c_int_p = POINTER(c_int)
c_float_p = POINTER(c_float)
c_double_p = POINTER(c_double)

MODULE_DIR = os.path.dirname(os.path.abspath(__file__))
PAWC = CDLL(os.path.join(MODULE_DIR, "pawpy.so"))

PAWC.get_energy.argtypes = [c_void_p, c_int, c_int, c_int]
PAWC.get_energy.restype = c_double

PAWC.read_wavefunctions_from_str.argtypes = (c_char_p, c_double_p)
PAWC.read_wavefunctions_from_str.restype = c_void_p

PAWC.read_wavefunctions.argtypes = (c_char_p, c_double_p)
PAWC.read_wavefunctions.restype = c_void_p

PAWC.get_projector_list.argtypes = [c_int, c_int_p, c_int_p] + [c_double_p]*6 + [c_double]
PAWC.get_projector_list.restype = c_void_p

PAWC.overlap_setup_real.argtypes = [c_void_p, c_void_p, c_int_p, c_int_p,
                c_double_p, c_double_p] + 4*[c_int_p] + 3*[c_int]
PAWC.overlap_setup_real.restype = None

PAWC.realspace_state_ri.argtypes = [c_int, c_int, c_void_p, c_void_p, c_int_p, c_int_p, c_double_p]
PAWC.write_realspace_state_ri_return.argtypes = [c_char_p, c_char_p] + PAWC.realspace_state_ri.argtypes
PAWC.write_realspace_state_ri_noreturn.argtypes = PAWC.write_realspace_state_ri_return.argtypes
PAWC.write_realspace_state_ncl_ri.argtypes = [c_char_p, c_char_p] + PAWC.write_realspace_state_ri_return.argtypes
PAWC.realspace_state_ri.restype = c_double_p
PAWC.write_realspace_state_ri_return.restype = c_double_p
PAWC.write_realspace_state_ri_noreturn.restype = None
PAWC.write_realspace_state_ncl_ri.restype = None

PAWC.write_density_return.argtypes = [c_char_p, c_void_p, c_void_p, c_int_p, c_int_p, c_double_p]
PAWC.write_density_noreturn.argtypes = PAWC.write_density_return.argtypes
PAWC.write_density_return.restype = c_double_p
PAWC.write_density_noreturn.restype = None

PAWC.setup_projections.argtypes = [c_void_p, c_void_p, c_int, c_int, c_int_p, c_int_p, c_double_p]

PAWC.project_realspace_state.argtypes = [c_int, c_int, c_void_p, c_void_p, c_void_p, c_int_p, c_int_p, c_double_p, c_int_p, c_double_p]
PAWC.project_realspace_state.restype = c_void_p

PAWC.pseudoprojection.argtypes = [c_void_p, c_void_p, c_int]
PAWC.pseudoprojection.restype = c_double_p

PAWC.compensation_terms.argtypes = [c_int] + [c_void_p]*2 + [c_int]*5 + [c_int_p]*6 \
									+ [c_int_p, c_double_p]*2 + [c_int_p]
PAWC.compensation_terms.restype = c_double_p

PAWC.expand_symm_wf.argtypes = [c_void_p, c_int, c_int_p, c_double_p,\
								c_double_p, c_double_p, c_int_p]
PAWC.expand_symm_wf.restype = c_void_p

PAWC.get_occs.argtypes = [c_void_p]
PAWC.get_nband.argtypes = [c_void_p]
PAWC.get_nwk.argtypes = [c_void_p]
PAWC.get_nspin.argtypes = [c_void_p]
PAWC.get_encut.argtypes = [c_void_p]
PAWC.is_ncl.argtypes = [c_void_p]
PAWC.get_occs.restype = c_double_p
PAWC.get_nband.restype = c_int
PAWC.get_nwk.restype = c_int
PAWC.get_nspin.restype = c_int
PAWC.get_encut.restype = c_double
PAWC.is_ncl.restype = c_int

PAWC.free_ptr.argtypes = [c_void_p]
PAWC.free_ptr.restype = None

PAWC.free_ppot_list.argtypes = [c_void_p, c_int]
PAWC.free_ppot_list.restype = None

PAWC.free_pswf.argtypes = [c_void_p]
PAWC.free_pswf.restype = None

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

def cfunc_call(func, outsize, *args):
	"""
	converts args to C types and passes them to the C
	function func
	"""
	cargs = []
	for argtype, arg in zip(func.argtypes, args):
		if type(arg) == list or type(arg) == tuple:
			arg = np.array(arg)
		if type(arg) == np.ndarray:
			if argtype == c_double_p:
				cargs.append(numpy_to_cdouble(arg))
			elif argtype == c_float_p:
				cargs.append(numpy_to_cfloat(arg))
			elif argtype == c_int_p:
				cargs.append(numpy_to_cint(arg))
			else:
				raise PAWpyError("cfunc_call: invalid array type %s" % repr(arg.dtype))
		elif argtype == c_void_p:
			cargs.append(arg)
		elif argtype == c_int:
			cargs.append(c_int(arg))
		elif argtype == c_float:
			cargs.append(c_float(arg))
		elif argtype == c_double:
			cargs.append(c_double(arg))
		elif type(arg) == str:
			cargs.append(arg.encode('utf-8'))
		else:
			raise PAWpyError("cfunc_call: unsupported argument type %s" % type(arg))
	res = func(*cargs)
	if outsize == None:
		return res
	if func.restype == c_double_p:
		return cdouble_to_numpy(res, outsize)
	elif func.restype == POINTER(c_float):
		return cfloat_to_numpy(res, outsize)
	elif func.restype == c_int_p:
		return cint_to_numpy(res, outsize)
	return res

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

def cdouble_to_numpy(arr, length):
	"""
	Convert a pointer to length doubles
	in C to a numpy array of np.float64.
	Frees the pointer.
	"""
	arr = cast(arr, c_double_p)
	newarr = np.zeros(length)
	for i in range(length):
		newarr[i] = arr[i]
	PAWC.free_ptr(arr)
	return newarr

def cfloat_to_numpy(arr, length):
	"""
	Convert a pointer to length floats
	in C to a numpy array of np.float64.
	Frees the pointer.
	"""
	arr = cast(arr, c_int_p)
	newarr = np.zeros(length)
	for i in range(length):
		newarr[i] = arr[i]
	PAWC.free_ptr(arr)
	return newarr

def numpy_to_cdouble(arr):
	"""
	Convert a numpy array to
	a C double array.
	"""
	newarr = (c_double * len(arr))()
	for i in range(len(arr)):
		newarr[i] = arr[i]
	newarr = cast(newarr, c_double_p)
	return newarr

def numpy_to_cfloat(arr):
	"""
	Convert a numpy array to
	a C float array
	"""
	newarr = (c_float * len(arr))()
	for i in range(len(arr)):
		newarr[i] = arr[i]
	newarr = cast(newarr, POINTER(c_float))
	return newarr

def numpy_to_cint(arr):
	"""
	Convert a numpy array to
	a C int array.
	Casts each element to int
	"""
	newarr = (c_int * len(arr))()
	for i in range(len(arr)):
		newarr[i] = int(arr[i])
	newarr = cast(newarr, c_int_p)
	return newarr

def el(site):
	"""
	Return the element symbol of a pymatgen
	site object
	"""
	return site.specie.symbol
