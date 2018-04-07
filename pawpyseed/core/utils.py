
"""
\file
Main utilities file for the Python portion of the code.
This files stores 1) the PAWC ctypes module, which
contains all of the C functions used to read and analyze
PAW wavefunctions, 2) converter functions than transfer
data from numpy arrays to C pointers and vice versa,
3) cfunc_call, which is used to conveniently call ctypes
functions, and 4) a few other utilties employed mainly
by the wavefunction classes.
"""

import numpy as np
from ctypes import *

import os

MODULE_DIR = os.path.dirname(os.path.abspath(__file__))
PAWC = CDLL(os.path.join(MODULE_DIR, "pawpy.so"))

PAWC.read_wavefunctions.argtypes = (c_char_p, POINTER(c_double))
PAWC.read_wavefunctions.restype = c_void_p

PAWC.get_projector_list.argtypes = [c_int, POINTER(c_int), POINTER(c_int)] + [POINTER(c_double)]*6
PAWC.get_projector_list.restype = c_void_p

PAWC.overlap_setup.argtypes = [c_void_p, c_void_p, c_void_p, POINTER(c_int), POINTER(c_int),
				POINTER(c_double), POINTER(c_double)] + 4*[POINTER(c_int)] + 3*[c_int]
PAWC.overlap_setup.restype = None

PAWC.pseudoprojection.restype = POINTER(c_double)
PAWC.compensation_terms.restype = POINTER(c_double)
PAWC.get_occs.restype = POINTER(c_double)
PAWC.get_nband.restype = c_int
PAWC.get_nwk.restype = c_int
PAWC.get_nspin.restype = c_int

PAWC.free_ptr.restype = None
PAWC.free_ppot_list.restype = None
PAWC.free_pswf.restype = None

class PAWpyError(Exception):
	"""
	Class for handling errors that occur during execution
	of Python functions in pawpyseed
	"""
	def __init__(self, msg):
		self.msg = msg

def cfunc_call(func, outsize = None, *args):
	"""
	converts args to C types and passes them to the C
	function func
	"""
	cargs = []
	for arg in args:
		if type(arg) == list:
			arg = np.array(arg)
		if type(arg) == np.ndarray:
			if arg.dtype == np.float64:
				cargs.append(numpy_to_cdouble(arg))
			elif arg.dtype == np.float32:
				cargs.append(numpy_to_cfloat(arg))
			elif arg.dtype == np.int32 or arg.dtype == int:
				cargs.append(numpy_to_cint(arg))
			else:
				raise PAWpyError("cfunc_call: invalid array type %s" % repr(arg.dtype))
		elif type(arg) == int:
			cargs.append(arg)
		elif type(arg) == float:
			cargs.append(arg)
		elif type(arg) == str:
			cargs.append(arg.encode('utf-8'))
	res = func(*cargs)
	if outsize == None:
		return res
	if func.restype == POINTER(c_double):
		return cdouble_to_numpy(res, outsize)
	elif func.restype == POINTER(c_float):
		return cfloat_to_numpy(res, outsize)
	elif func.restype == POINTER(c_int):
		return cint_to_numpy(res, outsize)
	return res

def check_spin(spin, nspin):
	"""
	Utility to check if the spin input parameter to single_band_projection
	and similar functions is allowed given nspin of the wavefunction object
	being analyzed. Returns a new value of spin if spin must be changed,
	raises an error is spin is not allowed.
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
	arr = cast(arr, POINTER(c_double))
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
	arr = cast(arr, POINTER(c_int))
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
	newarr = cast(newarr, POINTER(c_double))
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
	newarr = cast(newarr, POINTER(c_int))
	return newarr

def el(site):
	"""
	Return the element symbol of a pymatgen
	site object
	"""
	return site.specie.symbol
