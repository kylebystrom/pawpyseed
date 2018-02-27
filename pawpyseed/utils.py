import numpy as np
from ctypes import *

import os

MODULE_DIR = os.path.dirname(os.path.abspath(__file__))
PAWC = CDLL(os.path.join(MODULE_DIR, "pawpy.so"))

PAWC.read_wavefunctions.restype = POINTER(None)
PAWC.get_projector_list.restype = POINTER(None)
PAWC.read_wavefunctions.restype = POINTER(None)
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
	return newarr

def numpy_to_cfloat(arr):
	"""
	Convert a numpy array to
	a C float array
	"""
	newarr = (c_float * len(arr))()
	for i in range(len(arr)):
		newarr[i] = arr[i]
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
	return newarr

def el(site):
	"""
	Return the element symbol of a pymatgen
	site object
	"""
	return site.specie.symbol
