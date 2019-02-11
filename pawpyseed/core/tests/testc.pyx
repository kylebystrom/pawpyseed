# cython : profile=True
# cython : language_level=3

from pawpyseed.core.tests cimport testc_extern as tc
from pawpyseed.core cimport pawpyc

from cpython cimport array
from libc.stdlib cimport malloc, free
from libc.stdio cimport FILE
from pymatgen.core.structure import Structure
from pymatgen.io.vasp.outputs import Vasprun
from monty.io import zopen
import numpy as np
from numpy.testing import assert_almost_equal
cimport numpy as np
import time
import sys
from libc.stdint cimport uintptr_t
from pawpyseed.core.symmetry import *

cpdef fft_check(str wavecar, np.ndarray[double, ndim=1] kpt_weights,
	np.ndarray[int, ndim=1] fftgrid):

	cdef double[::1] kws = kpt_weights
	cdef int[::1] dimv = fftgrid
	return tc.fft_check(wavecar.encode('utf-8'), &kws[0], &dimv[0]);

cpdef proj_check(pawpyc.CWavefunction wf):
	for b in range(wf.nband):
		for k in range(wf.nwk * wf.nspin):
			print("testing %d %d" % (b, k))
			sys.stdout.flush()
			tc.proj_check(b, k, wf.wf_ptr, &wf.dimv[0],
				&wf.nums[0], &wf.coords[0])
	print("FINISHED PROJ CHECK")