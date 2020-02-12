# cython : language_level=3

from pawpyseed.core cimport pawpyc_extern as ppc
cimport numpy as np 
import numpy as np

cdef class PWFPointer:

    cdef ppc.pswf_t* ptr
    cdef readonly np.ndarray kpts
    cdef readonly np.ndarray weights
    cdef readonly np.ndarray band_props

    @staticmethod
    cdef PWFPointer from_pointer_and_kpts(ppc.pswf_t* ptr,
        structure, kpts, band_props, allkpts, weights, symprec,
        time_reversal_symmetry)

cdef class PseudoWavefunction:

    cdef ppc.pswf_t* wf_ptr
    cdef readonly int nband
    cdef readonly int nwk
    cdef readonly int nspin
    cdef readonly int ncl
    cdef readonly np.ndarray kws
    cdef readonly np.ndarray kpts

cdef class CWavefunction(PseudoWavefunction):

    cdef int[::1] dimv
    cdef int[::1] fdimv
    cdef int gridsize

    cdef int[::1] nums
    cdef double[::1] coords
    cdef int number_projector_elements
    cdef readonly int projector_owner

cdef class CProjector:

    cdef public CWavefunction wf
    cdef public CWavefunction basis

    cdef int[::1] M_R
    cdef int[::1] M_S
    cdef int[::1] N_R
    cdef int[::1] N_S
    cdef int[::1] N_RS_R
    cdef int[::1] N_RS_S

    cdef int num_M_R
    cdef int num_M_S 
    cdef int num_N_R
    cdef int num_N_S
    cdef int num_N_RS_R
    cdef int num_N_RS_S

cdef class CMomentumMatrix:

    cdef public CWavefunction wf
    cdef np.ndarray ggrid
    cdef ppc.density_ft_elem_t* elem_density_transforms
    cdef public double momentum_encut
    cdef np.ndarray gdim
    cdef np.ndarray gbounds
    cdef np.ndarray grid3d
