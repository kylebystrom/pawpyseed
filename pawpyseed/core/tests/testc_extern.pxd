# cython : language_level=3
from libc.stdio cimport FILE


cdef extern from "tests/tests.h":

    cdef int fft_check(char* wavecar, double* kpt_weights, int* fftg)
    cdef void proj_check(int BAND_NUM, int KPOINT_NUM,
        pswf_t* wf, int* fftg, int* labels, double* coords)
    

cdef extern from "utils.h":

    ctypedef struct  funcset_t:
        int l
        double* proj
        double** proj_spline
        double* aewave
        double** aewave_spline
        double* pswave
        double** pswave_spline
        double* diffwave
        double** diffwave_spline
        double* kwave
        double** kwave_spline
        double* smooth_diffwave
        double** smooth_diffwave_spline
        double* dense_kwave
        double** dense_kwave_spline
    ctypedef struct  ppot_t:
        int num_projs
        int total_projs
        int lmax
        funcset_t* funcs
        double rmax
        double wave_rmax
        double* pspw_overlap_matrix
        double* aepw_overlap_matrix
        double* diff_overlap_matrix
        int proj_gridsize
        int wave_gridsize
        int num_cart_gridpts
        double* wave_grid
        double* kwave_grid
        double* proj_grid
        double* smooth_grid
        double* dense_kgrid
    ctypedef struct  projection_t:
        int num_projs
        int total_projs
        int* ns
        int* ls
        int* ms
        double complex* overlaps
    ctypedef struct  band_t:
        int n
        int num_waves
        double occ
        double N
        double energy
        float complex* Cs
        double complex* CRs
        float complex* CAs
        projection_t* projections
        projection_t* up_projections
        projection_t* down_projections
        projection_t* wave_projections
    ctypedef struct  rayleigh_set_t:
        int l
        double complex* terms
    ctypedef struct  kpoint_t:
        short int up
        int num_waves
        int* Gs
        double* k
        double weight
        int num_bands
        band_t** bands
        rayleigh_set_t** expansion
    ctypedef struct  pswf_t:
        double encut
        int num_elems
        int* num_projs
        int num_sites
        ppot_t* pps
        int* G_bounds
        kpoint_t** kpts
        int nspin
        int nband
        int nwk
        double* lattice
        double* reclattice
        int* fftg
        int is_ncl
        int wp_num
        int num_aug_overlap_sites
        double* dcoords
        double complex** overlaps
    ctypedef struct  projgrid_t:
        double complex* values
    ctypedef struct  real_proj_t:
        int l
        int m
        int func_num
        double complex* values
    ctypedef struct  real_proj_site_t:
        int index
        int elem
        int num_projs
        int total_projs
        int num_indices
        int gridsize
        double rmax
        double* coord
        int* indices
        double* paths
        real_proj_t* projs
    cdef void affine_transform(double* out, double* op, double* inv)
    cdef void rotation_transform(double* out, double* op, double* inv)
    cdef int min(int a, int b)
    cdef int max(int a, int b)
    cdef void vcross(double* res, double* top, double* bottom)
    cdef double dot(double* x1, double* x2)
    cdef double mag(double* x1)
    cdef double determinant(double* m)
    cdef double dist_from_frac(double* coords1, double* coords2, double* lattice)
    cdef void frac_to_cartesian(double* coord, double* lattice)
    cdef void cartesian_to_frac(double* coord, double* reclattice)
    cdef void min_cart_path(double* coord, double* center, double* lattice, double* path, double* r)
    cdef void trilinear_interpolate_values(double complex* x, double* frac, int* fftg, double complex* values)
    cdef double complex trilinear_interpolate(double complex* c, double* frac, int* fftg)
    cdef void free_projection_list(projection_t* projlist, int num)
    cdef void clean_wave_projections(pswf_t* wf)
    cdef void free_kpoint(kpoint_t* kpt, int num_elems, int num_sites, int wp_num, int* num_projs)
    cdef void free_ppot(ppot_t* pp)
    cdef void free_real_proj(real_proj_t* proj)
    cdef void free_real_proj_site(real_proj_site_t* site)
    cdef void free_pswf(pswf_t* wf)
    cdef void free_ptr(void* ptr)
    cdef void free_real_proj_site_list(real_proj_site_t* sites, int length)
    cdef void free_ppot_list(ppot_t* pps, int length)
    cdef double* get_occs(pswf_t* wf)
    cdef int get_nband(pswf_t* wf)
    cdef int get_nwk(pswf_t* wf)
    cdef int get_nspin(pswf_t* wf)
    cdef int is_ncl(pswf_t* wf)
    cdef double get_encut(pswf_t* wf)
    cdef double get_energy(pswf_t* wf, int band, int kpt, int spin)
    cdef double get_occ(pswf_t* wf, int band, int kpt, int spin)
    cdef void set_num_sites(pswf_t* wf, int nsites)
    cdef double legendre(int l, int m, double x)
    cdef void legendre_coeff(double* ptr, int l, int m)
    cdef double* legendre_product(int l1, int l2, int m1, int m2)
    cdef double fac(int n)
    cdef double complex Ylm(int l, int m, double theta, double phi)
    cdef double complex Ylm2(int l, int m, double costheta, double phi)
    cdef double proj_interpolate(double r, double rmax, int size, double* x,
        double* proj, double** proj_spline)
    cdef double wave_interpolate(double r, int size, double* x, double* f,
        double** wave_spline)
    cdef double complex proj_value_helper(double r, double rmax, int size,
        double* pos, double* x, double* f, double** s, int l, int m)
    cdef double complex proj_value(funcset_t funcs, double* x, int m, double rmax,
        int size, double* ion_pos, double* pos, double* lattice)
    cdef double complex smooth_wave_value(funcset_t funcs, double* x, int m, double rmax,
        int size, double* ion_pos, double* pos, double* lattice)
    cdef double complex wave_value(funcset_t funcs, int size, double* x, int m,
            double* ion_pos, double* pos, double* lattice)
    cdef double complex wave_value2(double* x, double* wave, double** spline, int size,
        int l, int m, double* pos)
    cdef void setup_site(real_proj_site_t* sites, ppot_t* pps, int num_sites, int* site_nums,
        int* labels, double* coords, double* lattice, int* fftg, int pr0_pw1)
    cdef double** spline_coeff(double* x, double* y, int N)
    cdef double spline_integral(double* x, double* a, double** s, int size)
    cdef void frac_from_index(int index, double* coord, int* fftg)
    cdef double sph_bessel(double k, double r, int l)
    cdef double sbf(double x, int l)
    cdef pswf_t* expand_symm_wf(pswf_t* rwf, int num_kpts, int* maps,
        double* ops, double* drs, double* kws, int* trs)
    cdef void CHECK_ALLOCATION(void* ptr)
    cdef void ALLOCATION_FAILED()
    cdef void CHECK_STATUS(int status)
    