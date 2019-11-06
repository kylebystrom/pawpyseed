# cython : language_level=3
from libc.stdio cimport FILE


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
    

cdef extern from "projector.h":

    cdef ppot_t* get_projector_list(int num_els, int* labels, int* ls, double* wave_grids,
        double* projectors, double* aewaves, double* pswaves, double* rmaxs, double grid_encut)
    cdef real_proj_site_t* projector_values(int num_sites, int* labels, double* coords,
        double* lattice, double* reclattice, ppot_t* pps, int* fftg)
    cdef real_proj_site_t* smooth_pw_values(int num_N, int* Nlst, int* labels, double* coords,
        double* lattice, double* reclattice, ppot_t* pps, int* fftg)
    cdef void onto_projector_helper(band_t* band, double complex* x, real_proj_site_t* sites,
        int num_sites, double* lattice, double* reclattice, double* kpt, int num_cart_gridpts,
        int* fftg, projection_t* projections)
    cdef void get_aug_freqs_helper(band_t* band, double complex* x, real_proj_site_t* sites,
        int num_sites, double* lattice, double* reclattice, double* kpt, int num_cart_gridpts,
        int* fftg, projection_t* projections)
    cdef void onto_projector(kpoint_t* kpt, int band_num, real_proj_site_t* sites, int num_sites,
        int* G_bounds, double* lattice, double* reclattice, int num_cart_gridpts, int* fftg)
    cdef void onto_projector_ncl(kpoint_t* kpt, int band_num, real_proj_site_t* sites, int num_sites,
        int* G_bounds, double* lattice, double* reclattice, int num_cart_gridpts, int* fftg)
    cdef void onto_smoothpw(kpoint_t* kpt, int band_num, real_proj_site_t* sites, int num_sites,
        int* G_bounds, double* lattice, double* reclattice, int num_cart_gridpts, int* fftg)
    cdef void get_aug_freqs(kpoint_t* kpt, int band_num, real_proj_site_t* sites, int num_sites,
        int* G_bounds, double* lattice, double* reclattice, int num_cart_gridpts, int* fftg)
    cdef void add_num_cart_gridpts(ppot_t* pp_ptr, double* lattice, int* fftg)
    cdef void make_pwave_overlap_matrices(ppot_t* pp_ptr)
    cdef void setup_projections(pswf_t* wf, ppot_t* pps, int num_elems,
        int num_sites, int* fftg, int* labels, double* coords)
    cdef void overlap_setup_real(pswf_t* wf_R, pswf_t* wf_S,
        int* labels_R, int* labels_S, double* coords_R, double* coords_S,
        int* N_R, int* N_S, int* N_RS_R, int* N_RS_S, int num_N_R, int num_N_S, int num_N_RS)
    cdef void overlap_setup_recip(pswf_t* wf_R, pswf_t* wf_S,
        int* labels_R, int* labels_S, double* coords_R, double* coords_S,
        int* N_R, int* N_S, int* N_RS_R, int* N_RS_S, int num_N_R, int num_N_S, int num_N_RS)
    cdef void compensation_terms(double complex* overlap, int BAND_NUM, pswf_t* wf_S, pswf_t* wf_R,
        int num_M, int num_N_R, int num_N_S, int num_N_RS,
        int* M_R, int* M_S, int* N_R, int* N_S, int* N_RS_R, int* N_RS_S,
        int* proj_labels, double* proj_coords, int* ref_labels, double* ref_coords,
        int* fft_grid, int spin_flip)
    cdef void compensation_terms_recip(double complex* overlap, int BAND_NUM, pswf_t* wf_S, pswf_t* wf_R,
        int num_M, int num_N_R, int num_N_S, int num_N_RS,
        int* M_R, int* M_S, int* N_R, int* N_S, int* N_RS_R, int* N_RS_S,
        int* proj_labels, double* proj_coords, int* ref_labels, double* ref_coords,
        int* fft_grid, int spin_flip)
    cdef double* besselt(double* r, double* k, double* f, double encut, int N, int l)
    

cdef extern from "pseudoprojector.h":

    cdef void vc_pseudoprojection(pswf_t* wf_ref, pswf_t* wf_proj, int BAND_NUM, double* results)
    cdef void pseudoprojection(double complex* projections, pswf_t* wf_ref, pswf_t* wf_proj, int BAND_NUM,
                            int flip_spin)
    

cdef extern from "reader.h":

    ctypedef struct  WAVECAR:
        int type
        FILE* fp
        char* start
        char* curr
    cdef WAVECAR* wcopen(char* f, int type)
    cdef void wcseek(WAVECAR* wc, long loc)
    cdef void wcread(void* ptr0, long size, long nmemb, WAVECAR* wc)
    cdef void wcclose(WAVECAR* wc)
    cdef void setup(int nspin, int nwk, int nband,
        double* nb1, double* nb2, double* nb3, int* np, double ecut,
        double* lattice, double* reclattice)
    cdef pswf_t* read_wavecar(WAVECAR* wc, double* kpt_weights)
    cdef pswf_t* read_wavefunctions(char* filename, double* kpt_weights)
    cdef pswf_t* read_wavefunctions_from_str(char* start, double* kpt_weights)
    cdef kpoint_t** read_one_band(int* G_bounds, double* kpt_weights, int* ns, int* nk, int* nb, int BAND_NUM, char* filename)
    

cdef extern from "density.h":

    cdef void realspace_state(double complex* x, int BAND_NUM, int KPOINT_NUM,
        pswf_t* wf, int* fftg, int* labels, double* coords)
    cdef void remove_phase(double complex* x, int KPOINT_NUM, pswf_t* wf, int* fftg)
    cdef void ae_state_density(double* P, int BAND_NUM, int KPOINT_NUM, pswf_t* wf,
        int* fftg, int* labels, double* coords)
    cdef void ncl_realspace_state(double complex* x, int BAND_NUM, int KPOINT_NUM,
        pswf_t* wf, int* fftg, int* labels, double* coords)
    cdef void ae_chg_density(double* P, pswf_t* wf, int* fftg, int* labels, double* coords)
    cdef void ncl_ae_chg_density(double* P, pswf_t* wf, int* fftg, int* labels, double* coords)
    cdef void project_realspace_state(double complex* projs, int BAND_NUM, pswf_t* wf, pswf_t* wf_R,
        int* fftg, int* labels, double* coords, int* labels_R, double* coords_R)
    cdef void write_realspace_state_ncl_ri(char* filename1, char* filename2,
        char* filename3, char* filename4, int BAND_NUM, int KPOINT_NUM,
        pswf_t* wf, int* fftg, int* labels, double* coords)
    cdef double* realspace_state_ri(int BAND_NUM, int KPOINT_NUM, pswf_t* wf, int* fftg,
            int* labels, double* coords)
    cdef void write_volumetric(char* filename, double* x, int* fftg, double scale)
    cdef double* write_realspace_state_ri_return(char* filename1, char* filename2, int BAND_NUM, int KPOINT_NUM,
        pswf_t* wf, int* fftg, int* labels, double* coords)
    cdef double* write_density_return(char* filename, pswf_t* wf,
        int* fftg, int* labels, double* coords)
    cdef void write_realspace_state_ri_noreturn(char* filename1, char* filename2, int BAND_NUM, int KPOINT_NUM,
        pswf_t* wf, int* fftg, int* labels, double* coords)
    cdef void write_density_noreturn(char* filename, pswf_t* wf,
        int* fftg, int* labels, double* coords)
    

cdef extern from "sbt.h":

    ctypedef struct  sbt_descriptor_t:
            double kmin
            double kappamin
            double rmin
            double rhomin
            double drho
            double dt
            int N
            double complex** mult_table
            double* ks
            double* rs
            int lmax
    cdef sbt_descriptor_t* spherical_bessel_transform_setup(double encut, double enbuf, int lmax, int N,
        double* r, double* ks)
    cdef double* wave_spherical_bessel_transform(sbt_descriptor_t* d, double* f, int l)
    cdef double* inverse_wave_spherical_bessel_transform(sbt_descriptor_t* d, double* f, int l)
    cdef void free_sbt_descriptor(sbt_descriptor_t* d)
    

cdef extern from "linalg.h":

    cdef void fft3d(double complex* x, int* G_bounds, double* lattice,
        double* kpt, int* Gs, float complex* Cs, int num_waves, int* fftg)
    cdef void fwd_fft3d(double complex* x, int* G_bounds, double* lattice,
        double* kpt, int* Gs, float complex* Cs, int num_waves, int* fftg)
    

cdef extern from "radial.h":

    cdef double complex offsite_wave_overlap(double* dcoord, double* r1, double* f1, double** spline1, int size1,
        double* r2, double* f2, double** spline2, int size2,
        double* lattice, int l1, int m1, int l2, int m2)
    cdef double complex reciprocal_offsite_wave_overlap(double* dcoord,
        double* k1, double* f1, double** s1, int size1,
        double* k2, double* f2, double** s2, int size2,
        double* lattice, int l1, int m1, int l2, int m2)
    

cdef extern from "momentum.h":

    ctypedef struct  transform_spline_t:
        double* transform
        double** spline
    ctypedef struct  density_ft_t:
        int n1
        int l1
        int m1
        int n2
        int l2
        int m2
        int size
        double* ks
        transform_spline_t* transforms
    ctypedef struct  density_ft_elem_t:
        int num_densities
        int total_projs
        density_ft_t* densities
    cdef void free_transform_spline_list(transform_spline_t* transforms, int num_transforms)
    cdef void free_density_ft_list(density_ft_t* densities, int total_projs)
    cdef void free_density_ft_elem_list(density_ft_elem_t* elems, int num_elems)
    cdef float complex pseudo_momentum(int* GP, int* G_bounds, double* lattice,
        int* G1s, float complex* C1s, int num_waves1,
        int* G2s, float complex* C2s, int num_waves2, int* fftg)
    cdef void mul_partial_waves(double* product, int size, double* r, double* f1, double* f2)
    cdef void make_rho(double* rho, int size, double* grid, double* aewave1, double* pswave1, double* aewave2, double* pswave2)
    cdef density_ft_t spher_transforms(int size, double* r, double* f, int l1, int m1, int l2, int m2, double encut)
    cdef double complex spher_momentum(density_ft_t densities, double* G)
    cdef density_ft_elem_t get_transforms(ppot_t pp, double encut)
    cdef density_ft_elem_t* get_all_transforms(pswf_t* wf, double encut)
    cdef double complex get_momentum_matrix_element(pswf_t* wf, int* labels, double* coords,
                                                int b1, int k1, int s1,
                                                int b2, int k2, int s2,
                                                int* GP, density_ft_elem_t* elems)
    cdef void get_momentum_matrix(double complex* matrix, int numg, int* igall,
                            pswf_t* wf, int* labels, double* coords,
                            int band1, int kpt1, int spin1,
                            int band2, int kpt2, int spin2,
                            density_ft_elem_t* transforms_list,
                            double encut)
    cdef void momentum_grid_size(pswf_t* wf, double* nb1max, double* nb2max, double* nb3max,
                            int* npmax, double encut)
    cdef int get_momentum_grid(int* igall, pswf_t* wf, double nb1max, double nb2max, double nb3max, double encut)
    cdef void grid_bounds(int* G_bounds, int* gdim, int* igall, int num_waves)
    cdef void list_to_grid_map(int* grid, int* G_bounds, int* gdim, int* igall, int num_waves)
    cdef void fill_grid(float complex* x, int* Gs, float complex* Cs, int* fftg, int numg)
    cdef void fullwf_reciprocal(double complex* Cs, int* igall, pswf_t* wf, int numg,
        int band_num, int kpt_num, int* labels, double* coords)
    cdef double complex kwave_value(double* x, double* wave, double** spline, int size,
        int l, int m, double* pos)
    cdef double complex quick_overlap(int* dG, double complex* C1s, double complex* C2s, int numg,
        int* Gs, int* gmap, int* G_bounds, int* gdim)
    