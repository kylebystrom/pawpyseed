#ifndef HF_H
#define HF_H

typedef struct radial_set {
	int l; ///< angular momentum quantum number
	int v; ///< valence level
	int X; ///< number of basis functions/wavefunction solutions per l
	int* occs; ///< (X) number of electrons in a shell
	double* h; ///< single particle hamiltonian (XxX)
	double*** ee; ///< electron-electron repulsion terms 4x(XxX)x(XxX) 0 2 4 6
	double* es; ///< (X)
	double* Ps; ///< (XxX)
	double*** yks; ///< (XxX)x(4)x(N) 0 2 4 6
	double** bfs;
	double* DM; ///< density matrix (XxX)
	double E; ///< total energy
} radial_set_t;

typedef struct awf {
	int Z; ///< atomic number/number of electrons
	int L; ///< number of l quantums numbers, ie lmax+1
	int X; ///< number of basis functions/wavefunction solutions per l
	int XT; ///< L*X
	int N; ///< grid size
	double*** yks; ///< (XTxXT)x(7)x(N), k = 0 1 2 3 4 5 6, indexed by (l1*X+n1)*XT+(l2*X_n2),l2,k,rindex
	double***** J; ///< (7)x(L)x(L)x(XxX)x(XxX) indexed by l1,l2,n1a*X+n1b,n2a*X+n2b <n1al1,n2al2||n1bl1,n2bl2>
	double***** K; ///< (7)x(L)x(L)x(XxX)x(XxX) indexed by l1,l2,n1a*X+n1b,n2a*X+n2b <n1al1,n2al2||n2bl2,n1bl1>
	radial_set_t* wfs; ///< L radial_set_t wavefunctions
	double* r;
	double E; ///< total energy
} awf_t;

double get_yk(double*** yks, int X, int XT, int l1, int n1, int l2, int n2, int k, int rindex);

void set_yk(double*** yks, double* nums, int X, int XT, int l1, int n1, int l2, int n2, int k);

double get_coul(double***** J_or_K, int k, int X, int l1, int n1a, int n1b, int l2, int n2a, int n2b);

void set_coul(double***** J_or_K, double num, int k, int X, int l1, int n1a, int n1b, int l2, int n2a, int n2b);

awf_t* construct_basis(int Z, int N, int maxN, int maxL, double* r);

awf_t* setup_H(int N, int maxN, int maxL, double* r);

void assign_occs(awf_t* wf);

awf_t* setup(int Z, int N, int maxN, int maxL, double* r, double** P0s);

void make_density_matrix(radial_set_t* wf);

void calc_energy(awf_t* wf);

void solve(awf_t* wf, int maxsteps);

double get_E(awf_t* wf) {
	return wf->E;
}

double* get_occs(awf_t* wf, int l) {
	return wf->wfs[l].occs;
}

double* get_Ps(awf_t* wf, int l) {
	return wf->wfs[l].Ps;
}

#endif

