#include <math.h>
#include <mkl.h>
#include "condon_shortley.h"

typedef radial_set {
	int X; ///< number of radial functions;
	double** Ps; ///< XxN
	double** Pinvrs; ///< XxN
} radial_set_t;

typedef struct radial_set {
	int l; ///< angular momentum quantum number
	int v; ///< valence level
	int X; ///< number of basis functions/wavefunction solutions per l
	int* occs; ///< (X) number of electrons in a shell
	double* h; ///< single particle hamiltonian (XxX)
	double*** ee; ///< electron-electron repulsion terms 4x(XxX)x(XxX) 0 2 4 6
	double* es; ///< (X)
	double** Ps; ///< (L)x(XxX)
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

double laguerre(double x, int n, int alpha) {
	double l0 = 1;
	if (n == 0) return l0;
	double l1 = 1 + alpha - x;
	double lcurr = l1;
	for (int k = 1; k < n; k++) {
		l1 = lcurr;
		lcurr = ((2*k+1+alpha-x) * l1 - (k+alpha) * l0) / (k+1);
		l0 = l1;
	}
}

double hradial() {
	pow(4.0/n/n/n * fac(n-l-1)/fac(n+l), 0.5) * exp(-r/2/n) * pow(r/n, l), * laguerre(r, n-l-1, 2*l+1);
}

double get_yk(double*** yks, int X, int XT, int l1, int n1, int l2, int n2, int k, int rindex) {
	return yks[(l1*X+n1)*XT+(l2*X_n2)][k][rindex];
}

void set_yk(double*** yks, double* nums, int X, int XT, int l1, int n1, int l2, int n2, int k) {
	yks[(l1*X+n1)*XT+(l2*X_n2)][k][rindex] = nums;
}

double get_coul(double**** J_or_K, int k, int X, int l1, int n1a, int n1b, int l2, int n2a, int n2b) {
	return J_or_K[k][l1][l2][n1a*X+n1b][n2a*X+n2b];
}

void set_coul(double**** J_or_K, double num, int k, int X, int l1, int n1a, int n1b, int l2, int n2a, int n2b) {
	J_or_K[k][l1][l2][n1a*X+n1b][n2a*X+n2b] =  num;
}

double* discrete_sp_hamiltonian(int Z, int l, int size, double* r) {
	double* h = (double*) calloc(size * size * sizeof(double));

	for (int i = 0; i < N - 1; i++) {
		dR = r[i+1] - r[i];
		h[N*i+i] = 1/(r[i+1]-r[i]) + 1/(r[i]-r[i-1])
					- 2*Z/r[i] + l*(l+1)/r[i]/r[i]
					- 2/r[i]*yks[i*X+i];
		h[N*i+i+1] = -1/(r[i+1]-r[i]);
		h[N*(i+1)+i] = h[N*i+i+1];
		h[N*i+i-1] = -1/(r[i]-r[i-1]);
		h[N*(i-1)+i] = h[N*i+i-1];
	}
}

awf_t* construct_basis(int Z, int N, int maxN, int maxL, double* r) {
	int X = maxN;
	int L = maxL + 1;
	int XT = X * L;
	radial_set_t* wfs = (radial_set_t*) malloc(L * sizeof(radial_set_t));

	for (int l = 0; l < L; l++) {
		double* h = (double*) calloc(X * X, sizeof(double));
		int* ls = (int*) malloc(X * sizeof(int));
		int* occs = (int*) calloc(X, sizeof(int));
		double** bfs = (double**) malloc(X * sizeof(double*));
		double*** splines = (double***) malloc(X * sizeof(double**));
		for (int n = 1; n <= maxN; n++) {
			double* bf = (double*) malloc(N * sizeof(double*));
			for (int j = 0; j < N; j++) {
				bf[j] = r[j] * hradial(n, l, r[j]);
			}
			bfs[n-1] = bf;
			splines[n-1] = spline_coeff(r, bf, N);
			h[X*n-X+n-1] = -0.5*Z*Z/(n+l)/(n+l);
		}
		wfs[l].bfs = bfs;
		wfs[l].l = l;
		wfs[l].v = 0;
		wfs[l].X = X;
		wfs[l].h = h;
		wfs[l].occs = occs;
	}
	double*** yks = (double***) malloc(XT*XT * sizeof(double**));
	int b1, b2, n1, n2, l1, l2;
	for (int l1 = 0; l1 < L; l1++) {
		for (int l2 = 0; l2 < L; l2++) {
			for (int n1 = 0; n1 < X; n1++) {
				for (int n2 = 0; n2 < X; n2++) {
					yks[i] = (double**) malloc(7 * sizeof(double*));
					for (int k = 0; k <= 6; k++) {
						set_yk(yks, yk(k, N, r, wfs[l1].bfs[n1], wfs[l2].bfs[n2]),
							X, XT, l1, n1, l2, n2, k);
					}
				}
			}
		}
	}
	double***** J = (double*****) malloc(7 * sizeof(double****));
	double***** K = (double*****) malloc(7 * sizeof(double****));
	for (int k = 0; k < 7; k++) {
		J[k] = (double****) malloc(L * sizeof(double***));
		K[k] = (double****) malloc(L * sizeof(double***));
		for (int l1 = 0; l1 < L; l1++) {
			J[k][l1] = (double***) malloc(L * sizeof(double**));
			K[k][l1] = (double***) malloc(L * sizeof(double**));
			for (int l2 = 0; l2 < L; l2++) {
				J[k][l1][l2] = (double**) malloc(X * X * sizeof(double*));
				K[k][l1][l2] = (double**) malloc(X * X * sizeof(double*));
				for (int n1a = 0; n1a < X; n1a++) {
					for (int n1b = 0; n1b < X; n1b++) {
						J[k][l1][l2][n1a*X+n1b] = (double*) malloc(X * X * sizeof(double));
						K[k][l1][l2][n1a*X+n1b] = (double*) malloc(X * X * sizeof(double));
						for (int n2a = 0; n2a < X; n2a++) {
							for (int n2b = 0; n2b < X; n2b++) {
								double* integrandJ = (double*) malloc(N * sizeof(double));
								double* integrandK = (double*) malloc(N * sizeof(double));
								for (int rnum = 0; rnum < N; rnum++) {
									integrandJ[rnum] = get_yk(yks, X, XT, l2, n2a, l2, n2b, k, rnum)
										* wfs[l1].bfs[n1a] * wfs[l1].bfs[n1a] / r[rnum];
									integrandK[rnum] = get_yk(yks, X, XT, l2, n2a, l1, n1b, k, rnum)
										* wfs[l1].bfs[n1a] * wfs[l2].bfs[n2b] / r[rnum];
								}
								double** splineJ = spline_coeff(r, integrandJ, N);
								set_coul(J, spline_integral(r, integrandJ, splineJ, N),
									k, X, XT, l1, n1a, n1b, l2, n2a, n2b);
								double** splineK = spline_coeff(r, integrandK, N);
								set_coul(K, spline_integral(r, integrandK, splineK, N)
									* pow((4*l1+2)*(4*l2+2), -0.5) * cg_coeff[k][l1][l2][l1][l2],
									k, X, XT, l1, n1a, n1b, l2, n2a, n2b);
								free(splineJ[0]);
								free(splineJ[1]);
								free(splineJ[2]);
								free(splineJ);
								free(integrandJ);
								free(splineK[0]);
								free(splineK[1]);
								free(splineK[2]);
								free(splineK);
								free(integrandK);
							}
						}
					}
				}
			}
		}
	}

	awf_t* wf = (awf_t*) malloc(sizeof(awf_t));
	wf->wfs = wfs;
	wf->L = L;
	wf->X = X;
	wf->N = N;
	wf->J = J;
	wf->K = K;
	wf->Z = 1;
	wf->yks = yks;
	wf->r = r;
	double E = 0;

	return wf;
}

awf_t* setup_H(int N, int maxN, int maxL, double* r) {
	awf_t* wf = construct_basis(1, N, maxN, maxL, r);
	for (int l = 0; l < wf->L; l++) {
		radial_set_t* rwf = wf->wfs[l];
		int X = rwf->X;
		for (int b = 0; b < X; b1++) {
			rwf->Ps[b*X+b] = 1.0;
		}
	}
	wf->wfs[0].occs[0] = 1;
	wf->Z = 1;
	return wf;
}

void assign_occs(awf_t* wf) {
	for (int l = 0; l < wf->L; l++) {
		for (int b = 0; b < wf->wfs[l].X; b++) {
			wf->wfs[l].occs[b] = 0;
		}
		wf->wfs[l].v = 0;
	}
	int reme = wf->Z;
	int nume = 0;
	int minl = 0;
	double teste = 0;
	double mine = 10000;
	while (reme > 0) {
		for (int l = 0; l < wf->L; l++) {
			teste = wf->wfs[l].es[wf->wfs[l].v];
			if (teste < mine) {
				mine = teste;
				minl = l;
			}
		}
		if (reme > 4*minl+2) {
			reme -= 4*minl+2;
			wf->wfs[minl].occs[wf->wfs[minl].v] = 4*minl+2;
			wf->wfs[minl].v++;
		} else {
			wf->wfs[minl].occs[wf->wfs[minl].v] = reme;
			wf->wfs[minl].v++;
			reme = 0;
		}
	}
}

awf_t* setup(int Z, int N, int maxN, int maxL, double* r, double** P0s) {
	awf_t* wf = construct_basis(Z, N, maxN, maxL, r);
	for (int l = 0; l < wf->L; l++) {
		wf->wfs[l].Ps = P0s[l];
	}
	wf->Z = Z;
	assign_occs(awf_t* wf);
	return wf;
}

void make_density_matrix(radial_set_t* wf) {
	int X = wf->X;
	wf->DM = (double*) malloc(X*X * sizeof(double));
	for (int i = 0; i < X; i++) {
		wf->DM[i*X+i] = wf->occs[i];
	}
	double* temp = (double*) malloc(X*X * sizeof(double));

	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, X, X, X, 1, wf->DM, X, wf->Ps, X, 0, temp, X);
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, X, X, X, 1, wf->Ps, X, temp, X, 0, wf->DM, X);
	free(temp);
}

void solve(awf_t* wf, int maxsteps) {

	int X = wf->X;
	for (int step = 0; step < maxsteps; step++) {
		for (int l = 0; l < wf->L; l++) {

			radial_set_t rwf = wf->wfs[l];

			double* hamiltonian = (double*) malloc(X * X * sizeof(double));

			double* es = rwf.es;
			/*make_density_matrix(rwf.Ps);
			for (int p = 0; p < X; P++) {
				if (wf->wfs[l].occs[p] != 0) {
					double coul = 0;
					for (int mu = 0; mu < X; mu++) {
						for (int nu = 0; nu < X; nu++) {
							coul += (1-wf->wfs[l].occs[p]) * Ps[X*mu+p] * Ps[X*nu+p]
								* get_coul(wf->J, 0, X, l, mu, nu, l, mu, nu);
							for (int k = 2; k < 7; k+=2) {
								coul -= get_coul(wf->K, k, X, l, mu, nu, l, mu, nu);
							}
						}
					}
				}
			}*/
			double J, K;
			for (int mu = 0; mu < X; mu++) {
				for (int nu = 0; nu < X; nu++) {
					hamiltonian[X*mu+nu] = rwf.h[X*mu+nu];
					for (int lj = 0; lj < wf->L; lj++) {
						double* DM = wf->wfs[lj].DM;
						double** ee = wf->wfs[l].ee[lj];
						for (int lambda = 0; lambda < X; lambda++) {
							for (int sigma = 0; sigma < X; sigma++) {
								hamiltonian[X*mu+nu] += DM[lambda*X+sigma] *
									get_coul(wf->J, 0, X, l, mu, nu, lj, lambda, sigma);
								for (int k = 0; k < 7; k++) {
									hamiltonian[X*mu+nu] -= DM[lambda*X+sigma] *
										get_coul(wf->K, k, X, l, mu, nu, lj, lambda, sigma);
								}
							}
						}
					}
					double* DM = wf->wfs[l].DM;
					double** ee = wf->wfs[l].ee[l];
					for (int lambda = 0; lambda < X; lambda++) {
						for (int sigma = 0; sigma < X; sigma++) {
							J = ee[mu*X+nu][lambda*X+sigma];
							hamiltonian[X*mu+nu] += DM[lambda*X+sigma] *
								();
						}
					}
				}
			}
			LAPACKE_dsyevd(LAPACKE_ROW_MAJOR, 'V', 'U', X, hamiltonian, X, rwf.es);
			free(rwf.Ps);
			rwf.Ps = hamiltonian;
		}
	}
}

double* yk(int k, int size, double* r, double* P1, double* P2) {
	int size1 = rnum+1;
	int size2 = size-rnum;

	double* integrand1 = (double*) malloc(size * sizeof(double));
	double* integrand2 = (double*) malloc(size * sizeof(double));

	for (int i = 0; i < size; i++) {
		integrand1[i] = P1[i] * P2[i] * pow(r[i], k-1); //might need to change to k
	}
	double** spline1 = spline_coeff(grid, integrand1, size);
	for (int i = 0; i < size; i++) {
		integrand2[i] = P1[i] * P2[i] * pow(r[i], -k-2); //might need to change to -k-1
	}
	double** spline2 = spline_coeff(grid, integrand2, size);
	integrals1[0] = r[0] * P1[0];
	integrals2[size-1] = 0;
	double* a = integrand1, b = spline1[0], c = spline1[1], d = spline1[2];
	double* a = integrand2, b = spline2[0], c = spline2[1], d = spline2[2];
	double dx=0;
	for (int i = 0; i < size-1; i++) {
		dx = r[i+1] - r[i];
		integrals1[i+1] = integral1s[i] + dx * (a[i] + dx * (b[i]/2 + dx * (c[i]/3 + d[i]*dx/4)));
		integrals1[i+1] *= pow(r[i], -k);
		j = size - i - 2;
		dx = r[j] - r[j-1];
		integrals2[j] = integrals2[j+1] + dx * (e[j] + dx * (f[j]/2 + dx * (g[i]/3 + h[j]*dx/4)));
		integrals2[j] *= power(r[j], k+1);
	}

	for (int i = 0; i < size; i++) {
		intgrals1[i] += integrals2[i];
	}

	free(integrals2);
	free(integrand1);
	free(integrand2);
	free(spline1[0]);
	free(spline1[1]);
	free(spline1[2]);
	free(spline2[0]);
	free(spline2[1]);
	free(spline2[2]);
	return integrals1;

}

double* sp_hamiltonian() {

}

void construct_hamiltonian(awf_t* wf, double* H) {
	int N = wf->N;
	int X = wf->X;
	double* r = wf->r;
	for (int i = 0; i < X; i++) {
		for (j = 0; j < N; j++) {
			H
		}
	}
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < X; j++) {
			H[N*i+i] += 2*wf->occs[j]/r[j] * yks[X*j+j][0][i];
		}
		H[N*i+i] = 1/(r[i+1]-r[i]) + 1/(r[i]-r[i-1])
					- 2*Z/r[i] + l*(l+1)/r[i]/r[i]
					- 2/r[i]*yks[i*X+i];
		H[N*i+i+1] = -1/(r[i+1]-r[i]);
		H[N*(i+1)+i] = H[N*i+i+1];
		H[N*i+i-1] = -1/(r[i]-r[i-1]);
		H[N*(i-1)+i] = H[N*i+i-1];
	}

	for (int i = 0; i < N; i++) {
		H[N*i+i] = 1/(r[i+1]-r[i]) + 1/(r[i]-r[i-1])
					- 2*Z/r[i] + l*(l+1)/r[i]/r[i]
					- 2/r[i]*ykbfs[i*X+i];
		H[N*i+i+1] = -1/(r[i+1]-r[i]);
		H[N*(i+1)+i] = H[N*i+i+1];
		H[N*i+i-1] = -1/(r[i]-r[i-1]);
		H[N*(i-1)+i] = H[N*i+i-1];
	}
}

integral += dx * (a[i] + dx * (b[i]/2 + dx * (c[i]/3 + d[i]*dx/4)));