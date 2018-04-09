
typedef struct awf {
	int Z; ///< atomic number/number of electrons
	int L; ///< number of l quantums numbers, ie lmax+1
	int X; ///< number of basis functions/wavefunction solutions per l
	int XT; ///< L*X
	int N; ///< number of electrons
	int G; ///< size of radial grid
	double*** yks; ///< (XTxXT)x(7)x(G), k = 0 1 2 3 4 5 6, indexed by (l1*X+n1)*XT+(l2*X_n2),l2,k,rindex
	double***** J; ///< (7)x(L)x(L)x(XxX)x(XxX) indexed by l1,l2,n1a*X+n1b,n2a*X+n2b <n1al1,n2al2||n1bl1,n2bl2>
	double***** K; ///< (7)x(L)x(L)x(XxX)x(XxX) indexed by l1,l2,n1a*X+n1b,n2a*X+n2b <n1al1,n2al2||n2bl2,n1bl1>
	radial_set_t* wfs; ///< L radial_set_t wavefunctions
	double* r;
	double E; ///< total energy
} awf_t;

awf_t* construct_basis(int Z, int N, int maxN, int maxL, double* r) {
	int X = maxN;
	int L = maxL + 1;
	for (int l = 0; l < L; l++) {
		XT += 2*l+1;
	}
	int XT = X * L;

	double* h = (double*) calloc(XT * XT, sizeof(double));
	int* ls = (int*) malloc(XT * sizeof(int));
	int* ms = (int*) malloc(XT * sizeof(int));
	double** bfs = (double**) malloc(X * sizeof(double*));
	double*** splines = (double***) malloc(X * sizeof(double**));
	int nq = 0, t = 0;
	for (int n = 0; n < maxN; n++) {
		for (int l = 0; l < L; l++) {
			for (int m = -l; m <= l; m++) {
				nq = n+l+1;
				double* bf = (double*) malloc(N * sizeof(double*));
				for (int j = 0; j < N; j++) {
					bf[j] = r[j] * hradial(nq, l, r[j]);
				}
				bfs[t] = bf;
				splines[t] = spline_coeff(r, bf, N);
				h[XT*t+t] = -0.5*Z*Z/nq/nq;
				t++;
			}
		}
	}
	for (int i = 0; i < XT; i++) {
		l = i / X;
		n = i % X + 1;
	}
	wf->bfs = bfs;
	wf->v = 0;
	wf->X = X;
	wf->h = h;
	
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
									k, X, l1, n1a, n1b, l2, n2a, n2b);
								double** splineK = spline_coeff(r, integrandK, N);
								set_coul(K, spline_integral(r, integrandK, splineK, N)
									* pow((4*l1+2)*(4*l2+2), -0.5) * cg_coeff[k][l1][l2][l1][l2],
									k, X, l1, n1a, n1b, l2, n2a, n2b);
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
	wf->J = J;
	wf->K = K;
	wf->Z = 1;
	wf->yks = yks;
	wf->r = r;
	double E = 0;

	return wf;
}

void make_density_matrix(awf_t* wf) {

	double* temp = (double*) calloc(X*(N), sizeof(double));
	int X = wf->X;
	for (int i = 0; i < X; i++) {
		for (int j = 0; j < N; j++) {
			temp[i*(N)+j] = wf->Ps[i*X+j];
		}
	}

	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, X, X, N, 1, wf->DM, N, wf->Ps, X, 0, wf->DM, X);
	free(temp);
}

void calc_energy(awf_t* wf) {
	wf->E = 0;
	double* temp = (double*) malloc(X*X * sizeof(double));
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, X, X, X, 1,
		wf->DM, X, wf->h, X, 0, temp, X);
	for (int i = 0; i < N; i++) {
		wf->E += temp[i*X+i];
		wf->E += wf.es[i];
	}
	free(temp);
	wf->E = wf->E * 0.5;
}

awf_t* setup_H(int G, int maxN, int maxL, double* r) {
	awf_t* wf = construct_basis(1, G, maxN, maxL, r);
	int X = wf->X;
	for (int b = 0; b < X; b1++) {
		wf->Ps[b*X+b] = 1.0;
	}
	wf->Z = 1;
	wf->N = 1;
	wg->G = G;
	return wf;
}

awf_t* setup(int Z, int N, int G, int maxN, int maxL, double* r, double** P0s) {
	awf_t* wf = construct_basis(Z, N, maxN, maxL, r);
	wf->Ps = P0s[l];
	wf->DM = (double*) malloc(X*X * sizeof(double));
	wf->Z = Z;
	wf->N = N;
	wf->G = G;
	return wf;
}

void solve(awf_t* wf, in nsteps, double tol) {
	double dE = 10000;
	int step = 0;

	while (step < nsteps && dE > tol) {
		make_density_matrix(wf);

		for (int mu = 0; mu < X, mu++) {
			for (int nu = 0; nu < X; nu++) {
				F[mu*X+nu] += h[mu*X+nu];
				for (int lambda = 0; lambda < X; lambda++) {
					for (int sigma = 0; sigma < X; sigma++) {
						F[mu*X+nu] += DM[lambda*X+sigma] * 
							(EE[mu*X+nu][lambda*X+sigma] - EE[mu*X+sigma][lambda*X+nu]);
					}
				}
			}
		}

		LAPACKE_dsyevd(LAPACK_ROW_MAJOR, 'V', 'U', X, F, X, wf.es);
		free(wf.Ps);
		wf.Ps = F;

		dE = wf->E;
		calc_energy(wf);

		dE = fabs(wf->e - dE);
		step++;
	}
}