

float complex pseudo_momentum(int* GP, int* G1_bounds, int* G2_bounds, double* lattice,
	int* G1s, float complex* C1s, int num_waves1,
	int* G2s, float complex* C2s, int num_waves2, int* fftg) {
	// sum(u_1'*(k'+G+G') u_2(k+G))
	// NOTE: NEED TO CHECK THAT G2s+GP IS NOT TOO BIG OR TOO SMALL

	float complex* x = mkl_calloc(fftg[0]*fftg[1]*fftg[2], sizeof(float complex), 64);

	float complex total = 0;
	int gridsize = fftg[0] * fftg[1] * fftg[2];
	for (int w = 0; w < gridsize; w++) {
		x[w] = 0;
	}
	int g1, g2, g3;
	for (int w = 0; w < num_waves; w++) {
		g1 = (G1s[3*w+0]+fftg[0]) % fftg[0];
		g2 = (G1s[3*w+1]+fftg[1]) % fftg[1];
		g3 = (G1s[3*w+2]+fftg[2]) % fftg[2];
		x[g1*fftg[1]*fftg[2] + g2*fftg[2] + g3] = C1s[w];
	}
	for (int w = 0; w < num_waves; w++) {
		if ( G2s[3*w+0]+GP[0] >= G2_bounds[0] && G2s[3*w+0]+GP[0] <= G2_bounds[1] 
		  && G2s[3*w+1]+GP[1] >= G2_bounds[2] && G2s[3*w+1]+GP[1] <= G2_bounds[3]
		  && G2s[3*w+2]+GP[2] >= G2_bounds[4] && G2s[3*w+2]+GP[2] <= G2_bounds[5]) {
			g1 = (G2s[3*w+0]+GP[0] +fftg[0]) % fftg[0];
			g2 = (G2s[3*w+1]+GP[1] +fftg[1]) % fftg[1];
			g3 = (G2s[3*w+2]+GP[2] +fftg[2]) % fftg[2];
			x[g1*fftg[1]*fftg[2] + g2*fftg[2] + g3] *= conj(C2s[w]);
		} else {
			x[g1*fftg[1]*fftg[2] + g2*fftg[2] + g3] = 0;
		}
	}
	for (int w = 0; w < fftg[0] * fftg[1] * fftg[2]; w++) {
		total += x[w];
	}
	return total;
}

void mul_partial_waves(double* product, int size, double* r, double* f1, double* f2) {
	for (int i = 0; i < size; i++) {
		product[i] = f1[i] * f2[i] / r[i];
	}
}

double** spher_transforms(int size, double* r, double* f1, double** s1, int l1, int m1,
												  double* f2, double** s2, int l2, int m2) {

	double* ks = (double*) calloc(size, sizeof(double));
	sbt_descriptor_t* d = spherical_bessel_transform_setup(encut, 0, l1+l2, size, r, ks);

	double* product = (double*) calloc(size, sizeof(double));
	mul_partial_waves(product, size, r, f1, f2);
	double** transforms = (double**) malloc(l1+l2+1-abs(l1-l2) * sizeof(double*));

	for (int L = abs(l1-l2); L <= l1+l2; L++) {
		transforms[L-abs(l1-l2)] = wave_spherical_bessel_transform(d, product, L);

		//total += SBTFACS[lx][ly][(L-abs(l1-l2))/2][lx+mx][my] * integral
	}

	return transforms;

	// < f1 | G | f2 >
	/*
	double k, kprime, val1, val2;
	for (int i = 0; i < size2; i++) {
		k = k2[i];
		kprime = k2[i] + k;
		val2 = f2[i];
		val1 = wave_interpolate(kprime, size1, k1, f1, s1);
		integral += val1 * val2;
	}
	*/
}

double complex spher_momentum(density_ft_t densities, int* GP) {
	// NOTE: maybe change GP->double so it can contain delta_k component?

	int l1=densities.l1, l2=densities.l2,
		m1=densities.m1, m2=densities.m2;
	transform_spline_t* transforms = densities.tranforms;
	
	double complex total = 0;

	for (int L = abs(l1-l2); L <= l1+l2; L++) {
		double* f = transforms[L-abs(l1-l2)].transform;
		double** s = transforms[L-abs(l1-l2)].spline;
		frac_to_cartesian(GP, lattice);

		double complex sph_val;
		if (magG == 0) {
			sph_val = Ylm(L, m2-m1, 0, 0);
		}
		else {
			double theta = 0, phi = 0;
			theta = acos(G[2]/magG);
			if (r - fabs(G[2]) == 0) phi = 0;
			else phi = acos(G[0] / pow(G[0]*G[0] + G[1]*G[1], 0.5));
			if (G[1] < 0) phi = 2*PI - phi;
			sph_val = Ylm(L, m2-m1, theta, phi);
		}

		//NOTE: issue with scaling r^-1
		total += SBTFACS[lx][ly][(L-abs(l1-l2))/2][lx+mx][my]
				 * 4 * PI * cpow(I, L) * wave_interpolate(mag(G), size, k, f, s);
	}

	return total;
}

double*** get_transforms(ppot_t pp) {

	double*** transforms_list = (double***) malloc(pp.total_projs * pp.total_projs * sizeof(double**));

	for (int n1 = 0; n1 < pp.num_projs; n1++) {
		int l1 = pp.funcs[n1].l;
		funcset_t func1 = pp.funcs[n1];
		for (int m1 = -l1; m1 <= l1; m1++) {
			for (int n2 = 0; n2 < pp.num_projs; n2++) {
				int l2 = pp.funcs[n2].l;
				funcset_t func2 = pp.funcs[n2];
				for (int m2 = -l2; m2 <= l2; m2++) {
					double** aetransforms = spher_transforms(pp.wave_gridsize, pp.wave_grid,
														func1.aewave, func1.aewave_spline, l1, m1,
														func2.aewave, func2.aewave_spline, l2, m2);
					double** pstransforms = spher_transforms(pp.wave_gridsize, pp.wave_grid,
														func1.pswave, func1.pswave_spline, l1, m1,
														func2.pswave, func2.pswave_spline, l2, m2);
					transforms_list[i*pp.total_projs+j] = aetransforms - pstransforms;
				}
			}
		}
	}

	return transforms_list;
}

double complex get_momentum_matrix_element(pswf_t* wf, int b1, int k1, int s1,
													   int b2, int k2, int s2,
													   int* GP, density_ft_elem_t* transforms_list) {
	band_t* band1 = wf->kpts[k1+s1*wf->nwk]->bands[b1];
	band_t* band2 = wf->kpts[k2+s2*wf->nwk]->bands[b2];
	pseudo_momentum(GP, G_bounds, wf->lattice, band1->Gs, band1->Cs, band1->num_waves,
											band2->Gs, band2->Cs, band2->num_waves, wf->fftg);

	for (int s = 0; s < num_sites; s++) {
		density_ft_elem_t transforms = transforms_list[labels[s]];
		phase = cexp(2 * PI * I * dot(GP, coords[s]));
		for (int i = 0; i < transforms.total_projs; i++) {
			for (int j = 0; j < transforms.total_projs; j++) {
				spher_momentum(transforms.densities[i*tranforms.total_projs+j], GP)
				* conj(band1->projections[s].overlaps[i]) * band2->projections[s].overlaps[j]
				* phase;
			}
		}

	}
}

double complex** get_momentum_matrix(pswf_t* wf, int b1, int k1, int s1,
												int b2, in k2, int s2,
												density_ft_elem_t* transforms_list,
												double encut) {
	get_momentum_grid(wf, encut);
	for (int i = 0; i < numg; i++) {
		int* GP = igall[3*i];
		get_momentum_matrix_element(wf, b1,k1,s1, b2,k2,s2, GP, transforms_list);
	}

	return matrix;
}


int* get_momentum_grid(pswf_t* wf, double encut) {

	double nb1max, nb2max, nb3max, npmax;
	setup(wf->nspin, wf->nwk, wf->nband,
		  &nb1max, &nb2max, &nb3max, &npmax, encut,
		  wf->lattice, wf->reclattice);
	int* igall = malloc(3*npmax*sizeof(int));
	if (igall == NULL) {
	    ALLOCATION_FAILED();
	}

	int ncnt = -1;
	for (int ig3 = 0; ig3 <= 2 * nb3max; ig3++) {
		int ig3p = ig3;
		if (ig3 > nb3max) ig3p = ig3 - 2 * nb3max - 1;
		for (int ig2 = 0; ig2 <= 2 * nb2max; ig2++) {
			int ig2p = ig2;
			if (ig2 > nb2max) ig2p = ig2 - 2 * nb2max - 1;
			for (int ig1 = 0; ig1 <= 2 * nb1max; ig1++) {
				int ig1p = ig1;
				if (ig1 > nb1max) ig1p = ig1 - 2 * nb1max - 1;
				double sumkg[3];
				for (int j = 0; j < 3; j++) {
					sumkg[j] = (kx+ig1p) * b1[j]
								+ (ky+ig2p) * b2[j]
								+ (kz+ig3p) * b3[j];
				}
				double gtot = mag(sumkg);
				double etot = pow(gtot,2.0) / c;
				//printf("%lf %lf\n", etot, gtot);
				if (etot <= encut) {
					ncnt++;
					igall[ncnt*3+0] = ig1p;
					igall[ncnt*3+1] = ig2p;
					igall[ncnt*3+2] = ig3p;
					if (ig1p < wf->G_bounds[0]) wf->G_bounds[0] = ig1p;
					else if (ig1p > wf->G_bounds[1]) wf->G_bounds[1] = ig1p;
					if (ig2p < wf->G_bounds[2]) wf->G_bounds[2] = ig2p;
					else if (ig2p > wf->G_bounds[3]) wf->G_bounds[3] = ig2p;
					if (ig3p < wf->G_bounds[4]) wf->G_bounds[4] = ig3p;
					else if (ig3p > wf->G_bounds[5]) wf->G_bounds[5] = ig3p;
				}
			}
		}
	}
	ncnt++;

	return igall;
}
