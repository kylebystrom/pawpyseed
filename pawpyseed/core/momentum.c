#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <mkl.h>
#include "utils.h"
#include "reader.h"
#include "sbt.h"
#include "momentum.h"
#include "gaunt.h"

#define c 0.262465831
#define PI 3.14159265358979323846

void free_transform_spline_list(transform_spline_t* transforms, int num_transforms) {
	for (int i = 0; i < num_transforms; i++) {
		free(transforms[i].transform);
		free(transforms[i].spline[0]);
		free(transforms[i].spline[1]);
		free(transforms[i].spline[2]);
		free(transforms[i].spline);
	}
	free(transforms);
}

void free_density_ft_list(density_ft_t* densities, int total_projs) {
	int total_densities = total_projs * total_projs;
	for (int i = 0; i < total_densities; i++) {
		int num_transforms = (densities[i].l1 + densities[i].l2
							- abs(densities[i].l1 - densities[i].l2)) / 2 + 1;
		free(densities[i].ks);
		free_transform_spline_list(densities[i].transforms, num_transforms);
	}
	free(densities);
}

void free_density_ft_elem_list(density_ft_elem_t* elems, int num_elems) {
	for (int i = 0; i < num_elems; i++) {
		free_density_ft_list(elems[i].densities, elems[i].total_projs);
	}
	free(elems);
}

float complex pseudo_momentum(int* GP, int* G_bounds, double* lattice,
	int* G1s, float complex* C1s, int num_waves1,
	int* G2s, float complex* C2s, int num_waves2, int* fftgrid) {
	// sum(u_1'*(k'+G+G') u_2(k+G)) = <u_1 |   >
	// NOTE: NEED TO CHECK THAT G2s+GP IS NOT TOO BIG OR TOO SMALL

	int fftg[3];
	fftg[0] = fftgrid[0]*2;
	fftg[1] = fftgrid[1]*2;
	fftg[2] = fftgrid[2]*2;
	float complex* x = (float complex*) mkl_calloc(fftg[0]*fftg[1]*fftg[2], sizeof(float complex), 64);

	float complex total = 0;
	int gridsize = fftg[0] * fftg[1] * fftg[2];
	for (int w = 0; w < gridsize; w++) {
		x[w] = 0;
	}
	int g1 = 0, g2 = 0, g3 = 0;
	for (int w = 0; w < num_waves1; w++) {
		g1 = (G1s[3*w+0]+fftg[0]) % fftg[0];
		g2 = (G1s[3*w+1]+fftg[1]) % fftg[1];
		g3 = (G1s[3*w+2]+fftg[2]) % fftg[2];
		x[g1*fftg[1]*fftg[2] + g2*fftg[2] + g3] = conj(C1s[w]);
	}
	for (int w = 0; w < num_waves2; w++) {
		if ( G2s[3*w+0]+GP[0] >= G_bounds[0] && G2s[3*w+0]+GP[0] <= G_bounds[1] 
		  && G2s[3*w+1]+GP[1] >= G_bounds[2] && G2s[3*w+1]+GP[1] <= G_bounds[3]
		  && G2s[3*w+2]+GP[2] >= G_bounds[4] && G2s[3*w+2]+GP[2] <= G_bounds[5]) {
			g1 = (G2s[3*w+0]+GP[0] +fftg[0]) % fftg[0];
			g2 = (G2s[3*w+1]+GP[1] +fftg[1]) % fftg[1];
			g3 = (G2s[3*w+2]+GP[2] +fftg[2]) % fftg[2];
			total += C2s[w] * x[g1*fftg[1]*fftg[2] + g2*fftg[2] + g3];
		}
		//else {
		//	x[g1*fftg[1]*fftg[2] + g2*fftg[2] + g3] = 0;
		//}
	}
	//for (int w = 0; w < gridsize; w++) {
	//	total += x[w];
	//}

	mkl_free(x);
	return total;
}

void mul_partial_waves(double* product, int size, double* r, double* f1, double* f2) {
	for (int i = 0; i < size; i++) {
		product[i] = f1[i] * f2[i] / r[i];
	}
}

void make_rho(double* rho, int size, double* grid, double* aewave1, double* pswave1, double* aewave2, double* pswave2) {

	for (int i = 0; i < size; i++) {
		rho[i] = (aewave1[i] * aewave2[i] - pswave1[i] * pswave2[i]) / grid[i];
	}
}

density_ft_t spher_transforms(int size, double* r, double* f, int l1, int m1, int l2, int m2, double encut) {

	density_ft_t density;
	density.l1 = l1;
	density.l2 = l2;
	density.m1 = m1;
	density.m2 = m2;
	density.size = size;
	density.transforms = (transform_spline_t*) malloc(((l1+l2 - abs(l1-l2))/2 + 1) * sizeof(transform_spline_t));

	double* ks = (double*) calloc(size, sizeof(double));
	sbt_descriptor_t* d = spherical_bessel_transform_setup(1e7, 0, l1+l2, size, r, ks);
	//printf("MINK %lf\n", ks[0]);

	density.ks = ks;

	for (int L = abs(l1-l2); L <= l1+l2; L+=2) {
		density.transforms[(L-abs(l1-l2))/2].transform = wave_spherical_bessel_transform(d, f, L);
		density.transforms[(L-abs(l1-l2))/2].spline = spline_coeff(ks,
			density.transforms[(L-abs(l1-l2))/2].transform, size);
		//total += SBTFACS[lx][ly][(L-abs(l1-l2))/2][lx+mx][my] * integral
	}

	free_sbt_descriptor(d);

	return density;

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

double complex spher_momentum(density_ft_t densities, double* G) {
	// NOTE: maybe change G->double so it can contain delta_k component?

	int l1=densities.l1, l2=densities.l2,
		m1=densities.m1, m2=densities.m2;
	transform_spline_t* transforms = densities.transforms;
	
	double complex total = 0;

	int lx, ly, mx, my;
	if (l1 < l2) {
		lx = l2;
		ly = l1;
		mx = m2;
		my = m1;
	}
	else {
		lx = l1;
		ly = l2;
		mx = m1;
		my = m2;
	}
	if (my < 0) {
		mx = -mx;
		my = -my;
	}

	for (int L = abs(l1-l2); L <= l1+l2; L+=2) {
		int size = densities.size;
		double* k = densities.ks;
		double* f = transforms[(L-abs(l1-l2))/2].transform;
		double** s = transforms[(L-abs(l1-l2))/2].spline;

		double complex sph_val = 0;
		double magG = mag(G);
		if (magG == 0) {
			if (L == 0 && m1 == m2)
				sph_val = Ylm(L, m2-m1, 0, 0);
			else
				sph_val = 0;
		}
		else {
			double theta = 0, phi = 0;
			theta = acos(G[2]/magG);
			if (magG - fabs(G[2]) == 0) phi = 0;
			else phi = acos(G[0] / pow(G[0]*G[0] + G[1]*G[1], 0.5));
			if (G[1] < 0) phi = 2*PI - phi;
			sph_val = Ylm(L, m2-m1, theta, phi);
		}

		//NOTE: issue with scaling r^-1
		//if (L == 0) {
		//	printf("%lf %lf\n", SBTFACS[lx][ly][(L-abs(l1-l2))/2][lx+mx][my]);
		//}
		total += SBTFACS[lx][ly][(L-abs(l1-l2))/2][lx+mx][my] * sph_val
				 * 4 * PI * cpow(I, L) * pow(-1, m2) * wave_interpolate(magG, size, k, f, s);

		//total += 8 * cpow(I, -L) * pow(-1, m1) * SBTFACS[lx][ly][(L-abs(l1-l2))/2][lx+mx][my]
		//			* sph_val * wave_interpolate(magG, size, k, f, s);

	}

	return total;
}

density_ft_elem_t get_transforms(ppot_t pp, double encut) {

	density_ft_elem_t elem;
	elem.total_projs = pp.total_projs;
	elem.num_densities = pp.total_projs * pp.total_projs;
	elem.densities = (density_ft_t*) malloc(elem.num_densities * sizeof(density_ft_t));

	int i = 0;
	for (int n1 = 0; n1 < pp.num_projs; n1++) {
		int l1 = pp.funcs[n1].l;
		funcset_t func1 = pp.funcs[n1];
		for (int m1 = -l1; m1 <= l1; m1++) {
			int j = 0;
			for (int n2 = 0; n2 < pp.num_projs; n2++) {
				int l2 = pp.funcs[n2].l;
				funcset_t func2 = pp.funcs[n2];
				for (int m2 = -l2; m2 <= l2; m2++) {

					double* rho = (double*) malloc(pp.wave_gridsize * sizeof(double));
					make_rho(rho, pp.wave_gridsize, pp.wave_grid, func1.aewave,
							func1.pswave, func2.aewave, func2.pswave);
					elem.densities[i*pp.total_projs+j] = spher_transforms(pp.wave_gridsize,
															pp.wave_grid,
															rho, l1, m1, l2, m2, encut);
					elem.densities[i*pp.total_projs+j].n1 = n1;
					elem.densities[i*pp.total_projs+j].n2 = n2;
					free(rho);
					j++;
				}
			}
			i++;
		}
	}

	return elem;
}

density_ft_elem_t* get_all_transforms(pswf_t* wf, double encut) {
	density_ft_elem_t* elems = (density_ft_elem_t*) malloc(wf->num_elems * sizeof(density_ft_elem_t));
	for (int i = 0; i < wf->num_elems; i++) {
		elems[i] = get_transforms(wf->pps[i], encut);
	}
	return elems;
}

double complex get_momentum_matrix_element(pswf_t* wf, int* labels, double* coords,
											int b1, int k1, int s1,
											int b2, int k2, int s2,
											int* GP, density_ft_elem_t* elems) {
	// Find < psi_{b1,k1,s1} | GP + k1 - k2 | psi_{b2,k2,s2} >
	kpoint_t* kpoint1 = wf->kpts[k1+s1*wf->nwk];
	kpoint_t* kpoint2 = wf->kpts[k2+s2*wf->nwk];
	band_t* band1 = kpoint1->bands[b1];
	band_t* band2 = kpoint2->bands[b2];

	double complex total = pseudo_momentum(GP, wf->G_bounds, wf->lattice, kpoint1->Gs, band1->Cs, kpoint1->num_waves,
											kpoint2->Gs, band2->Cs, kpoint2->num_waves, wf->fftg);

	double G[3];
	G[0] = GP[0];
	G[1] = GP[1];
	G[2] = GP[2];

	double Gcart[3];
	Gcart[0] = kpoint1->k[0] - kpoint2->k[0] + G[0];
	Gcart[1] = kpoint1->k[1] - kpoint2->k[1] + G[1];
	Gcart[2] = kpoint1->k[2] - kpoint2->k[2] + G[2];
	frac_to_cartesian(Gcart, wf->reclattice);

	//if ((GP[1] == 0) && (GP[2] == 0)) {
	//	printf("CHECK PSEUDO MOM %lf %lf %lf %lf %lf\n", creal(total), cimag(total),
	//		Gcart[0], Gcart[1], Gcart[2]);
	//}

	double complex phase;
	int ni, nj;
	for (int s = 0; s < wf->num_sites; s++) {
		density_ft_elem_t elem = elems[labels[s]];
		phase = cexp(2 * PI * I * dot(G, coords+s*3));
		for (int i = 0; i < elem.total_projs; i++) {
			for (int j = 0; j < elem.total_projs; j++) {
				//printf("spher_momentum %d %d %d %lf %lf %lf\n", s, i, j, Gcart[0], Gcart[1], Gcart[2]);
				if(0) {
				//if (mag(Gcart) <= 1e-10) {
					ppot_t pp = wf->pps[labels[s]];
					ni = elem.densities[i*elem.total_projs+j].n1;
					nj = elem.densities[i*elem.total_projs+j].n2;
					if (elem.densities[i*elem.total_projs+j].l1 == elem.densities[i*elem.total_projs+j].l2
						&& elem.densities[i*elem.total_projs+j].m1 == elem.densities[i*elem.total_projs+j].m2) {
						total += (pp.aepw_overlap_matrix[pp.num_projs*ni+nj]
							- pp.pspw_overlap_matrix[pp.num_projs*ni+nj])
							* conj(band1->projections[s].overlaps[i]) * band2->projections[s].overlaps[j];
					}
				}
				else {
					//if (elem.densities[i*elem.total_projs+j].l1 == elem.densities[i*elem.total_projs+j].l2
					//	&& elem.densities[i*elem.total_projs+j].m1 == elem.densities[i*elem.total_projs+j].m2){
					total += spher_momentum(elem.densities[i*elem.total_projs+j], Gcart)
						* conj(band1->projections[s].overlaps[i]) * band2->projections[s].overlaps[j]
						* phase;
				}
				//printf("done\n");
			}
		}
	}

	return total;

}

void get_momentum_matrix(double complex* matrix, int numg, int* igall,
						pswf_t* wf, int* labels, double* coords,
						int band1, int kpt1, int spin1,
						int band2, int kpt2, int spin2,
						density_ft_elem_t* transforms_list,
						double encut) {
	//double complex* matrix = (double complex*) malloc(ncnt * sizeof(double complex));
	// NEED TO GIVE BACK INFO ABOUT ncnt
	for (int i = 0; i < numg; i++) {
		//if (i%100==0)
		//	printf("ayo %d %d\n", i, numg);
		int* GP = igall + 3*i;
		matrix[i] = get_momentum_matrix_element(wf, labels, coords, band1, kpt1, spin1,
												band2, kpt2, spin2, GP, transforms_list);
	}
}

void momentum_grid_size(pswf_t* wf, double* nb1max, double* nb2max, double* nb3max,
						int* npmax, double encut) {
	setup(wf->nspin, wf->nwk, wf->nband,
		  nb1max, nb2max, nb3max, npmax, encut,
		  wf->lattice, wf->reclattice);
}

void grid_bounds(int* G_bounds, int* gdim, int* igall, int num_waves) {
	int* G = NULL;
	for (int w = 0; w < num_waves; w++) {
		G = igall + 3*w;
		if (G[0] < G_bounds[0]) G_bounds[0] = G[0];
		else if (G[0] > G_bounds[1]) G_bounds[1] = G[0];
		if (G[1] < G_bounds[2]) G_bounds[2] = G[1];
		else if (G[1] > G_bounds[3]) G_bounds[3] = G[1];
		if (G[2] < G_bounds[4]) G_bounds[4] = G[2];
		else if (G[2] > G_bounds[5]) G_bounds[5] = G[2];
	}
	printf("%d %d %d %d %d %d %d\n", G_bounds[0], G_bounds[1], G_bounds[2], G_bounds[3], G_bounds[4], G_bounds[5], num_waves);
	gdim[0] = G_bounds[1] - G_bounds[0] + 1;
	gdim[1] = G_bounds[3] - G_bounds[2] + 1;
	gdim[2] = G_bounds[5] - G_bounds[4] + 1;
	printf("%d %d %d\n", gdim[0], gdim[1], gdim[2]);
}

void list_to_grid_map(int* grid, int* G_bounds, int* gdim, int* igall, int num_waves) {
	int G[3];
	for (int w = 0; w < num_waves; w++) {
		G[0] = (igall[3*w+0]%gdim[0] + gdim[0]) % gdim[0];
		G[1] = (igall[3*w+1]%gdim[1] + gdim[1]) % gdim[1];
		G[2] = (igall[3*w+2]%gdim[2] + gdim[2]) % gdim[2];
		grid[G[0]*gdim[1]*gdim[2] + G[1]*gdim[2] + G[2]] = w;
	}
}

int get_momentum_grid(int* igall, pswf_t* wf, double nb1max, double nb2max, double nb3max, double encut) {

	double* b1 = wf->reclattice + 0;
	double* b2 = wf->reclattice + 3;
	double* b3 = wf->reclattice + 6;

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
					sumkg[j] = (ig1p) * b1[j]
								+ (ig2p) * b2[j]
								+ (ig3p) * b3[j];
				}
				double gtot = mag(sumkg);
				double etot = pow(gtot,2.0) / c;
				//printf("%lf %lf\n", etot, gtot);
				if (etot <= encut) {
					ncnt++;
					igall[ncnt*3+0] = ig1p;
					igall[ncnt*3+1] = ig2p;
					igall[ncnt*3+2] = ig3p;
					//if (ig1p < wf->G_bounds[0]) wf->G_bounds[0] = ig1p;
					//else if (ig1p > wf->G_bounds[1]) wf->G_bounds[1] = ig1p;
					//if (ig2p < wf->G_bounds[2]) wf->G_bounds[2] = ig2p;
					//else if (ig2p > wf->G_bounds[3]) wf->G_bounds[3] = ig2p;
					//if (ig3p < wf->G_bounds[4]) wf->G_bounds[4] = ig3p;
					//else if (ig3p > wf->G_bounds[5]) wf->G_bounds[5] = ig3p;
				}
			}
		}
	}
	ncnt++;

	return ncnt;
}

void fill_grid(float complex* x, int* Gs, float complex* Cs, int* fftg, int numg) {

	int gridsize = fftg[0] * fftg[1] * fftg[2];
	for (int w = 0; w < gridsize; w++) {
		x[w] = 0;
	}
	int g1 = 0, g2 = 0, g3 = 0;
	for (int w = 0; w < numg; w++) {
		g1 = (Gs[3*w+0]%fftg[0] + fftg[0]) % fftg[0];
		g2 = (Gs[3*w+1]%fftg[1] + fftg[1]) % fftg[1];
		g3 = (Gs[3*w+2]%fftg[2] + fftg[2]) % fftg[2];
		x[g1*fftg[1]*fftg[2] + g2*fftg[2] + g3] = Cs[w];
	}
	
}

void fullwf_reciprocal(double complex* Cs, int* igall, pswf_t* wf, int numg,
	int band_num, int kpt_num, int* labels, double* coords) {

	int G[3];
	int g1, g2, g3;
	double magG;
	double Gcart[3];
	int size;
	double* k;
	double* f;
	double** spline;
	ppot_t pp;
	int p;
	double complex sph_val;
	double complex phase;
	int l;
	double cart_coord[3];
	double inv_sqrt_vol = pow(determinant(wf->lattice), -0.5);

	int fftg[3];
	fftg[0] = wf->G_bounds[1] - wf->G_bounds[0] + 1;
	fftg[1] = wf->G_bounds[3] - wf->G_bounds[2] + 1;
	fftg[2] = wf->G_bounds[5] - wf->G_bounds[4] + 1;
	int gridsize = fftg[0] * fftg[1] * fftg[2];
	float complex* x = (float complex*) malloc(gridsize * sizeof(float complex));
	kpoint_t* kpt = wf->kpts[kpt_num];
	band_t* band = wf->kpts[kpt_num]->bands[band_num];
	fill_grid(x, kpt->Gs, band->Cs, fftg, kpt->num_waves);

	for (int w = 0; w < numg; w++) {
		G[0] = igall[3*w+0];
		G[1] = igall[3*w+1];
		G[2] = igall[3*w+2];
		g1 = (G[0]+fftg[0]) % fftg[0];
		g2 = (G[1]+fftg[1]) % fftg[1];
		g3 = (G[2]+fftg[2]) % fftg[2];
		if (   G[0] >= wf->G_bounds[0] && G[0] <= wf->G_bounds[1]
			&& G[1] >= wf->G_bounds[2] && G[1] <= wf->G_bounds[3]
			&& G[2] >= wf->G_bounds[4] && G[2] <= wf->G_bounds[5]) {
			Cs[w] += x[g1*fftg[1]*fftg[2] + g2*fftg[2] + g3];
		}

		Gcart[0] = -G[0] - kpt->k[0];
		Gcart[1] = -G[1] - kpt->k[1];
		Gcart[2] = -G[2] - kpt->k[2];
		frac_to_cartesian(Gcart, wf->reclattice);
		magG = mag(Gcart);

		for (int s = 0; s < wf->num_sites; s++) {
			pp = wf->pps[labels[s]];
			p = 0;
			size = pp.wave_gridsize;
			k = pp.kwave_grid;

			//cart_coord[0] = coords[3*s+0];
			//cart_coord[1] = coords[3*s+1];
			//cart_coord[2] = coords[3*s+2];
			//frac_to_cartesian(cart_coord, wf->lattice);
			phase = G[0] * coords[3*s+0] + G[1] * coords[3*s+1] + G[2] * coords[3*s+2];
			phase = cexp(-2 * PI * I * phase);

			for (int n = 0; n < band->projections[s].num_projs; n++) {
				l = pp.funcs[n].l;
				for (int m = -l; m <= l; m++) {
					f = pp.funcs[n].kwave;
					spline = pp.funcs[n].kwave_spline;

					Cs[w] += band->projections[s].overlaps[p]
							* kwave_value(k, f, spline, size, l, m, Gcart)
							* cpow(I, l) * phase
							* 4 * PI * inv_sqrt_vol;
					p++;
				}
			} 
		}
	}

	free(x);
}

double complex kwave_value(double* x, double* wave, double** spline, int size,
	int l, int m, double* pos) {
	double r = mag(pos);

	double radial_val = wave_interpolate(r, size, x, wave, spline);

	if (r==0) return Ylm(l, m, 0, 0) * radial_val;
	double theta = 0, phi = 0;
	theta = acos(pos[2]/r);
	if (r - fabs(pos[2]) == 0) phi = 0;
	else phi = acos(pos[0] / pow(pos[0]*pos[0] + pos[1]*pos[1], 0.5));
	if (pos[1] < 0) phi = 2*PI - phi;
	double complex sph_val = Ylm(l, m, theta, phi);
	return radial_val * sph_val;
}

double complex quick_overlap(int* dG, double complex* C1s, double complex* C2s, int numg,
	int* Gs, int* gmap, int* G_bounds, int* gdim) {

	int wp;
	int GP[3];
	double complex total = 0;
	for (int w = 0; w < numg; w++) {
		GP[0] = Gs[3*w+0] + dG[0];
		GP[1] = Gs[3*w+1] + dG[1];
		GP[2] = Gs[3*w+2] + dG[2];

		if (   GP[0] >= G_bounds[0] && GP[0] <= G_bounds[1]
			&& GP[1] >= G_bounds[2] && GP[1] <= G_bounds[3]
			&& GP[2] >= G_bounds[4] && GP[2] <= G_bounds[5]) {
			GP[0] = (GP[0]%gdim[0] + gdim[0]) % gdim[0];
			GP[1] = (GP[1]%gdim[1] + gdim[1]) % gdim[1];
			GP[2] = (GP[2]%gdim[2] + gdim[2]) % gdim[2];
			wp = gmap[GP[0]*gdim[1]*gdim[2] + GP[1]*gdim[2] + GP[2]];
			if (wp >= 0) {
				total += conj(C1s[wp]) * C2s[w];
			}
		}
	}

	return total;

}
