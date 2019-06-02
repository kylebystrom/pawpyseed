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

float complex pseudo_momentum(int* GP, int* G_bounds, double* lattice,
	int* G1s, float complex* C1s, int num_waves1,
	int* G2s, float complex* C2s, int num_waves2, int* fftg) {
	// sum(u_1'*(k'+G+G') u_2(k+G)) = <u_1 |   >
	// NOTE: NEED TO CHECK THAT G2s+GP IS NOT TOO BIG OR TOO SMALL

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
		x[g1*fftg[1]*fftg[2] + g2*fftg[2] + g3] = C1s[w];
	}
	for (int w = 0; w < num_waves2; w++) {
		if ( G2s[3*w+0]+GP[0] >= G_bounds[0] && G2s[3*w+0]+GP[0] <= G_bounds[1] 
		  && G2s[3*w+1]+GP[1] >= G_bounds[2] && G2s[3*w+1]+GP[1] <= G_bounds[3]
		  && G2s[3*w+2]+GP[2] >= G_bounds[4] && G2s[3*w+2]+GP[2] <= G_bounds[5]) {
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

void make_rho(double* rho, int size, double* aewave1, double* pswave1, double* aewave2, double* pswave2) {
	for (int i = 0; i < size; i++) {
		rho[i] = aewave1[i] * aewave2[i] - pswave1[i] * pswave2[i];
	}
}

density_ft_t spher_transforms(int size, double* r, double* f, int l1, int m1, int l2, int m2, double encut) {

	density_ft_t density;
	density.l1 = l1;
	density.l2 = l2;
	density.m1 = m1;
	density.m2 = m2;
	density.transforms = (transform_spline_t*) malloc((l1+l2+1 - abs(l1-l2)) * sizeof(transform_spline_t));

	double* ks = (double*) calloc(size, sizeof(double));
	sbt_descriptor_t* d = spherical_bessel_transform_setup(encut, 0, l1+l2, size, r, ks);

	density.ks = ks;

	for (int L = abs(l1-l2); L <= l1+l2; L++) {
		density.transforms[L-abs(l1-l2)].transform = wave_spherical_bessel_transform(d, f, L);
		density.transforms[L-abs(l1-l2)].spline = spline_coeff(ks, density.transforms[L-abs(l1-l2)].transform, size);

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

double complex spher_momentum(density_ft_t densities, double* G, double* lattice) {
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

	for (int L = abs(l1-l2); L <= l1+l2; L++) {
		int size = densities.size;
		double* k = densities.ks;
		double* f = transforms[L-abs(l1-l2)].transform;
		double** s = transforms[L-abs(l1-l2)].spline;
		frac_to_cartesian(G, lattice);

		double complex sph_val;
		double magG = mag(G);
		if (magG == 0) {
			sph_val = Ylm(L, m2-m1, 0, 0);
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
		total += SBTFACS[lx][ly][(L-abs(l1-l2))/2][lx+mx][my] * sph_val
				 * 4 * PI * cpow(I, L) * wave_interpolate(mag(G), size, k, f, s);
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
					make_rho(rho, pp.wave_gridsize, func1.aewave,
							func1.pswave, func2.aewave, func2.pswave);
					elem.densities[i*pp.total_projs+j] = spher_transforms(pp.wave_gridsize, pp.wave_grid,
															rho, l1, m1, l2, m2, encut);
					j++;
				}
			}
			i++;
		}
	}

	return elem;
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
	G[0] = kpoint1->k[0] - kpoint2->k[0] + GP[0];
	G[1] = kpoint1->k[1] - kpoint2->k[1] + GP[1];
	G[2] = kpoint1->k[2] - kpoint2->k[2] + GP[2];

	double complex phase;
	for (int s = 0; s < wf->num_sites; s++) {
		density_ft_elem_t elem = elems[labels[s]];
		phase = cexp(2 * PI * I * dot(G, coords+s*3));
		for (int i = 0; i < elem.total_projs; i++) {
			for (int j = 0; j < elem.total_projs; j++) {
				total += spher_momentum(elem.densities[i*elem.total_projs+j], G, wf->lattice)
						* conj(band1->projections[s].overlaps[i]) * band2->projections[s].overlaps[j]
						* phase;
			}
		}
	}

	return total;

}

double complex* get_momentum_matrix(pswf_t* wf, int* labels, double* coords,
									int band1, int kpt1, int spin1,
									int band2, int kpt2, int spin2,
									density_ft_elem_t* transforms_list,
									double encut) {
	double nb1max, nb2max, nb3max;
	int npmax;
	setup(wf->nspin, wf->nwk, wf->nband,
		  &nb1max, &nb2max, &nb3max, &npmax, encut,
		  wf->lattice, wf->reclattice);
	int* igall = malloc(3*npmax*sizeof(int));
	if (igall == NULL) {
	    ALLOCATION_FAILED();
	}

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

	double complex* matrix = (double complex*) malloc(ncnt * sizeof(double complex));
	// NEED TO GIVE BACK INFO ABOUT ncnt
	for (int i = 0; i < ncnt; i++) {
		int* GP = igall + 3*i;
		matrix[i] = get_momentum_matrix_element(wf, labels, coords, band1, kpt1, spin1,
												band2, kpt2, spin2, GP, transforms_list);
	}

	return matrix;
}


int* get_momentum_grid(pswf_t* wf, double encut) {

	double nb1max, nb2max, nb3max;
	int npmax;
	setup(wf->nspin, wf->nwk, wf->nband,
		  &nb1max, &nb2max, &nb3max, &npmax, encut,
		  wf->lattice, wf->reclattice);
	int* igall = malloc(3*npmax*sizeof(int));
	if (igall == NULL) {
	    ALLOCATION_FAILED();
	}

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
