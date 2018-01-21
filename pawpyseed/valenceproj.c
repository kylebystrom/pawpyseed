#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <omp.h>
#include <time.h>

#define PI 3.14159265359
#define c 0.262465831

void vcross(double* res, double* top, double* bottom) {
	res[0] = top[1] * bottom[2] - top[2] * bottom[1];
	res[1] = top[2] * bottom[0] - top[0] * bottom[2];
	res[2] = top[0] * bottom[1] - top[1] * bottom[0];
}

double dot(double* x1, double* x2) {
	return x1[0] * x2[0] + x1[1] * x2[1] + x1[2] * x2[2];
}

double mag(double* x1) {
	return pow(dot(x1, x1), 0.5);
}

double determinant(double* m) {
	return m[0] * m[4] * m[8]
		+  m[1] * m[5] * m[6]
		+  m[2] * m[3] * m[7]
		-  m[2] * m[4] * m[6]
		-  m[1] * m[3] * m[8]
		-  m[0] * m[5] * m[7];
}

typedef struct band {
	int n;
	int num_waves;
	double occ;
	double N;
	double complex energy;
	int* Gs;
	float complex* Cs;
	double complex* C_grid;
} band_t;

typedef struct kpoint {
	short int up;
	int* Gs;
	double* k;
	double weight;
	int num_bands;
	band_t** bands;
} kpoint_t;

void free_kpoint(kpoint_t* kpt) {
	for (int b = 0; b < kpt->num_bands; b++) {
		band_t* curr_band = kpt->bands[b];
		free(curr_band->Cs);
		//free(curr_band->Gs);
		//free(curr_band->C_grid);
		free(curr_band);
	}
	//printf("ya");
	free(kpt->Gs);
	free(kpt->bands);
	//free(kpt->k);
	free(kpt);
}
/*
double kband_norm(kpoint_t* kpt, int band_index) {
	band_t* curr_band = kpt->bands[band_index];
	double total = 0;
	for (int i = 0; i < curr_band->num_waves; i++) {
		total += creal(conj(curr_band->Cs[i]) * curr_band->Cs[i]);
	}
	curr_band->N = pow(total, -0.5);
	return total;
}

*/
void ALLOCATION_FAILED() {
	printf("ALLOCATION FAILED");
	exit(-1);
}
/*
double complex kband_overlap(kpoint_t* kpt1, kpoint_t* kpt2, int band_index) {
	double k_overlap = 1;
	double dkx = kpt1->k[0] - kpt2->k[0];
	double dky = kpt1->k[1] - kpt2->k[1];
	double dkz = kpt1->k[2] - kpt2->k[2];
	double expx = cexp(I * 2 * PI * dkx) - 1;
	double expy = cexp(I * 2 * PI * dky) - 1;
	double expz = cexp(I * 2 * PI * dkz) - 1;
	double complex total = 0;
	band_t* band1 = kpt1->bands[band_index];
	band_t* band2 = kpt2->bands[band_index];
	int* G1s = band1->Gs;
	int* G2s = band2->Gs;
	double complex* C1s = band1->Cs;
	double complex* C2s = band2->Cs;
	for (int i = 0; i < band1->num_waves; i++) {
		for (int j = 0; j < band2->num_waves; j++) {
			double overlap = C1s[i] * conj(C2s[j]);
			if (abs(dkx) > 0.000001)
				overlap *= expx / I / (dkx + G1s[i * 3] - G2s[j * 3]) / 2 / PI;
			else if (abs(G1s[i * 3] - G2s[j * 3]) != 0)
				overlap = 0;
			if (abs(dky) > 0.000001)
				overlap *= expy / I / (dky + G1s[i * 3 + 1] - G2s[j * 3 + 1]) / 2 / PI;
			else if (G1s[i * 3 + 1] - G2s[j * 3 + 1] != 0)
				overlap = 0;
			if (abs(dkz) > 0.000001)
				overlap *= expz / I / (dkz + G1s[i * 3 + 2] - G2s[j * 3 + 2]) / 2 / PI;
			else if (G1s[i * 3 + 2] - G2s[j * 3 + 2] != 0)
				overlap = 0;
			total += overlap;
		}
	}
	return total;
}

double complex band_overlap(band_t* band1, band_t* band2, double* kpt1, double* kpt2) {
	double dkx = kpt1[0] - kpt2[0];
	double dky = kpt1[1] - kpt2[1];
	double dkz = kpt1[2] - kpt2[2];
	double complex expx = cexp(I * 2 * PI * dkx) - 1;
	double complex expy = cexp(I * 2 * PI * dky) - 1;
	double complex expz = cexp(I * 2 * PI * dkz) - 1;
	double complex total = 0;
	int* G1s = band1->Gs;
	int* G2s = band2->Gs;
	double complex* C1s = band1->Cs;
	double complex* C2s = band2->Cs;
	for (int i = 0; i < band1->num_waves; i++) {
			for (int j = 0; j < band2->num_waves; j++) {
				double complex overlap = C1s[i] * conj(C2s[j]);
				if (abs(dkx) > 0.000001)
					overlap *= expx / I / (dkx + G1s[i * 3] - G2s[j * 3]) / 2 / PI;
				else if (abs(G1s[i * 3] - G2s[j * 3]) != 0)
					overlap = 0;
				if (abs(dky) > 0.000001)
					overlap *= expy / I / (dky + G1s[i * 3 + 1] - G2s[j * 3 + 1]) / 2 / PI;
				else if (G1s[i * 3 + 1] - G2s[j * 3 + 1] != 0)
					overlap = 0;
				if (abs(dkz) > 0.000001)
					overlap *= expz / I / (dkz + G1s[i * 3 + 2] - G2s[j * 3 + 2]) / 2 / PI;
				else if (G1s[i * 3 + 2] - G2s[j * 3 + 2] != 0)
					overlap = 0;
				total += overlap;
			}
		}
	return total;
}

double normalization_factor(kpoint_t** kpts, int num_kpts, int band_index) {
	double total = 0;
	for (int i = 0; i < num_kpts; i++) {
		total += kband_norm(kpts[i], band_index);
		for (int j = i+1; j < num_kpts; j++) {
			total += 2 * creal(kband_overlap(kpts[i], kpts[j], band_index));
		}
	}
	return pow(total, -0.5);
}
*/
void setup(char* filename, int* pnrecl, int* pnspin, int* pnwk, int* pnband,
	double* nb1, double* nb2, double* nb3, double* ecut,
	double* lattice, double* reclattice) {

	double readin[3];
	FILE* fp = fopen(filename, "rb");
	double* ptr = readin;
	fread(ptr,24,1,fp);
	int nrecl = (int) round(readin[0]);
	int nspin = (int) round(readin[1]);
	int nprec = (int) round(readin[2]);
	fclose(fp);
	fp = fopen(filename, "rb");
	double readarr[nrecl/8];
	ptr = readarr;
	FILE* fp0 = fp;
	fseek(fp0,1*nrecl,0);
	fread(ptr,nrecl,1,fp0);
	int nwk = (int) round(readarr[0]);
	int nband = (int) round(readarr[1]);
	double encut = readarr[2];
	printf("encut %lf\n", encut);
	for (int i = 0; i < 9; i++) {
		lattice[i] = readarr[i+3];
	}
	vcross(reclattice+0, lattice+3, lattice+6);
	vcross(reclattice+3, lattice+6, lattice+0);
	vcross(reclattice+6, lattice+0, lattice+3);
	double* b1 = reclattice;
	double* b2 = reclattice+3;
	double* b3 = reclattice+6;
	double Vcell = determinant(lattice);
	for (int i = 0; i < 9; i++) {
		reclattice[i] *= 2.0 * PI / Vcell;
	}
	double magb1 = mag(reclattice+0);
	double magb2 = mag(reclattice+3);
	double magb3 = mag(reclattice+6);
	double vtemp[3];
	double vmag, sinphi123, phi123;
	
	double phi12 = acos(dot(reclattice+0, reclattice+3) / (magb1 * magb2));
	vcross(vtemp, reclattice+0, reclattice+3);
	vmag = mag(vtemp);
	sinphi123 = dot(reclattice+6, vtemp) / (vmag * magb3);
	double nb1maxA = pow(encut*c,0.5) / (magb1 * fabs(sin(phi12))) + 1;
	double nb2maxA = pow(encut*c,0.5) / (magb2 * fabs(sin(phi12))) + 1;
	double nb3maxA = pow(encut*c,0.5) / (magb3 * fabs(sinphi123)) + 1;
	printf("%lf %lf %lf %lf %lf %lf\n", sin(phi12), pow(encut*c,0.5), phi12, vmag, sinphi123, magb3*abs(sinphi123));
	int npmaxA = (int) round(4.0/3.0*PI*nb1maxA*nb2maxA*nb3maxA);
	
	double phi13 = acos(dot(reclattice+0, reclattice+6) / (magb1 * magb3));
	vcross(vtemp, reclattice+0, reclattice+6);
	vmag = mag(vtemp);
	sinphi123 = dot(reclattice+3, vtemp) / (vmag * magb2);
	double nb1maxB = pow(encut*c,0.5) / (magb1 * fabs(sin(phi13))) + 1;
	double nb2maxB = pow(encut*c,0.5) / (magb2 * fabs(sinphi123)) + 1;
	double nb3maxB = pow(encut*c,0.5) / (magb3 * fabs(sin(phi13))) + 1;
	int npmaxB = (int) round(4.0/3.0*PI*nb1maxB*nb2maxB*nb3maxB);

	double phi23 = acos(dot(reclattice+6, reclattice+3) / (magb3 * magb2));
	vcross(vtemp, reclattice+3, reclattice+6);
	vmag = mag(vtemp);
	sinphi123 = dot(reclattice+0, vtemp) / (vmag * magb1);
	double nb1maxC = pow(encut*c,0.5) / (magb1 * fabs(sinphi123)) + 1;
	double nb2maxC = pow(encut*c,0.5) / (magb2 * fabs(sin(phi23))) + 1;
	double nb3maxC = pow(encut*c,0.5) / (magb3 * fabs(sin(phi23))) + 1;
	int npmaxC = (int) round(4.0/3.0*PI*nb1maxC*nb2maxC*nb3maxC);

	double nb1max = fmax(nb1maxA, fmax(nb1maxB, nb1maxC));
	double nb2max = fmax(nb2maxA, fmax(nb2maxB, nb2maxC));
	double nb3max = fmax(nb3maxA, fmax(nb3maxB, nb3maxC));
	printf("%lf %lf %lf\n", nb1maxA, nb2maxA, nb3maxA);
	printf("%lf %lf %lf\n", nb1maxB, nb2maxB, nb3maxB);
	printf("%lf %lf %lf\n", nb1maxC, nb2maxC, nb3maxC);

	int npmax = npmaxA;
	if (npmaxB < npmax) npmax = npmaxB;
	if (npmaxC < npmax) npmax = npmaxC;

	printf("spin k band %d %d %d\n", nspin, nwk, nband);
	printf(" %lf %lf %lf\n", lattice[0], lattice[1], lattice[2]);
	printf(" %lf %lf %lf\n", lattice[3], lattice[4], lattice[5]);
	printf(" %lf %lf %lf\n", lattice[6], lattice[7], lattice[8]);
	printf(" %lf %lf %lf\n", reclattice[0], reclattice[1], reclattice[2]);
	printf(" %lf %lf %lf\n", reclattice[3], reclattice[4], reclattice[5]);
	printf(" %lf %lf %lf\n", reclattice[6], reclattice[7], reclattice[8]);
	printf("\n %lf %lf %lf %d %d %d\n", nb1max, nb2max, nb3max, npmaxA, npmaxB, npmaxC);

	*pnrecl = nrecl;
	*pnspin = nspin;
	*pnwk = nwk;
	*pnband = nband;
	*nb1 = nb1max;
	*nb2 = nb2max;
	*nb3 = nb3max;
	*ecut = encut;
}

kpoint_t** read_wavefunctions(int* G_bounds, double* kpt_weights, int* ns, int* nk, int* nb, char* filename) {

	clock_t start = clock();
		
	int nrecli, nspin, nwk, nband;
	double nb1max, nb2max, nb3max, encut;
	double lattice[9];
	double reclattice[9];
	setup(filename, &nrecli, &nspin, &nwk, &nband, &nb1max, &nb2max, &nb3max,
		&encut, lattice, reclattice);
	double* b1 = reclattice;
	double* b2 = reclattice+3;
	double* b3 = reclattice+6;

	long nrecl = (long)nrecli;

	double* occ = malloc(nband*sizeof(double));
	double* cener = malloc(nband*sizeof(double));
	if (occ == NULL || cener == NULL) {
		ALLOCATION_FAILED();
	}

	*ns = nspin;
	*nk = nwk * nspin;
	*nb = nband;

	kpoint_t** kpts = (kpoint_t**) malloc(nwk*nspin*sizeof(kpoint_t*));
	if (kpts == NULL) {
		ALLOCATION_FAILED();
	}

	double* kptr = (double*) malloc(nrecl);
	float complex* cptr = (float complex*) malloc(nrecl);
	
	//#pragma omp parallel for
	for (int iwk = 0; iwk < nwk * nspin; iwk++) {
		long irec = iwk * (long)(1 + nband);
		FILE* pfp = fopen(filename, "rb");
		kpoint_t* kpt = (kpoint_t*) malloc(sizeof(kpoint_t));
		kpt->num_bands = nband;
		band_t** bands = (band_t**) malloc(nband*sizeof(band_t*));
		kpt->bands = bands;
		if (kpt == NULL || bands == NULL) {
                	ALLOCATION_FAILED();
        	}
		//void* karr = malloc(nrecl);
		fseek(pfp, irec*nrecl+2*nrecl, SEEK_SET);
		//printf("important %d %d", (long)irec*nrecl+2*nrecl, SEEK_SET);
		fread(kptr, 8, nrecl/8, pfp);
		
		int nplane = (int) round(kptr[0]);
		int* igall = malloc(3*nplane*sizeof(int));
		if (igall == NULL) {
                	ALLOCATION_FAILED();
        	}
		double kx = kptr[1];
		double ky = kptr[2];
		double kz = kptr[3];
		for (int i = 0; i < nband; i++) {
			cener[i] = kptr[4+i*3];
			occ[i] = kptr[6+i*3];
			band_t* band = (band_t*) malloc(sizeof(band_t));
			band->num_waves = nplane;
			band->energy = kptr[4+i*3];
			band->occ = kptr[6+i*3];
			kpt->bands[i] = band;
		}

		//printf("\n");
		//printf("nplane %d %d\n", nplane, npmax);

		//printf("kpt %lf %lf %lf\n", kx, ky, kz);

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
					}
				}
			}
		}
		if (ncnt != nplane - 1) printf("ERROR %lf %lf %lf %lf", kx,ky,kz, c);
		//if (ncnt > npmax) printf("BIG ERROR");
		//printf("%d %d\n", ncnt, npmax);

		for (int iband = 0; iband < nband; iband++) {
			irec++;
			fseek(pfp, (long)irec*nrecl+2*(long)nrecl, SEEK_SET);
			fread(cptr, 8, nrecl/8, pfp);
			float complex* coeff = malloc(nplane*sizeof(float complex));
			for (int iplane = 0; iplane < nplane; iplane++) {
				coeff[iplane] = cptr[iplane];
			}
			kpt->bands[iband]->Cs = coeff;
			kpt->bands[iband]->Gs = igall;
		}
		
		//printf("iwk %d\n", iwk);
		kpt->weight = kpt_weights[iwk%nwk];
		kpt->Gs = igall;
		kpts[iwk] = kpt;
	}

	free(kptr);
	free(cptr);
	clock_t end = clock();
	printf("%lf seconds taken to read bulk\n", (double)(end-start) / CLOCKS_PER_SEC);

	return kpts;
}

kpoint_t** read_one_band(int* G_bounds, double* kpt_weights, int* ns, int* nk, int* nb, int BAND_NUM, char* filename) {
	clock_t start = clock();

	int nrecli=0, nspin=0, nwk=0, nband=0;
	double nb1max=0, nb2max=0, nb3max=0, encut=0;
	double lattice[9] = {0};
	double reclattice[9] = {0};
	setup(filename, &nrecli, &nspin, &nwk, &nband, &nb1max, &nb2max, &nb3max,
		&encut, lattice, reclattice);
	double* b1 = reclattice;
	double* b2 = reclattice+3;
	double* b3 = reclattice+6;

	long nrecl = (long) nrecli;

	double* kptr = (double*) malloc(nrecl);
	float complex* cptr = (float complex*) malloc(nrecl);

	kpoint_t** kpts = (kpoint_t**) malloc(nwk*nspin*sizeof(kpoint_t*));

	//#pragma omp parallel for
	for (int iwk = 0; iwk < nwk * nspin; iwk++) {
		int irec = iwk * (1 + nband);
		FILE* pfp = fopen(filename, "rb");
		kpoint_t* kpt = (kpoint_t*) malloc(sizeof(kpoint_t));
		kpt->num_bands = 1;
		band_t** bands = (band_t**) malloc(sizeof(band_t*));
		kpt->bands = bands;
		//void* karr = malloc(nrecl);
		fseek(pfp, (long)irec*nrecl+(long)2*nrecl, SEEK_SET);
		//printf("reading dat %d\n", fread(kptr, 8, nrecl/8, pfp));
		fread(kptr, 8, nrecl/8, pfp);

		int nplane = (int) round(kptr[0]);
		int* igall = malloc(3*nplane*sizeof(int));
		double kx = kptr[1];
		double ky = kptr[2];
		double kz = kptr[3];

		band_t* band = (band_t*) malloc(sizeof(band_t));
		band->num_waves = nplane;
		band->energy = kptr[4+BAND_NUM*3];
		band->occ = kptr[6+BAND_NUM*3];
		kpt->bands[0] = band;

		//printf("nplane %d %d\n", nplane, npmax);
		/*#pragma omp critical
		{
		printf("hm %lf %lf %lf\n", kptr[6+BAND_NUM*3], kptr[7+BAND_NUM*3], band->occ);
		printf("eek,iwk nrecl np e o kpt %d %d %d %lf %lf %lf %lf\n", iwk, nrecl, nplane, kx,ky,kz, band->energy);
		
		printf("got some max %lf %lf %lf\n\n", nb1max, nb2max, nb3max);
		}
		continue;*/
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
					}
				}
			}
		}
		if (ncnt != nplane - 1) printf("ERROR %d %d", ncnt, nplane);

		int iband = BAND_NUM;
		irec += 1 + BAND_NUM;
		fseek(pfp,(long) irec*nrecl+2*(long)nrecl, SEEK_SET);
		fread(cptr, 8, nrecl/8, pfp);
		float complex* coeff = malloc(nplane*sizeof(float complex));
		for (int iplane = 0; iplane < nplane; iplane++) {
			coeff[iplane] = cptr[iplane];
		}
		kpt->bands[0]->Cs = coeff;
		kpt->bands[0]->Gs = igall;
		kpt->weight = kpt_weights[iwk%nwk];
		kpt->Gs = igall;
		kpts[iwk] = kpt;
	}
	clock_t end = clock();
	free(kptr);
	free(cptr);
	printf("%lf seconds to read defect\n", (double)(end-start) / CLOCKS_PER_SEC);

	return kpts;
}

void get_band_projection(int BAND_NUM, int NUM_KPTS, int NUM_BANDS, kpoint_t** kpts, kpoint_t** kptspro, int* G_bounds, double* results) {

	clock_t start = clock();

	double* cband = (double*) calloc(NUM_KPTS, sizeof(double));
	double* vband = (double*) calloc(NUM_KPTS, sizeof(double));

	/*int MIN_GX = G_bounds[0];
	int MAX_GX = G_bounds[1];
	int MIN_GY = G_bounds[2];
	int MAX_GY = G_bounds[3];
	int MIN_GZ = G_bounds[4];
	int MAX_GZ = G_bounds[5];

	printf("%d %d %d %d %d %d\n", MIN_GX, MAX_GX, MIN_GY, MAX_GY, MIN_GZ, MAX_GZ);
	int RANGE_DGX = MAX_GX - MIN_GX;
	int RANGE_DGY = MAX_GY - MIN_GY;
	int RANGE_DGZ = MAX_GZ - MIN_GZ;

	int NUM_DGX = RANGE_DGX * 2 + 1;
	int NUM_DGY = RANGE_DGY * 2 + 1;
	int NUM_DGZ = RANGE_DGZ * 2 + 1;

	int RANGE_GX = RANGE_DGX + 1;
	int RANGE_GY = RANGE_DGY + 1;
	int RANGE_GZ = RANGE_DGZ + 1;*/

	printf("stuff %e %e\n", crealf(kpts[0]->bands[70]->Cs[19]), crealf(kptspro[0]->bands[0]->Cs[19]));

	#pragma omp parallel for 
	for (int b = 0; b < NUM_BANDS; b++)
	{
		//printf("occ %d %lf\n", b, kpts[0]->bands[b]->occ);
		for (int kpt_num = 0; kpt_num < NUM_KPTS; kpt_num++)
		{
			float complex curr_overlap = 0;
			float complex* C1s = kptspro[kpt_num]->bands[0]->Cs;
			float complex* C2s = kpts[kpt_num]->bands[b]->Cs;
			int num_waves = kpts[kpt_num]->bands[b]->num_waves;
			for (int w = 0; w < num_waves; w++)
			{
				curr_overlap += C1s[w] * conj(C2s[w]);
			}
			#pragma omp critical
			{
				if (kpts[kpt_num]->bands[b]->occ > 0.5)
					vband[kpt_num] += creal((double) (curr_overlap * conj(curr_overlap)));
				else
					cband[kpt_num] += creal((double) (curr_overlap * conj(curr_overlap)));
			}
		}
	}

	/*
	#pragma omp parallel for
	for (int b = 0; b < NUM_BANDS; b++) {
		//if (b == BAND_NUM) continue;
		for (int kpt_num = 0; kpt_num < NUM_KPTS; kpt_num++) {
			double complex curr_overlap = 0;
			double complex* grid1 = kptspro[kpt_num]->bands[0]->C_grid;
			double complex* grid2 = kpts[kpt_num]->bands[b]->C_grid;
			for (int Gz = MIN_GZ; Gz < MAX_GZ + 1; Gz++) {
				for (int Gy = MIN_GY; Gy < MAX_GY + 1; Gy++) {
					for (int Gx = MIN_GX; Gx < MAX_GX + 1; Gx++) {
						int index = (Gz - MIN_GZ) * RANGE_GY * RANGE_GX +
									(Gy - MIN_GY) * RANGE_GX +
									(Gx - MIN_GX);
						curr_overlap += grid1[index] * conj(grid2[index]);
					}
				}
			}
			#pragma omp critical
			{
				//printf("kkbbri %d %d %d %lf %lf\n", kpt_num, BAND_NUM, b, creal(curr_overlap), cimag(curr_overlap));
				if (kpts[kpt_num]->bands[b]->occ > 0.5) vband[kpt_num] += creal(curr_overlap * conj(curr_overlap));
				else cband[kpt_num] += creal(curr_overlap * conj(curr_overlap));
			}
		}
	}
	*/

	double ctotal = 0.0;
	double vtotal = 0.0;
	for (int kpt_num = 0; kpt_num < NUM_KPTS; kpt_num++) {
		//printf("%lf %lf\n", cband[kpt_num], kpts[kpt_num]->weight);
		ctotal += cband[kpt_num] * kpts[kpt_num]->weight;
		vtotal += vband[kpt_num] * kpts[kpt_num]->weight;
	}

	printf("%lf\n", creal(kptspro[0]->bands[0]->energy));
	printf("c %lf\n", ctotal);
	printf("v %lf\n", vtotal);

	free(vband);
	free(cband);
	results[0] = vtotal;
	results[1] = ctotal;

	clock_t end = clock();
	printf("%lf seconds for band projection\n", (double)(end - start) / CLOCKS_PER_SEC);

}

double* read_and_project(int BAND_NUM, double* kpt_weights, char* bulkfile, char* defectfile) {
	printf("%lf\n", kpt_weights[0]);
	printf("%lf\n", kpt_weights[1]);
	printf("%lf\n", kpt_weights[5]);
	int* G_bounds = (int*) malloc(6*sizeof(double));
	double* results = (double*) malloc(2*sizeof(double));
	int NUM_SPINS, NUM_KPTS, NUM_BANDS;
	kpoint_t** kptspro = read_one_band(G_bounds, kpt_weights, &NUM_SPINS, &NUM_KPTS, &NUM_BANDS, BAND_NUM, defectfile);
	kpoint_t** kptsref = read_wavefunctions(G_bounds, kpt_weights, &NUM_SPINS, &NUM_KPTS, &NUM_BANDS, bulkfile);
	get_band_projection(BAND_NUM, NUM_KPTS, NUM_BANDS, kptsref, kptspro, G_bounds, results);
	for (int kpt_num = 0; kpt_num < NUM_KPTS; kpt_num++) {
		//printf("%d\n", kpt_num);
		free_kpoint(kptsref[kpt_num]);
		free_kpoint(kptspro[kpt_num]);
	}
	free(kptsref);
	free(kptspro);
	free(G_bounds);
	return results;
}
