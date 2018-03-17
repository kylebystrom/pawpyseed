

double kinetic(double* R, double** states, int* ls, int size; int Z) {
	double kin = 0;
	for (int n = 0; n < Z; n++) {
		double* C = states[n];
		int l = ls[n];
		for (int i = 1; i < size-1; i++) {
			d1 = (C[i]-C[i-1]) / (R[i]-R[i-1]);
			d2 = (C[i+1]-C[i]) / (R[i+1]-R[i]);
			diff = (R[i+1]-R[i-1])/2;
			psipp = 1/R[i] * (d1+d2) + (d2-d1)/diff;
			kin += (R[i] * R[i] * psipp + l*(l+1)/2) * C[i] * C[i] * diff;
		}
	}
	return kin;
}

double elec_nuc(double* R, double** states, int* ls, int size, int Z) {
	double en = 0;
	for (int n = 0; n < Z; n++) {
		
	}
}

void atom_gga(int Z, )