
typedef struct transform_spline {
	double* transform;
	double** spline;
} transform_spline_t;

typedef struct density_ft {
	int n1;
	int l1;
	int m1;
	int n2;
	int l2;
	int m2;
	int size;
	transform_spline_t* transforms; // size: sum(l1,l2)-diff(l1,l2)
} density_ft_t;

typedef struct density_ft_elem {
	int num_densities;
	int total_projs;
	density_ft_t* densities;
} density_ft_elem_t;