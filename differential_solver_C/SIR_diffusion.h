#pragma once
typedef void (*rhs_grid_fn)(
    double* out,              
    const double* state,      
    size_t n_fields,
    size_t rows, 
    size_t cols,
    void* userdata
);

typedef void (*rhs_graph_fn)(
    double* out,             
    const double* state,      
    size_t n_fields,
    size_t vec_size,
    void* userdata
);

double* euler_solve_normal(const double *x0, size_t rows, size_t cols,
			   size_t n_fields,
                           double dx, double dy, double dt, size_t steps, 
			   size_t snapshots, rhs_grid_fn rhs, void*userdata);
double* euler_solve_graph(const double *L, const double *x0,
                          size_t vec_size, size_t n_fields, 
			  double dt, size_t steps, size_t snapshots, 
			  rhs_graph_fn rhs, void* userdata);
void free_array(double* p);

double* euler_solve_graph_csr(const int *indptr, const int *indices, const double *data,
                              const double *x_0, size_t vec_size, size_t n_fields,
                              double dt, size_t steps, size_t snapshots,
                              rhs_graph_fn rhs, void *userdata);

