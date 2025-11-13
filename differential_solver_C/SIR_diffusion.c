#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <string.h>

/*----------------------------HELPERS---------------------------------------*/

static inline size_t idx2D(size_t i, size_t j, size_t cols) { return i*cols + j; }
static inline size_t idx3D(size_t f, size_t i, size_t j, size_t rows, size_t cols) { return f*cols*rows+i*cols+j; } /* block offset for field f */

static void fill_targets(size_t steps, size_t snapshots, size_t *targets) {
    for (size_t k = 0; k < snapshots; ++k) {
        long long num = (long long)(k + 1) * (long long)steps;
        size_t t = (size_t)(num / (long long)snapshots);
        if (t < 1) t = 1;
        if (t > steps) t = steps;
        targets[k] = t;
    }
}

int mirror(int k, size_t N) {
    if (k < 0){        
	    return -k;
    }	    
    if ((size_t)k >= N){
	    return N-((size_t)k-N)-2;
    }
    return k;
}

void csr_matvec_buf(const int *indptr,
                    const int *indices,
                    const double *data,
                    const double *x,
                    size_t vec_size,
                    double *y){
    for (size_t i = 0; i < vec_size; ++i) {
        double s = 0.0;
        int row_start = indptr[i];
        int row_end   = indptr[i+1];
        for (int p = row_start; p < row_end; ++p){
            s += data[p] * x[(size_t)indices[p]];
        }
        y[i] = s;
    }
}

void matvec_buf(const double *A, 
		const double *x,             
		double *y, 
		size_t n){
    for (size_t i = 0; i < n; ++i) {
        double s = 0.0;
        for (size_t k = 0; k < n; ++k){
	       	s += A[i*n+k]*x[k];
	}
        y[i] = s;
    }
}


/*---------------------------RHS PREDECLARATION---------------------------*/

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

/*----------------------------START OF SOLVER-----------------------------*/

void finite_laplacian(const double *x, double *lap_buf,
                      size_t rows, size_t cols, double dx, double dy)
{
    const double inv_dx2 = 1.0/(dx*dx);
    const double inv_dy2 = 1.0/(dy*dy);
    for (int i = 0; i < (int)rows; ++i) {
        for (int j = 0; j < (int)cols; ++j) {
            size_t ip = (size_t)mirror(i+1, rows);
            size_t im = (size_t)mirror(i-1, rows);
            size_t jp = (size_t)mirror(j+1, cols);
            size_t jm = (size_t)mirror(j-1, cols);
            lap_buf[idx2D((size_t)i, (size_t)j, cols)] =
                ( x[idx2D(ip,(size_t)j,cols)] - 2.0*x[idx2D((size_t)i,(size_t)j,cols)] + x[idx2D(im,(size_t)j,cols)] ) * inv_dx2
              + ( x[idx2D((size_t)i,jp,cols)] - 2.0*x[idx2D((size_t)i,(size_t)j,cols)] + x[idx2D((size_t)i,jm,cols)] ) * inv_dy2;
        }
    }
}

void euler_step_graph(const double* L, 
		double* x, 
		double* lap_buf,
		double* rhs_buf,
		size_t vec_size,
		size_t n_fields,
		double dt,
		rhs_graph_fn rhs,
		void* userdata){
	rhs(rhs_buf, x, n_fields, vec_size, userdata);
	for (size_t f = 0; f<n_fields; f++){
		const double* x_f = x + f*vec_size; 
		const double* L_f = L + f*vec_size*vec_size;
		matvec_buf(L_f, x_f, lap_buf, vec_size);
		for (size_t i = 0; i < vec_size; ++i){
			x[idx2D(f, i, vec_size)] += -lap_buf[i]*dt + rhs_buf[idx2D(f, i, vec_size)]*dt;
		}
	}
}

void euler_step_graph_csr(const int *indptr,
                          const int *indices,
                          const double *data,
                          double *x,               
                          double *rhs_buf,         
                          double *lap_buf,           
                          size_t vec_size,
                          size_t n_fields,
                          double dt,
                          rhs_graph_fn rhs,
                          void *userdata)
{
    rhs(rhs_buf, x, n_fields, vec_size, userdata);

    for (size_t f = 0; f < n_fields; ++f) {
        const double *x_f = x + f*vec_size;
        csr_matvec_buf(indptr, indices, data, x_f, vec_size, lap_buf);       
        for (size_t i = 0; i < vec_size; ++i) {
            x[idx2D(f, i, vec_size)] += dt * ( -lap_buf[i] + rhs_buf[idx2D(f, i, vec_size)] );
        }
    }
}



void euler_step_normal(double* x,
		double* lap_buf,
		double* rhs_buf,
		size_t rows,
		size_t cols,
		size_t n_fields,
		double dx,
		double dy,
		double dt,
		rhs_grid_fn rhs,
		void* userdata){
	rhs(rhs_buf, x, n_fields, rows, cols, userdata);
	for (size_t f = 0; f < n_fields; f++){
		const double* x_f = x + f*rows*cols; 
		finite_laplacian(x_f, lap_buf, rows, cols, dx, dy);
		for (size_t i = 0; i < rows; i++){
			for(size_t j = 0; j < cols; ++j){
				x[idx3D(f, i, j, rows, cols)] += 
					lap_buf[idx2D(i,j,cols)]*dt + rhs_buf[idx3D(f, i, j, rows, cols)]*dt;
			}
		}
	}
}

double* euler_solve_graph_csr(const int *indptr, 
		              const int *indices, 
			      const double *data,
                              const double *x_0, 
			      size_t vec_size, 
			      size_t n_fields,
                              double dt, 
			      size_t steps, 
			      size_t snapshots,
                              rhs_graph_fn rhs, 
			      void *userdata){
    	const size_t N = vec_size*n_fields;
    	
	if (snapshots > steps){ 
		snapshots = steps;
	}

    	double *state = malloc(N*sizeof(double));
    	memcpy(state, x_0, N*sizeof(double));
    	
	double *snaps = malloc(snapshots*N*sizeof(double));
    	size_t *targets = malloc(snapshots*sizeof(size_t));
    	fill_targets(steps, snapshots, targets);

   	double *rhs_buf = malloc(N*sizeof(double));
	double *lap_buf = malloc(vec_size*sizeof(double));

    	size_t next_k = 0;
    	for (size_t step = 1; step <= steps; ++step){
        	euler_step_graph_csr(indptr, indices, data, state, rhs_buf,
				     lap_buf, vec_size, n_fields, dt, rhs, userdata);
        	if (next_k < snapshots && step == targets[next_k]){
            		memcpy(snaps + next_k*N, state, N*sizeof(double));
            		++next_k;
        	}
   	}
    	free(rhs_buf); 
	free(targets); 
	free(state);
    	return snaps;
}



double* euler_solve_normal(
        const double* x_0,
        size_t rows,
        size_t cols,
	size_t n_fields,
        double dx,
        double dy,
        double dt,
        size_t steps,
        size_t snapshots,
	rhs_grid_fn rhs,
       	void* userdata)   
{
    	const size_t N = rows * cols * n_fields;
	
	if (snapshots > steps){ 
	    snapshots = steps;
    	}

    	double* state = malloc(N * sizeof(double));
    	memcpy(state, x_0, N * sizeof(double));

    	double* snaps = malloc(snapshots * N * sizeof(double));
    	size_t* targets = malloc(snapshots * sizeof(size_t));
    
   	fill_targets(steps, snapshots, targets);
	double* lap_buf = malloc(rows*cols*sizeof(double));
	double* rhs_buf = malloc(N*sizeof(double));

    	size_t next_k = 0;
    	for (size_t step = 1; step <= steps; ++step) {
        	euler_step_normal(state, lap_buf, rhs_buf, rows, cols, n_fields, dx, dy, dt, rhs, userdata);
        	if (next_k < snapshots && step == targets[next_k]) {
            		memcpy(snaps + next_k * N, state, N * sizeof(double));
            		++next_k;
        	}
    	}

    	free(targets);
    	free(state);
	free(lap_buf);
	free(rhs_buf);
    	return snaps; 
}

double* euler_solve_graph(
        const double* L,
        const double* x_0,
        size_t vec_size,
	size_t n_fields,
        double dt,
        size_t steps,
        size_t snapshots,
	rhs_graph_fn rhs,
	void* userdata)
{
    	const size_t N = vec_size*n_fields;

    	if (snapshots > steps) snapshots = steps;

    	double* state = malloc(N * sizeof(double));
    	memcpy(state, x_0, N * sizeof(double));

    	double* snaps = malloc(snapshots * N * sizeof(double));
    	size_t* targets = malloc(snapshots * sizeof(size_t));
    
   	fill_targets(steps, snapshots, targets);

	double* lap_buf = malloc(vec_size * sizeof(double));
	double* rhs_buf = malloc(N * sizeof(double));

    	size_t next_k = 0;
    	for (size_t step = 1; step <= steps; ++step) {
        	euler_step_graph(L, state, lap_buf, rhs_buf, vec_size, n_fields, dt, rhs, userdata);
        	if (next_k < snapshots && step == targets[next_k]) {
            		memcpy(snaps + next_k * N, state, N * sizeof(double));
            		++next_k;
        	}
    	}

    	free(targets);
    	free(state);
	free(lap_buf);
	free(rhs_buf);
	return snaps; 
}


void free_array(double* p){
	free(p);
}
