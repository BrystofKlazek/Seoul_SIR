#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <string.h>
#include <math.h>

/*----------------------------HELPERS---------------------------------------*/

static size_t idx2D(size_t i, size_t j, size_t cols) { 
	return i*cols + j; 
}

static size_t idx3D(size_t f, size_t i, size_t j, size_t rows, size_t cols) { 
	return f*cols*rows+i*cols+j; 
} 

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

void csr_matvec(const int *indptr,
                const int *indices,
                const double *data,   
                const double *x,
                size_t vec_size,
                double *y)
{
    for (size_t i = 0; i < vec_size; ++i) {
        double s = 0.0;
        int row_start = indptr[i];
        int row_end   = indptr[i+1];

        for (int p = row_start; p < row_end; ++p) {
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

static double interp_hourly_weight(
    const double *hourly_weights,
    double t_hours,
    size_t i,
    size_t j,
    size_t vec_size,
    size_t n_hours
) {
    double period = (double)n_hours;

    // bring t_hours into [0, period) using modulo
    double t = fmod(t_hours, period);
    if (t < 0.0) {
        t += period;  
    }

    // now t ∈ [0, period), interpret as between hour h0 and h1 = (h0+1) mod n_hours
    size_t h0    = (size_t)t;        // 0 <= h0 <= n_hours-1
    double alpha = t - (double)h0;   // fractional part in [0,1)
    size_t h1    = (h0 + 1) % n_hours;

    double w0 = hourly_weights[idx3D(h0, i, j, vec_size, vec_size)];
    double w1 = hourly_weights[idx3D(h1, i, j, vec_size, vec_size)];

    // linear interpolation on the circle of hours
    return (1.0 - alpha) * w0 + alpha * w1;
}



/*-----------------------------MATRIX BUILDING---------------------------*/
/*
static void build_laplacian_at_time(
    double *L,
    const double *hourly_weights,
    double t_hours,
    size_t vec_size,
    size_t n_fields,
    size_t n_hours
) {
    for (size_t f = 0; f < n_fields; ++f) {
        for (size_t i = 0; i < vec_size; ++i) {
            double deg = 0.0;

            for (size_t j = 0; j < vec_size; ++j) {
                double w = interp_hourly_weight(
                    hourly_weights,
                    t_hours,
                    i, j,
                    vec_size,
                    n_hours
                );

                if (i != j) {
                    L[idx3D(f, i, j, vec_size, vec_size)] = -w;
                }

                deg += w;
            }

            L[idx3D(f, i, i, vec_size, vec_size)] = deg;
        }
    }
}
*/
static void build_laplacian_at_time(
    double *L,
    const double *hourly_weights,
    double t_hours,
    size_t vec_size,
    size_t n_fields,
    size_t n_hours
) {
    for (size_t f = 0; f < n_fields; ++f) {
        for (size_t i = 0; i < vec_size; ++i) {
            double deg = 0.0;

            for (size_t j = 0; j < vec_size; ++j) {
                if (i == j) {
                    // skip self-flows: they do not contribute to diffusion
                    continue;
                }

                double w = interp_hourly_weight(
                    hourly_weights,
                    t_hours,
                    i, j,
                    vec_size,
                    n_hours
                );

                // off-diagonal i -> j
                L[idx3D(f, i, j, vec_size, vec_size)] = -w;
                deg += w;
            }

            // diagonal: sum of outgoing *to other nodes only*
            L[idx3D(f, i, i, vec_size, vec_size)] = deg;
        }
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


void euler_step_graph(
    const double *hourly_weights,
    size_t n_hours,
    double t_hours,
    double *L_buf,
    double *x,
    double *lap_buf,
    double *rhs_buf,
    size_t vec_size,
    size_t n_fields,
    double dt,
    rhs_graph_fn rhs,
    void *userdata
) {
    build_laplacian_at_time(
        L_buf,
        hourly_weights,
        t_hours,
        vec_size,
        n_fields,
        n_hours
    );

    rhs(rhs_buf, x, n_fields, vec_size, userdata);

    for (size_t f = 0; f < n_fields; ++f) {
        double       *x_f   = x       + f * vec_size;             
        const double *rhs_f = rhs_buf + f * vec_size;              
        const double *L_f   = L_buf   + f * vec_size * vec_size;    

        matvec_buf(L_f, x_f, lap_buf, vec_size);

        for (size_t i = 0; i < vec_size; ++i) {
            x_f[i] += dt * ( -lap_buf[i] + rhs_f[i] );
        }
    }
}


void euler_step_graph_csr(const int * const *indptr_arr,
                          const int * const *indices_arr,
                          const double * const *data_arr,
                          size_t step,              
                          double *x,               
                          double *rhs_buf,        
                          double *lap_buf,       
                          size_t vec_size,
                          size_t n_fields,
                          double dt,
                          rhs_graph_fn rhs,
                          void *userdata)
{
    const int    *indptr  = indptr_arr[step];
    const int    *indices = indices_arr[step];
    const double *data    = data_arr[step];

    rhs(rhs_buf, x, n_fields, vec_size, userdata);

    for (size_t f = 0; f < n_fields; ++f) {
        const double *x_f = x + f * vec_size;

        // L_step * x_f -> lap_buf
        csr_matvec(indptr, indices, data, x_f, vec_size, lap_buf);

        // explicit Euler update
        for (size_t i = 0; i < vec_size; ++i) {
            size_t idx = idx2D(f, i, vec_size);
            x[idx] += dt * ( -lap_buf[i] + rhs_buf[idx] );
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

double* euler_solve_graph(
    const double *hourly_weights,
    const double *x_0,
    size_t vec_size,
    size_t n_fields,
    size_t n_hours,
    double t0_hours,
    double dt,
    size_t steps,
    size_t snapshots,
    rhs_graph_fn rhs,
    void *userdata
) {
    const size_t N = vec_size * n_fields;

    if (snapshots > steps) {
        snapshots = steps;
    }

    double *state = malloc(N * sizeof(double));
    if (!state) return NULL;
    memcpy(state, x_0, N * sizeof(double));

    double *snaps = malloc(snapshots * N * sizeof(double));

    size_t *targets = malloc(snapshots * sizeof(size_t));
    fill_targets(steps, snapshots, targets);

    double *rhs_buf = malloc(N * sizeof(double));
    double *lap_buf = malloc(vec_size * sizeof(double));
    double *L_buf   = malloc(n_fields * vec_size * vec_size * sizeof(double));

    size_t next_k = 0;
    double t = t0_hours;

    for (size_t step_idx = 0; step_idx < steps; ++step_idx) {
        size_t step_num = step_idx + 1;  

        euler_step_graph(
            hourly_weights,
            n_hours,
            t,
            L_buf,
            state,
            lap_buf,
            rhs_buf,
            vec_size,
            n_fields,
            dt,
            rhs,
            userdata
        );

        t += dt;

        if (next_k < snapshots && step_num == targets[next_k]) {
            memcpy(snaps + next_k * N, state, N * sizeof(double));
            ++next_k;
        }
    }

    free(L_buf);
    free(rhs_buf);
    free(lap_buf);
    free(targets);
    free(state);
    return snaps;
}



double* euler_solve_graph_csr(const int **indptr,
                              const int **indices,
                              const double **data,
                              const double *x_0,
                              size_t vec_size,
                              size_t n_fields,
                              double dt,
                              size_t steps,
                              size_t snapshots,
                              rhs_graph_fn rhs,
                              void *userdata)
{
    const size_t N = vec_size * n_fields;

    if (snapshots > steps) {
        snapshots = steps;
    }

    double *state = malloc(N * sizeof(double));
    memcpy(state, x_0, N * sizeof(double));

    double *snaps = malloc(snapshots * N * sizeof(double));
    size_t *targets = malloc(snapshots * sizeof(size_t));
    fill_targets(steps, snapshots, targets);

    double *rhs_buf = malloc(N * sizeof(double));
    double *lap_buf = malloc(vec_size * sizeof(double));

    size_t next_k = 0;
    for (size_t step_idx = 0; step_idx < steps; ++step_idx) {
        size_t step_num = step_idx + 1;  // physical time step 1..steps

        euler_step_graph_csr(indptr, indices, data,
                             step_idx,       // index into CSR arrays
                             state, rhs_buf, lap_buf,
                             vec_size, n_fields,
                             dt, rhs, userdata);

        if (next_k < snapshots && step_num == targets[next_k]) {
            memcpy(snaps + next_k * N, state, N * sizeof(double));
            ++next_k;
        }
    }

    free(rhs_buf);
    free(lap_buf);
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

void free_array(double* p){
	free(p);
}
