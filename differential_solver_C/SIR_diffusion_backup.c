#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <string.h>

int mirror(int k, int N) {
    if (k < 0){        
	    return -k;
    }	    
    if (k >= N){
	    return N-((k-N))-2;
    }
    return k;
}


double* matrix_matrix(const double *A, 
		const double *B,  
		size_t rows_a,
	       	size_t cols_a,
		size_t cols_b) {
	double *AB = calloc(rows_a*cols_b, sizeof(double));
	for (size_t i = 0; i < rows_a; i++) {
		for (size_t j = 0; j < cols_b; j ++){
			double sum = 0.0;
			for (size_t k = 0; k < cols_a; k++){
				sum += A[i*cols_a + k] * B[k*cols_b+j];	
			}
			AB[i*cols_b + j] = sum;
		}		
	}
	return AB;
}

double* finite_laplacian(double*x,
		size_t rows_x,
		size_t cols_x,
		double dx,
		double dy){
	double* lap = calloc(rows_x*cols_x, sizeof(double));
	for (int i = 0; i < (int)rows_x; i++) {
		for (int j = 0; j < (int)cols_x; j++){
			size_t i_plus = (size_t)mirror(i+1, rows_x);
			size_t i_minus = (size_t)mirror(i-1, rows_x);
			size_t j_plus = (size_t)mirror(j+1, cols_x);
			size_t j_minus = (size_t)mirror(j-1, cols_x);
			lap[(size_t)i*cols_x + j] = (x[i_plus*cols_x+(size_t)j]-2*x[(size_t)i*cols_x+(size_t)j]+x[i_minus*cols_x+(size_t)j])/(dx*dx) +
				(x[(size_t)i*cols_x+j_plus]-2*x[(size_t)i*cols_x+(size_t)j]+x[(size_t)i*cols_x+j_minus])/(dy*dy);
		}
	}
	return lap;
}

void euler_step_graph(const double* L, 
		double* x, 
		size_t vec_size, 
		double dt){
	double *lap_M = matrix_matrix(L, x, vec_size, vec_size, 1);
	for (size_t i = 0; i < vec_size; i++){
		x[i] += -lap_M[i]*dt;
	}
	free(lap_M);
}

void euler_step_normal(double* x,
		size_t rows_x,
		size_t cols_x,
		double dx,
		double dy,
		double dt){
	double* lap = finite_laplacian(x, rows_x, cols_x, dx, dy);
	for (size_t i = 0; i < rows_x; i++){
		for(size_t j = 0; j < cols_x; j++){
			x[i*cols_x+j] += lap[i*cols_x+j]*dt;
		}
	}
	free(lap);
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

double* euler_solve_normal(
        const double* x_0,
        size_t rows_x,
        size_t cols_x,
        double dx,
        double dy,
        double dt,
        size_t steps,
        size_t snapshots           
	)   
{
    	const size_t N = rows_x * cols_x;
	
	if (snapshots > steps){ 
	    snapshots = steps;
    	}

    	double* state = malloc(N * sizeof(double));
    	memcpy(state, x_0, N * sizeof(double));

    	double* snaps = malloc(snapshots * N * sizeof(double));
    	size_t* targets = malloc(snapshots * sizeof(size_t));
    
   	fill_targets(steps, snapshots, targets);

    	size_t next_k = 0;
    	for (size_t step = 1; step <= steps; ++step) {
        	euler_step_normal(state, rows_x, cols_x, dx, dy, dt);
        	if (next_k < snapshots && step == targets[next_k]) {
            		memcpy(snaps + next_k * N, state, N * sizeof(double));
            		++next_k;
        	}
    	}

    	free(targets);
    	free(state);
    	return snaps; 
}

double* euler_solve_graph(
        const double* L,
        const double* x_0,
        size_t vec_size,
        double dt,
        size_t steps,
        size_t snapshots
        )
{
    const size_t N = vec_size;

    if (snapshots > steps) snapshots = steps;

    double* state = malloc(N * sizeof(double));
    memcpy(state, x_0, N * sizeof(double));

    double* snaps = malloc(snapshots * N * sizeof(double));
    size_t* targets = malloc(snapshots * sizeof(size_t));
    
    fill_targets(steps, snapshots, targets);

    size_t next_k = 0;
    for (size_t step = 1; step <= steps; ++step) {
        euler_step_graph(L, state, N, dt);
        if (next_k < snapshots && step == targets[next_k]) {
            memcpy(snaps + next_k * N, state, N * sizeof(double));
            ++next_k;
        }
    }

    free(targets);
    free(state);
    return snaps; 
}


void free_array(double* p){
	free(p);
}
