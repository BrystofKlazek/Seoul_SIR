#include <stdlib.h>
#include <stdio.h>

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
		int rows_a,
	       	int cols_a,
		int cols_b) {
	double *AB = calloc(rows_a*cols_b, sizeof(double));
	for (int i = 0; i < rows_a; i++) {
		for (int j = 0; j < cols_b; j ++){
			double sum = 0.0;
			for (int k = 0; k < cols_a; k++){
				sum += A[i*cols_a + k] * B[k*cols_b+j];	
			}
			AB[i*cols_b + j] = sum;
		}		
	}
	return AB;
}

double* finite_laplacian(double*x,
		int rows_x,
		int cols_x,
		double dx,
		double dy){
	double* lap = calloc(rows_x*cols_x, sizeof(double));
	for (int i = 0; i < rows_x; i++) {
		for (int j = 0; j < cols_x; j++){
			int i_plus = mirror(i+1, rows_x);
			int i_minus = mirror(i-1, rows_x);
			int j_plus = mirror(j+1, cols_x);
			int j_minus = mirror(j-1, cols_x);
			lap[i*cols_x + j] = (x[i_plus*cols_x+j]-2*x[i*cols_x+j]+x[i_minus*cols_x+j])/(dx*dx) +
				(x[i*cols_x+j_plus]-2*x[i*cols_x+j]+x[i*cols_x+j_minus])/(dy*dy);
		}
	}
	return lap;
}

void euler_step_graph(const double* L, 
		double* x, 
		int vec_size, 
		double dt){
	double *lap_M = matrix_matrix(L, x, vec_size, vec_size, 1);
	for (int i = 0; i < vec_size; i++){
		x[i] += -lap_M[i]*dt;
	}
	free(lap_M);
}

void euler_step_normal(double* x,
		int rows_x,
		int cols_x,
		double dx,
		double dy,
		double dt){
	double* lap = finite_laplacian(x, rows_x, cols_x, dx, dy);
	for (int i = 0; i < rows_x; i++){
		for(int j = 0; j < cols_x; j++){
			x[i*cols_x+j] += lap[i*cols_x+j]*dt;
		}
	}
	free(lap);
}

double* euler_solve_normal(const double* x_0,
		int rows_x,
		int cols_x,
		double dx,
		double dy,
		double dt,
		int steps,
		int snapshots){
	double* x = calloc(rows_x*cols_x, sizeof(double));
	for (int i = 0; i < rows_x; i++){
		for(int j = 0; j < cols_x; j++){
			x[i*cols_x+j]=x_0[i*cols_x+j];
		}
	}
	for (int step = 1; step <= steps; step++){
		euler_step_normal(x, rows_x, cols_x, dx, dy, dt);
	}
	return x;
}

double* euler_solve_graph(const double* L, 
		const double* x_0, 
		int vec_size,
		double dt,
		int steps,
		int snapshots){
	double* x = calloc(vec_size, sizeof(double));
	for (int i = 0; i < vec_size; i++){
		x[i] = x_0[i];
	}
	for (int step = 1; step <= steps; step++){
		euler_step_graph(L, x, vec_size, dt);
	}
	return x;
}

void free_array(double* p){
	free(p);
}
