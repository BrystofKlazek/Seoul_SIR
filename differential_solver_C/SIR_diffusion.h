#pragma once
double* euler_solve_normal(const double *x0, int rows, int cols,
                           double dx, double dy, double dt, int steps, int snapshots);
double* euler_solve_graph(const double *L, const double *x0,
                          int n, double dt, int steps, int snapshots);
void euler_step_normal(double *x, int rows, int cols,
                       double dx, double dy, double dt);

void euler_step_graph(const double *L, double *x, int n, double dt);

void free_array(double* p);
