# distutils: language = c
# cython: boundscheck=False, wraparound=False, cdivision=True
cimport numpy as cnp
import numpy as np
from libc.string cimport memcpy

ctypedef cnp.float64_t DTYPE_t

# ---- C declarations (aliased, so no name clashes) ----
cdef extern from "SIR_diffusion.h":
    double* c_euler_solve_normal "euler_solve_normal"(const double* x0, int rows, int cols,
                                                  double dx, double dy, double dt, int steps, int snapshots) nogil
    double* c_euler_solve_graph  "euler_solve_graph"(const double* L, const double* x0,
                                                 int n, double dt, int steps, int snapshots) nogil
    void    c_euler_step_normal  "euler_step_normal"(double* x, int rows, int cols,
                                                     double dx, double dy, double dt) nogil
    void    c_euler_step_graph   "euler_step_graph"(const double* L, double* x, int n, double dt) nogil
    void    c_free_array         "free_array"(double* p) nogil

# ---- Python-callable wrappers ----

def euler_solve_normal(cnp.ndarray x0, double dx, double dy, double dt, int steps, int snapshots):
    cdef int rows = <int>x0.shape[0]
    cdef int cols = <int>x0.shape[1]
    cdef const double* x0_ptr = <const double*> x0.data
    cdef double* out_ptr = c_euler_solve_normal(x0_ptr, rows, cols, dx, dy, dt, steps, snapshots)

    cdef cnp.ndarray out = np.empty((rows, cols), dtype=np.float64)
    memcpy(<void*>out.data, <const void*>out_ptr, rows*cols*sizeof(DTYPE_t))
    c_free_array(out_ptr)
    return out

def euler_solve_graph(cnp.ndarray L, cnp.ndarray x0, double dt, int steps, int snapshots):
    cdef int n = <int>x0.size
    cdef const double* L_ptr  = <const double*> L.data
    cdef const double* x0_ptr = <const double*> x0.data
    cdef double* out_ptr = c_euler_solve_graph(L_ptr, x0_ptr, n, dt, steps, snapshots)

    cdef cnp.ndarray out = np.empty((n,), dtype=np.float64)
    memcpy(<void*>out.data, <const void*>out_ptr, n*sizeof(DTYPE_t))
    c_free_array(out_ptr)
    return out
