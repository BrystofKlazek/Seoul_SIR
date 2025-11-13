# distutils: language = c
# cython: boundscheck=False, wraparound=False, cdivision=True
cimport numpy as cnp
import numpy as np
from libc.string cimport memcpy

ctypedef cnp.float64_t DTYPE_t

cnp.import_array()

# ---- C declarations (aliased, so no name clashes) ----
cdef extern from "SIR_diffusion.h":
    double* c_euler_solve_normal "euler_solve_normal"(const double* x0, size_t rows, int cols,
                                                      double dx, double dy, double dt,
                                                      size_t steps, int snapshots) nogil
    double* c_euler_solve_graph  "euler_solve_graph"(const double* L, const double* x0,
                                                     size_t n, double dt, int steps, int snapshots) nogil
    void    c_euler_step_normal  "euler_step_normal"(double* x, size_t rows, int cols,
                                                     double dx, double dy, double dt) nogil
    void    c_euler_step_graph   "euler_step_graph"(const double* L, double* x, size_t n, double dt) nogil
    void    c_free_array         "free_array"(double* p) nogil

# ---- Python-callable wrappers ----

def euler_solve_normal(cnp.ndarray x0, double dx, double dy, double dt, size_t steps, size_t snapshots):
    cdef size_t rows = <int>x0.shape[0]
    cdef size_t cols = <int>x0.shape[1]
    cdef size_t take = snapshots if snapshots <= steps else steps
    cdef size_t N = <size_t>(rows * cols)

    cdef const double* x0_ptr = <const double*> x0.data
    cdef double* out_ptr = c_euler_solve_normal(x0_ptr, rows, cols, dx, dy, dt, steps, snapshots)

    cdef cnp.ndarray out = np.empty((take, rows, cols), dtype=np.float64)
    memcpy(<void*>out.data, <const void*>out_ptr, <size_t>take * N * sizeof(DTYPE_t))
    c_free_array(out_ptr)
    return out

def euler_solve_graph(cnp.ndarray L, cnp.ndarray x0, double dt, size_t steps, size_t snapshots):
    cdef size_t n = <int>x0.size
    cdef size_t take = snapshots if snapshots <= steps else steps

    cdef const double* L_ptr  = <const double*> L.data
    cdef const double* x0_ptr = <const double*> x0.data
    cdef double* out_ptr = c_euler_solve_graph(L_ptr, x0_ptr, n, dt, steps, snapshots)

    cdef cnp.ndarray out = np.empty((take, n), dtype=np.float64)
    memcpy(<void*>out.data, <const void*>out_ptr, <size_t>take * <int>n * sizeof(DTYPE_t))
    c_free_array(out_ptr)
    return out

