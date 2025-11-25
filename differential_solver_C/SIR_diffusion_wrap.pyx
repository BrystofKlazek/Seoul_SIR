# distutils: language = c
# cython: boundscheck=False, wraparound=False, cdivision=True
cimport numpy as cnp
import numpy as np
from libc.string cimport memcpy
from libc.stddef cimport size_t
from libc.stdlib cimport malloc, free

ctypedef cnp.float64_t DTYPE_t
cnp.import_array()

cdef extern from "SIR_diffusion.h":
    ctypedef void (*rhs_grid_fn)(double* out, const double* state,
                                 size_t n_fields, size_t rows, size_t cols,
                                 void* userdata)
    ctypedef void (*rhs_graph_fn)(double* out, const double* state,
                                  size_t n_fields, size_t vec_size,
                                  void* userdata)

    double* euler_solve_normal(const double *x0, size_t rows, size_t cols,
                               size_t n_fields,
                               double dx, double dy, double dt,
                               size_t steps, size_t snapshots,
                               rhs_grid_fn rhs, void* userdata)

    double* euler_solve_graph(const double *hourly_weights,
                              const double *x0,
                              size_t vec_size, size_t n_fields,
                              size_t n_hours,
                              double t0_hours,
                              double dt,
                              size_t steps, size_t snapshots,
                              rhs_graph_fn rhs, void* userdata)

    double* euler_solve_graph_csr(const int **indptr_arr,
                                  const int **indices_arr,
                                  const double **data_arr,
                                  const double *x0,
                                  size_t vec_size, size_t n_fields,
                                  double dt, size_t steps, size_t snapshots,
                                  rhs_graph_fn rhs, void* userdata)

    void free_array(double* p)

cdef struct PyRHSContext:
    void* fn_obj
    void* params_obj


cdef void _rhs_grid_trampoline_gil(double* out, const double* state,
                                   size_t n_fields, size_t rows, size_t cols,
                                   void* userdata) with gil:
    cdef PyRHSContext* ctx = <PyRHSContext*>userdata
    cdef cnp.npy_intp d3[3]
    d3[0] = <cnp.npy_intp>n_fields
    d3[1] = <cnp.npy_intp>rows
    d3[2] = <cnp.npy_intp>cols
    py_fn = <object>ctx.fn_obj
    py_params = <object>ctx.params_obj
    state_arr = cnp.PyArray_SimpleNewFromData(3, &d3[0], cnp.NPY_FLOAT64, <void*>state)
    out_arr   = cnp.PyArray_SimpleNewFromData(3, &d3[0], cnp.NPY_FLOAT64, <void*>out)
    if py_params is None:
        py_fn(state_arr, out_arr)
    elif isinstance(py_params, dict):
        py_fn(state_arr, out_arr, **py_params)
    else:
        py_fn(state_arr, out_arr, params=py_params)

cdef void _rhs_grid_trampoline(double* out, const double* state,
                               size_t n_fields, size_t rows, size_t cols,
                               void* userdata) nogil:
    with gil:
        _rhs_grid_trampoline_gil(out, state, n_fields, rows, cols, userdata)


cdef void _rhs_graph_trampoline_gil(double* out, const double* state,
                                    size_t n_fields, size_t vec_size,
                                    void* userdata) with gil:
    cdef PyRHSContext* ctx = <PyRHSContext*>userdata
    cdef cnp.npy_intp d2_state[2]
    cdef cnp.npy_intp d2_out[2]
    d2_state[0] = <cnp.npy_intp>n_fields
    d2_state[1] = <cnp.npy_intp>vec_size
    d2_out[0]   = <cnp.npy_intp>n_fields
    d2_out[1]   = <cnp.npy_intp>vec_size
    py_fn = <object>ctx.fn_obj
    py_params = <object>ctx.params_obj
    state_arr = cnp.PyArray_SimpleNewFromData(2, &d2_state[0], cnp.NPY_FLOAT64, <void*>state)
    out_arr   = cnp.PyArray_SimpleNewFromData(2, &d2_out[0],   cnp.NPY_FLOAT64, <void*>out)
    if py_params is None:
        py_fn(state_arr, out_arr)
    elif isinstance(py_params, dict):
        py_fn(state_arr, out_arr, **py_params)
    else:
        py_fn(state_arr, out_arr, params=py_params)

cdef void _rhs_graph_trampoline(double* out, const double* state,
                                size_t n_fields, size_t vec_size,
                                void* userdata) nogil:
    with gil:
        _rhs_graph_trampoline_gil(out, state, n_fields, vec_size, userdata)


def euler_solve_normal_rd(cnp.ndarray x0,
                          double dx, double dy, double dt,
                          size_t steps, size_t snapshots,
                          rhs_fn, params=None):
    cdef size_t n_fields = <size_t>x0.shape[0]
    cdef size_t rows     = <size_t>x0.shape[1]
    cdef size_t cols     = <size_t>x0.shape[2]
    cdef size_t take     = snapshots if snapshots <= steps else steps
    cdef size_t Nfld     = rows*cols

    cdef PyRHSContext ctx
    ctx.fn_obj = <void*>rhs_fn
    ctx.params_obj = <void*>params

    cdef const double* x0_ptr = <const double*> x0.data
    cdef double* out_ptr = euler_solve_normal(x0_ptr, rows, cols, n_fields,
                                              dx, dy, dt, steps, snapshots,
                                              <rhs_grid_fn>_rhs_grid_trampoline,
                                              <void*>&ctx)

    cdef cnp.ndarray out = np.empty((take, n_fields, rows, cols), dtype=np.float64)
    memcpy(<void*>out.data, <const void*>out_ptr,
           <size_t>take * n_fields * Nfld * sizeof(DTYPE_t))
    free_array(out_ptr)
    return out


def euler_solve_graph_rd(cnp.ndarray hourly_weights,
                         cnp.ndarray x0,
                         double t0_hours,
                         double dt, size_t steps, size_t snapshots,
                         rhs_fn, params=None):
    # ensure dtypes / contiguity (turn to numpy, so it can be worked into the C code
    if hourly_weights.dtype != np.float64:
        hourly_weights = np.ascontiguousarray(hourly_weights, dtype=np.float64)
    if x0.dtype != np.float64:
        x0 = np.ascontiguousarray(x0, dtype=np.float64)

    cdef size_t n_hours  = <size_t>hourly_weights.shape[0]
    cdef size_t n        = <size_t>hourly_weights.shape[1]
    cdef size_t n2       = <size_t>hourly_weights.shape[2]

    cdef size_t n_fields = <size_t>x0.shape[0]
    cdef size_t n_x      = <size_t>x0.shape[1]

    cdef size_t take     = snapshots if snapshots <= steps else steps

    cdef PyRHSContext ctx
    ctx.fn_obj = <void*>rhs_fn
    ctx.params_obj = <void*>params

    cdef const double* H_ptr  = <const double*>hourly_weights.data
    cdef const double* x0_ptr = <const double*>x0.data

    cdef double* out_ptr = euler_solve_graph(
        H_ptr,
        x0_ptr,
        n, n_fields,
        n_hours,
        t0_hours,
        dt,
        steps,
        snapshots,
        <rhs_graph_fn>_rhs_graph_trampoline,
        <void*>&ctx
    )

    cdef cnp.ndarray out = np.empty((take, n_fields, n), dtype=np.float64)
    memcpy(<void*>out.data, <const void*>out_ptr,
           <size_t>take * n_fields * n * sizeof(DTYPE_t))
    free_array(out_ptr)
    return out

def euler_solve_graph_csr_rd(indptr_list, indices_list, data_list,
                             cnp.ndarray x0,
                             double dt, size_t steps, size_t snapshots,
                             rhs_fn, params=None):
    cdef size_t n_fields = <size_t>x0.shape[0]
    cdef size_t n        = <size_t>x0.shape[1]
    cdef size_t take     = snapshots if snapshots <= steps else steps

    cdef PyRHSContext ctx
    ctx.fn_obj = <void*>rhs_fn
    ctx.params_obj = <void*>params

    cdef Py_ssize_t n_steps = len(indptr_list)

    if steps > <size_t>n_steps:
        steps = <size_t>n_steps

    cdef const int   **indptr_ptrs  = <const int **>  malloc(n_steps * sizeof(const int *))
    cdef const int   **indices_ptrs = <const int **>  malloc(n_steps * sizeof(const int *))
    cdef const double**data_ptrs    = <const double**>malloc(n_steps * sizeof(const double *))


    cdef cnp.ndarray arr
    cdef Py_ssize_t i

    for i in range(n_steps):
        arr = <cnp.ndarray>indptr_list[i]
        indptr_ptrs[i] = <const int*>arr.data

        arr = <cnp.ndarray>indices_list[i]
        indices_ptrs[i] = <const int*>arr.data

        arr = <cnp.ndarray>data_list[i]
        data_ptrs[i] = <const double*>arr.data

    cdef const double* x0_ptr = <const double*> x0.data

    cdef double* out_ptr = euler_solve_graph_csr(
        indptr_ptrs, indices_ptrs, data_ptrs,
        x0_ptr,
        n, n_fields,
        dt, steps, snapshots,
        <rhs_graph_fn>_rhs_graph_trampoline,
        <void*>&ctx
    )

    free(<void*>indptr_ptrs)
    free(<void*>indices_ptrs)
    free(<void*>data_ptrs)

    cdef cnp.ndarray out = np.empty((take, n_fields, n), dtype=np.float64)
    memcpy(<void*>out.data, <const void*>out_ptr,
           <size_t>take * n_fields * n * sizeof(DTYPE_t))
    free_array(out_ptr)
    return out

