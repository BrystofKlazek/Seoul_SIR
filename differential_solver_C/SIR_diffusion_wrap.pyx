# distutils: language = c
# cython: boundscheck=False, wraparound=False, cdivision=True
cimport numpy as cnp
import numpy as np
from libc.string cimport memcpy
from libc.stddef cimport size_t

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

    double* euler_solve_graph(const double *L, const double *x0,
                              size_t vec_size, size_t n_fields,
                              double dt, size_t steps, size_t snapshots,
                              rhs_graph_fn rhs, void* userdata)

    # --- CSR variant (int indices) ---
    double* euler_solve_graph_csr(const int *indptr, const int *indices, const double *data,
                                  const double *x0,
                                  size_t vec_size, size_t n_fields,
                                  double dt, size_t steps, size_t snapshots,
                                  rhs_graph_fn rhs, void* userdata)

    void free_array(double* p)

cdef struct PyRHSContext:
    void* fn_obj
    void* params_obj

cdef void _rhs_grid_trampoline(double* out, const double* state,
                               size_t n_fields, size_t rows, size_t cols,
                               void* userdata):
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

cdef void _rhs_graph_trampoline(double* out, const double* state,
                                size_t n_fields, size_t vec_size,
                                void* userdata):
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

def euler_solve_graph_rd(cnp.ndarray L, cnp.ndarray x0,
                         double dt, size_t steps, size_t snapshots,
                         rhs_fn, params=None):
    cdef size_t n_fields = <size_t>L.shape[0]
    cdef size_t n        = <size_t>L.shape[1]
    cdef size_t take     = snapshots if snapshots <= steps else steps

    cdef PyRHSContext ctx
    ctx.fn_obj = <void*>rhs_fn
    ctx.params_obj = <void*>params

    cdef const double* L_ptr  = <const double*> L.data
    cdef const double* x0_ptr = <const double*> x0.data

    cdef double* out_ptr = euler_solve_graph(L_ptr, x0_ptr,
                                             n, n_fields,
                                             dt, steps, snapshots,
                                             <rhs_graph_fn>_rhs_graph_trampoline,
                                             <void*>&ctx)

    cdef cnp.ndarray out = np.empty((take, n_fields, n), dtype=np.float64)
    memcpy(<void*>out.data, <const void*>out_ptr,
           <size_t>take * n_fields * n * sizeof(DTYPE_t))
    free_array(out_ptr)
    return out

def euler_solve_graph_csr_rd(cnp.ndarray indptr, cnp.ndarray indices, cnp.ndarray data,
                             cnp.ndarray x0,
                             double dt, size_t steps, size_t snapshots,
                             rhs_fn, params=None):
    cdef size_t n_fields = <size_t>x0.shape[0]
    cdef size_t n        = <size_t>x0.shape[1]
    cdef size_t take     = snapshots if snapshots <= steps else steps

    cdef PyRHSContext ctx
    ctx.fn_obj = <void*>rhs_fn
    ctx.params_obj = <void*>params

    cdef const int* indptr_ptr  = <const int*> indptr.data
    cdef const int* indices_ptr = <const int*> indices.data
    cdef const double* data_ptr = <const double*> data.data
    cdef const double* x0_ptr   = <const double*> x0.data

    cdef double* out_ptr = euler_solve_graph_csr(indptr_ptr, indices_ptr, data_ptr,
                                                 x0_ptr,
                                                 n, n_fields,
                                                 dt, steps, snapshots,
                                                 <rhs_graph_fn>_rhs_graph_trampoline,
                                                 <void*>&ctx)

    cdef cnp.ndarray out = np.empty((take, n_fields, n), dtype=np.float64)
    memcpy(<void*>out.data, <const void*>out_ptr,
           <size_t>take * n_fields * n * sizeof(DTYPE_t))
    free_array(out_ptr)
    return out

