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

    double* euler_solve_graph(const double *L, const double *x0,
                              size_t vec_size, size_t n_fields,
                              double dt, size_t steps, size_t snapshots,
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

cdef struct PyRHSContext:
    void* fn_obj
    void* params_obj

# --- grid trampoline ---

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

# --- graph trampoline ---

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

def euler_solve_graph_csr_rd(indptr_list, indices_list, data_list,
                             cnp.ndarray x0,
                             double dt, size_t steps, size_t snapshots,
                             rhs_fn, params=None):
    """
    Run the CSR graph solver with time-dependent Laplacians.

    Parameters
    ----------
    indptr_list, indices_list, data_list : sequence of length n_steps
        Each element is a 1D NumPy array:
            indptr_list[k].dtype  == int32, shape (n+1,)
            indices_list[k].dtype == int32, shape (nnz_k,)
            data_list[k].dtype    == float64, shape (nnz_k,)
    x0 : ndarray, shape (n_fields, n)
        Initial state [S, I, R, ...].
    dt, steps, snapshots : as in the C solver.
    rhs_fn : Python RHS callback.
    params : optional parameters for rhs_fn.
    """
    cdef size_t n_fields = <size_t>x0.shape[0]
    cdef size_t n        = <size_t>x0.shape[1]
    cdef size_t take     = snapshots if snapshots <= steps else steps

    cdef PyRHSContext ctx
    ctx.fn_obj = <void*>rhs_fn
    ctx.params_obj = <void*>params

    # how many Laplacians do we have?
    cdef Py_ssize_t n_steps = len(indptr_list)
    if len(indices_list) != n_steps or len(data_list) != n_steps:
        raise ValueError("indptr_list, indices_list, data_list must have same length")

    # we must not ask the C solver for more steps than we have matrices for
    if steps > <size_t>n_steps:
        steps = <size_t>n_steps

    # allocate C arrays of pointers
    cdef const int   **indptr_ptrs  = <const int **>  malloc(n_steps * sizeof(const int *))
    cdef const int   **indices_ptrs = <const int **>  malloc(n_steps * sizeof(const int *))
    cdef const double**data_ptrs    = <const double**>malloc(n_steps * sizeof(const double *))

    if indptr_ptrs == NULL or indices_ptrs == NULL or data_ptrs == NULL:
        if indptr_ptrs  != NULL: free(<void*>indptr_ptrs)
        if indices_ptrs != NULL: free(<void*>indices_ptrs)
        if data_ptrs    != NULL: free(<void*>data_ptrs)
        raise MemoryError()

    cdef cnp.ndarray arr
    cdef Py_ssize_t i

    # fill pointer arrays from Python lists
    for i in range(n_steps):
        arr = <cnp.ndarray>indptr_list[i]
        if arr.ndim != 1:
            free(<void*>indptr_ptrs); free(<void*>indices_ptrs); free(<void*>data_ptrs)
            raise ValueError("each indptr must be 1D array")
        indptr_ptrs[i] = <const int*>arr.data

        arr = <cnp.ndarray>indices_list[i]
        if arr.ndim != 1:
            free(<void*>indptr_ptrs); free(<void*>indices_ptrs); free(<void*>data_ptrs)
            raise ValueError("each indices must be 1D array")
        indices_ptrs[i] = <const int*>arr.data

        arr = <cnp.ndarray>data_list[i]
        if arr.ndim != 1:
            free(<void*>indptr_ptrs); free(<void*>indices_ptrs); free(<void*>data_ptrs)
            raise ValueError("each data must be 1D array")
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

