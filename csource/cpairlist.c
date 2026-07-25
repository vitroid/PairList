#include <Python.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>
#include <limits.h>
#include <stdlib.h>
#include "pairlist.h"

static PyObject *pairs(PyObject *self, PyObject *args);
static PyObject *pairs2(PyObject *self, PyObject *args);
static PyObject *pairs_filtered(PyObject *self, PyObject *args);
static PyObject *pairs2_filtered(PyObject *self, PyObject *args);
static PyObject *filtered_iter_new(PyObject *self, PyObject *args);
static PyObject *filtered2_iter_new(PyObject *self, PyObject *args);
static PyObject *filtered_iter_next(PyObject *self, PyObject *args);
static PyObject *filtered_iter_close(PyObject *self, PyObject *args);

static PyMethodDef module_methods[] = {
  {"pairs", pairs, METH_VARARGS,
   "Rough neighbor list from fractional coordinates and a grid."},
  {"pairs2", pairs2, METH_VARARGS,
   "Rough hetero neighbor list from fractional coordinates and a grid."},
  {"pairs_filtered", pairs_filtered, METH_VARARGS,
   "Neighbor pairs within rc (triclinic). Args: rpos, gx,gy,gz, cell(3,3), rc[, distance=0]."},
  {"pairs2_filtered", pairs2_filtered, METH_VARARGS,
   "Hetero neighbor pairs within rc (triclinic). Args: rpos0, rpos1, gx,gy,gz, cell, rc[, distance=0]."},
  {"filtered_iter_new", filtered_iter_new, METH_VARARGS,
   "Create chunk iterator. Args: rpos, gx,gy,gz, cell, rc"},
  {"filtered2_iter_new", filtered2_iter_new, METH_VARARGS,
   "Create hetero chunk iterator. Args: rpos0, rpos1, gx,gy,gz, cell, rc"},
  {"filtered_iter_next", filtered_iter_next, METH_VARARGS,
   "Next chunk. Args: iter, chunk_size[, distance=0]. Returns None when done."},
  {"filtered_iter_close", filtered_iter_close, METH_VARARGS,
   "Release chunk iterator."},
  {NULL, NULL, 0, NULL}};

static struct PyModuleDef moduledef = {
  PyModuleDef_HEAD_INIT,
  "cpairlist",
  NULL,
  -1,
  module_methods,
  NULL,
  NULL,
  NULL,
  NULL
};

PyMODINIT_FUNC PyInit_cpairlist(void) {
  /* import_array1: on failure returns NULL with ImportError set (no PyErr_Print). */
  import_array1(NULL);
  return PyModule_Create(&moduledef);
}


/* Convert obj to a C-contiguous float64 array of shape (N, 3).
   On success, returns a new reference (caller must Py_DECREF).
   On failure, sets a Python exception and returns NULL. */
static PyArrayObject *as_rpos_array(PyObject *obj, const char *name) {
  /* NPY_ARRAY_IN_ARRAY: C-contiguous + aligned (read-only OK).
     NPY_DOUBLE: cast to float64 when needed. */
  PyObject *arr = PyArray_FROM_OTF(obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
  if (arr == NULL) {
    return NULL;
  }
  if (PyArray_NDIM((PyArrayObject *)arr) != 2) {
    PyErr_Format(PyExc_ValueError,
                 "%s must be a 2-dimensional array of shape (N, 3), "
                 "got ndim=%d",
                 name, PyArray_NDIM((PyArrayObject *)arr));
    Py_DECREF(arr);
    return NULL;
  }
  if (PyArray_DIM((PyArrayObject *)arr, 1) != 3) {
    PyErr_Format(PyExc_ValueError,
                 "%s must have shape (N, 3), got shape (%" NPY_INTP_FMT
                 ", %" NPY_INTP_FMT ")",
                 name, PyArray_DIM((PyArrayObject *)arr, 0),
                 PyArray_DIM((PyArrayObject *)arr, 1));
    Py_DECREF(arr);
    return NULL;
  }
  return (PyArrayObject *)arr;
}


static PyArrayObject *as_cell_array(PyObject *obj) {
  PyObject *arr = PyArray_FROM_OTF(obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
  if (arr == NULL) {
    return NULL;
  }
  if (PyArray_NDIM((PyArrayObject *)arr) != 2 ||
      PyArray_DIM((PyArrayObject *)arr, 0) != 3 ||
      PyArray_DIM((PyArrayObject *)arr, 1) != 3) {
    PyErr_SetString(PyExc_ValueError, "cell must have shape (3, 3)");
    Py_DECREF(arr);
    return NULL;
  }
  return (PyArrayObject *)arr;
}


/* ngrid[d] must be >= 1 (matches pairlist.determine_grid). */
static int check_ngrid(const int ngrid[3]) {
  for (int d = 0; d < 3; d++) {
    if (ngrid[d] < 1) {
      PyErr_Format(PyExc_ValueError,
                   "ngrid[%d] must be >= 1, got %d", d, ngrid[d]);
      return -1;
    }
  }
  return 0;
}


/* Convert a malloc'd pair buffer into a numpy array.
   Pair indices are C int; NPY_INT is the dtype that matches sizeof(int).
   npairs==0 with pairs_data==NULL → empty (0,2) array.
   On failure after taking ownership path, frees pairs_data. */
static PyObject *pairs_to_ndarray(int npairs, int *pairs_data) {
  npy_intp dims[2] = {npairs, 2};
  if (npairs == 0) {
    /* do not pass NULL to SimpleNewFromData */
    return PyArray_SimpleNew(2, dims, NPY_INT);
  }
  PyObject *narray = PyArray_SimpleNewFromData(2, dims, NPY_INT, pairs_data);
  if (narray == NULL) {
    free(pairs_data);
    return NULL;
  }
  PyArray_ENABLEFLAGS((PyArrayObject *)narray, NPY_ARRAY_OWNDATA);
  return narray;
}


static PyObject *dists_to_ndarray(int npairs, double *dists_data) {
  npy_intp dims[1] = {npairs};
  if (npairs == 0) {
    return PyArray_SimpleNew(1, dims, NPY_DOUBLE);
  }
  PyObject *narray =
      PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, dists_data);
  if (narray == NULL) {
    free(dists_data);
    return NULL;
  }
  PyArray_ENABLEFLAGS((PyArrayObject *)narray, NPY_ARRAY_OWNDATA);
  return narray;
}


static PyObject *pairs(PyObject *self, PyObject *args) {
  PyObject *rpos_obj;
  int ngrid[3];

  if (!PyArg_ParseTuple(args, "Oiii", &rpos_obj, &ngrid[0], &ngrid[1],
                        &ngrid[2])) {
    return NULL;
  }
  if (check_ngrid(ngrid) < 0) {
    return NULL;
  }
  long long ngrid_prod =
      (long long)ngrid[0] * (long long)ngrid[1] * (long long)ngrid[2];
  if (ngrid_prod > INT_MAX) {
    PyErr_SetString(PyExc_ValueError, "ngrid product is too large");
    return NULL;
  }

  PyArrayObject *rpos = as_rpos_array(rpos_obj, "rpos");
  if (rpos == NULL) {
    return NULL;
  }
  npy_intp n_atoms = PyArray_DIM(rpos, 0);
  if (n_atoms > INT_MAX) {
    Py_DECREF(rpos);
    PyErr_SetString(PyExc_ValueError, "too many atoms for cpairlist");
    return NULL;
  }

  int n = (int)n_atoms;
  double *a = (double *)PyArray_DATA(rpos);
  int *pairs_buf = NULL;
  int npairs;
  Py_BEGIN_ALLOW_THREADS
  npairs = Pairs(n, a, ngrid, &pairs_buf);
  Py_END_ALLOW_THREADS
  Py_DECREF(rpos);
  if (npairs < 0) {
    PyErr_SetString(PyExc_MemoryError,
                    "failed to allocate neighbor pair list");
    return NULL;
  }
  return pairs_to_ndarray(npairs, pairs_buf);
}


static PyObject *pairs2(PyObject *self, PyObject *args) {
  PyObject *rpos0_obj, *rpos1_obj;
  int ngrid[3];

  if (!PyArg_ParseTuple(args, "OOiii", &rpos0_obj, &rpos1_obj, &ngrid[0],
                        &ngrid[1], &ngrid[2])) {
    return NULL;
  }
  if (check_ngrid(ngrid) < 0) {
    return NULL;
  }
  long long ngrid_prod =
      (long long)ngrid[0] * (long long)ngrid[1] * (long long)ngrid[2];
  if (ngrid_prod > INT_MAX) {
    PyErr_SetString(PyExc_ValueError, "ngrid product is too large");
    return NULL;
  }

  PyArrayObject *rpos0 = as_rpos_array(rpos0_obj, "rpos0");
  if (rpos0 == NULL) {
    return NULL;
  }
  PyArrayObject *rpos1 = as_rpos_array(rpos1_obj, "rpos1");
  if (rpos1 == NULL) {
    Py_DECREF(rpos0);
    return NULL;
  }
  npy_intp n0_atoms = PyArray_DIM(rpos0, 0);
  npy_intp n1_atoms = PyArray_DIM(rpos1, 0);
  if (n0_atoms > INT_MAX || n1_atoms > INT_MAX) {
    Py_DECREF(rpos0);
    Py_DECREF(rpos1);
    PyErr_SetString(PyExc_ValueError, "too many atoms for cpairlist");
    return NULL;
  }

  int n0 = (int)n0_atoms;
  int n1 = (int)n1_atoms;
  double *a0 = (double *)PyArray_DATA(rpos0);
  double *a1 = (double *)PyArray_DATA(rpos1);
  int *pairs_buf = NULL;
  int npairs;
  Py_BEGIN_ALLOW_THREADS
  npairs = Pairs2(n0, a0, n1, a1, ngrid, &pairs_buf);
  Py_END_ALLOW_THREADS
  Py_DECREF(rpos0);
  Py_DECREF(rpos1);
  if (npairs < 0) {
    PyErr_SetString(PyExc_MemoryError,
                    "failed to allocate neighbor pair list");
    return NULL;
  }
  return pairs_to_ndarray(npairs, pairs_buf);
}


static PyObject *pairs_filtered(PyObject *self, PyObject *args) {
  PyObject *rpos_obj, *cell_obj;
  int ngrid[3];
  double rc;
  int want_dist = 0;

  if (!PyArg_ParseTuple(args, "OiiiOd|i", &rpos_obj, &ngrid[0], &ngrid[1],
                        &ngrid[2], &cell_obj, &rc, &want_dist)) {
    return NULL;
  }
  if (check_ngrid(ngrid) < 0) {
    return NULL;
  }
  long long ngrid_prod =
      (long long)ngrid[0] * (long long)ngrid[1] * (long long)ngrid[2];
  if (ngrid_prod > INT_MAX) {
    PyErr_SetString(PyExc_ValueError, "ngrid product is too large");
    return NULL;
  }

  PyArrayObject *rpos = as_rpos_array(rpos_obj, "rpos");
  if (rpos == NULL) {
    return NULL;
  }
  PyArrayObject *cell = as_cell_array(cell_obj);
  if (cell == NULL) {
    Py_DECREF(rpos);
    return NULL;
  }
  npy_intp n_atoms = PyArray_DIM(rpos, 0);
  if (n_atoms > INT_MAX) {
    Py_DECREF(rpos);
    Py_DECREF(cell);
    PyErr_SetString(PyExc_ValueError, "too many atoms for cpairlist");
    return NULL;
  }

  int n = (int)n_atoms;
  double *a = (double *)PyArray_DATA(rpos);
  double *c = (double *)PyArray_DATA(cell);
  int *pairs_buf = NULL;
  double *dists_buf = NULL;
  int npairs;
  Py_BEGIN_ALLOW_THREADS
  npairs = PairsFiltered(n, a, ngrid, c, rc, &pairs_buf,
                         want_dist ? &dists_buf : NULL);
  Py_END_ALLOW_THREADS
  Py_DECREF(rpos);
  Py_DECREF(cell);
  if (npairs < 0) {
    PyErr_SetString(PyExc_MemoryError,
                    "failed to allocate filtered pair list");
    return NULL;
  }

  PyObject *parr = pairs_to_ndarray(npairs, pairs_buf);
  if (parr == NULL) {
    free(dists_buf);
    return NULL;
  }
  if (!want_dist) {
    return parr;
  }
  PyObject *darr = dists_to_ndarray(npairs, dists_buf);
  if (darr == NULL) {
    Py_DECREF(parr);
    return NULL;
  }
  return Py_BuildValue("NN", parr, darr);
}


static PyObject *pairs2_filtered(PyObject *self, PyObject *args) {
  PyObject *rpos0_obj, *rpos1_obj, *cell_obj;
  int ngrid[3];
  double rc;
  int want_dist = 0;

  if (!PyArg_ParseTuple(args, "OOiiiOd|i", &rpos0_obj, &rpos1_obj, &ngrid[0],
                        &ngrid[1], &ngrid[2], &cell_obj, &rc, &want_dist)) {
    return NULL;
  }
  if (check_ngrid(ngrid) < 0) {
    return NULL;
  }
  long long ngrid_prod =
      (long long)ngrid[0] * (long long)ngrid[1] * (long long)ngrid[2];
  if (ngrid_prod > INT_MAX) {
    PyErr_SetString(PyExc_ValueError, "ngrid product is too large");
    return NULL;
  }

  PyArrayObject *rpos0 = as_rpos_array(rpos0_obj, "rpos0");
  if (rpos0 == NULL) {
    return NULL;
  }
  PyArrayObject *rpos1 = as_rpos_array(rpos1_obj, "rpos1");
  if (rpos1 == NULL) {
    Py_DECREF(rpos0);
    return NULL;
  }
  PyArrayObject *cell = as_cell_array(cell_obj);
  if (cell == NULL) {
    Py_DECREF(rpos0);
    Py_DECREF(rpos1);
    return NULL;
  }
  npy_intp n0_atoms = PyArray_DIM(rpos0, 0);
  npy_intp n1_atoms = PyArray_DIM(rpos1, 0);
  if (n0_atoms > INT_MAX || n1_atoms > INT_MAX) {
    Py_DECREF(rpos0);
    Py_DECREF(rpos1);
    Py_DECREF(cell);
    PyErr_SetString(PyExc_ValueError, "too many atoms for cpairlist");
    return NULL;
  }

  int n0 = (int)n0_atoms;
  int n1 = (int)n1_atoms;
  double *a0 = (double *)PyArray_DATA(rpos0);
  double *a1 = (double *)PyArray_DATA(rpos1);
  double *c = (double *)PyArray_DATA(cell);
  int *pairs_buf = NULL;
  double *dists_buf = NULL;
  int npairs;
  Py_BEGIN_ALLOW_THREADS
  npairs = Pairs2Filtered(n0, a0, n1, a1, ngrid, c, rc, &pairs_buf,
                          want_dist ? &dists_buf : NULL);
  Py_END_ALLOW_THREADS
  Py_DECREF(rpos0);
  Py_DECREF(rpos1);
  Py_DECREF(cell);
  if (npairs < 0) {
    PyErr_SetString(PyExc_MemoryError,
                    "failed to allocate filtered pair list");
    return NULL;
  }

  PyObject *parr = pairs_to_ndarray(npairs, pairs_buf);
  if (parr == NULL) {
    free(dists_buf);
    return NULL;
  }
  if (!want_dist) {
    return parr;
  }
  PyObject *darr = dists_to_ndarray(npairs, dists_buf);
  if (darr == NULL) {
    Py_DECREF(parr);
    return NULL;
  }
  return Py_BuildValue("NN", parr, darr);
}


/* ---- chunk iterator (capsule) ------------------------------------------- */

typedef struct {
  PairChunkIter *it;
  PyObject *keep0; /* owned refs keeping rpos/cell buffers alive */
  PyObject *keep1;
  PyObject *keep_cell;
} FilteredIterCapsule;

static const char *filtered_iter_name = "pairlist.FilteredIter";

static void filtered_iter_destructor(PyObject *capsule) {
  FilteredIterCapsule *wrap =
      (FilteredIterCapsule *)PyCapsule_GetPointer(capsule, filtered_iter_name);
  if (wrap == NULL) {
    PyErr_Clear();
    return;
  }
  PairChunkIter_free(wrap->it);
  Py_XDECREF(wrap->keep0);
  Py_XDECREF(wrap->keep1);
  Py_XDECREF(wrap->keep_cell);
  free(wrap);
}

static PyObject *filtered_iter_new(PyObject *self, PyObject *args) {
  PyObject *rpos_obj, *cell_obj;
  int ngrid[3];
  double rc;
  if (!PyArg_ParseTuple(args, "OiiiOd", &rpos_obj, &ngrid[0], &ngrid[1],
                        &ngrid[2], &cell_obj, &rc)) {
    return NULL;
  }
  if (check_ngrid(ngrid) < 0) {
    return NULL;
  }
  PyArrayObject *rpos = as_rpos_array(rpos_obj, "rpos");
  if (rpos == NULL) {
    return NULL;
  }
  PyArrayObject *cell = as_cell_array(cell_obj);
  if (cell == NULL) {
    Py_DECREF(rpos);
    return NULL;
  }
  npy_intp n_atoms = PyArray_DIM(rpos, 0);
  if (n_atoms > INT_MAX) {
    Py_DECREF(rpos);
    Py_DECREF(cell);
    PyErr_SetString(PyExc_ValueError, "too many atoms for cpairlist");
    return NULL;
  }

  FilteredIterCapsule *wrap =
      (FilteredIterCapsule *)calloc(1, sizeof(FilteredIterCapsule));
  if (wrap == NULL) {
    Py_DECREF(rpos);
    Py_DECREF(cell);
    return PyErr_NoMemory();
  }
  wrap->it = PairChunkIter_create((int)n_atoms, (double *)PyArray_DATA(rpos),
                                  ngrid, (double *)PyArray_DATA(cell), rc);
  if (wrap->it == NULL) {
    free(wrap);
    Py_DECREF(rpos);
    Py_DECREF(cell);
    PyErr_SetString(PyExc_MemoryError, "failed to create filtered iterator");
    return NULL;
  }
  wrap->keep0 = (PyObject *)rpos;
  wrap->keep_cell = (PyObject *)cell;
  PyObject *cap = PyCapsule_New(wrap, filtered_iter_name, filtered_iter_destructor);
  if (cap == NULL) {
    filtered_iter_destructor(NULL); /* can't call easily; manual cleanup */
    PairChunkIter_free(wrap->it);
    Py_DECREF(rpos);
    Py_DECREF(cell);
    free(wrap);
    return NULL;
  }
  return cap;
}

static PyObject *filtered2_iter_new(PyObject *self, PyObject *args) {
  PyObject *rpos0_obj, *rpos1_obj, *cell_obj;
  int ngrid[3];
  double rc;
  if (!PyArg_ParseTuple(args, "OOiiiOd", &rpos0_obj, &rpos1_obj, &ngrid[0],
                        &ngrid[1], &ngrid[2], &cell_obj, &rc)) {
    return NULL;
  }
  if (check_ngrid(ngrid) < 0) {
    return NULL;
  }
  PyArrayObject *rpos0 = as_rpos_array(rpos0_obj, "rpos0");
  if (rpos0 == NULL) {
    return NULL;
  }
  PyArrayObject *rpos1 = as_rpos_array(rpos1_obj, "rpos1");
  if (rpos1 == NULL) {
    Py_DECREF(rpos0);
    return NULL;
  }
  PyArrayObject *cell = as_cell_array(cell_obj);
  if (cell == NULL) {
    Py_DECREF(rpos0);
    Py_DECREF(rpos1);
    return NULL;
  }
  npy_intp n0 = PyArray_DIM(rpos0, 0);
  npy_intp n1 = PyArray_DIM(rpos1, 0);
  if (n0 > INT_MAX || n1 > INT_MAX) {
    Py_DECREF(rpos0);
    Py_DECREF(rpos1);
    Py_DECREF(cell);
    PyErr_SetString(PyExc_ValueError, "too many atoms for cpairlist");
    return NULL;
  }

  FilteredIterCapsule *wrap =
      (FilteredIterCapsule *)calloc(1, sizeof(FilteredIterCapsule));
  if (wrap == NULL) {
    Py_DECREF(rpos0);
    Py_DECREF(rpos1);
    Py_DECREF(cell);
    return PyErr_NoMemory();
  }
  wrap->it = PairChunkIter_create2((int)n0, (double *)PyArray_DATA(rpos0),
                                   (int)n1, (double *)PyArray_DATA(rpos1),
                                   ngrid, (double *)PyArray_DATA(cell), rc);
  if (wrap->it == NULL) {
    free(wrap);
    Py_DECREF(rpos0);
    Py_DECREF(rpos1);
    Py_DECREF(cell);
    PyErr_SetString(PyExc_MemoryError, "failed to create filtered iterator");
    return NULL;
  }
  wrap->keep0 = (PyObject *)rpos0;
  wrap->keep1 = (PyObject *)rpos1;
  wrap->keep_cell = (PyObject *)cell;
  PyObject *cap = PyCapsule_New(wrap, filtered_iter_name, filtered_iter_destructor);
  if (cap == NULL) {
    PairChunkIter_free(wrap->it);
    Py_DECREF(rpos0);
    Py_DECREF(rpos1);
    Py_DECREF(cell);
    free(wrap);
    return NULL;
  }
  return cap;
}

static PyObject *filtered_iter_next(PyObject *self, PyObject *args) {
  PyObject *cap;
  int chunk_size;
  int want_dist = 0;
  if (!PyArg_ParseTuple(args, "Oi|i", &cap, &chunk_size, &want_dist)) {
    return NULL;
  }
  if (chunk_size < 1) {
    PyErr_SetString(PyExc_ValueError, "chunk_size must be >= 1");
    return NULL;
  }
  FilteredIterCapsule *wrap =
      (FilteredIterCapsule *)PyCapsule_GetPointer(cap, filtered_iter_name);
  if (wrap == NULL || wrap->it == NULL) {
    return NULL;
  }

  int *pairs_buf = (int *)malloc(sizeof(int) * (size_t)chunk_size * 2u);
  double *dists_buf = NULL;
  if (pairs_buf == NULL) {
    return PyErr_NoMemory();
  }
  if (want_dist) {
    dists_buf = (double *)malloc(sizeof(double) * (size_t)chunk_size);
    if (dists_buf == NULL) {
      free(pairs_buf);
      return PyErr_NoMemory();
    }
  }

  int n;
  Py_BEGIN_ALLOW_THREADS
  n = PairChunkIter_next(wrap->it, chunk_size, pairs_buf, dists_buf);
  Py_END_ALLOW_THREADS
  if (n < 0) {
    free(pairs_buf);
    free(dists_buf);
    PyErr_SetString(PyExc_RuntimeError, "filtered_iter_next failed");
    return NULL;
  }
  if (n == 0) {
    free(pairs_buf);
    free(dists_buf);
    Py_RETURN_NONE;
  }

  /* shrink buffers to exact n before handing to numpy */
  if (n < chunk_size) {
    int *p2 = (int *)realloc(pairs_buf, sizeof(int) * (size_t)n * 2u);
    if (p2) {
      pairs_buf = p2;
    }
    if (dists_buf) {
      double *d2 = (double *)realloc(dists_buf, sizeof(double) * (size_t)n);
      if (d2) {
        dists_buf = d2;
      }
    }
  }

  PyObject *parr = pairs_to_ndarray(n, pairs_buf);
  if (parr == NULL) {
    free(dists_buf);
    return NULL;
  }
  if (!want_dist) {
    return parr;
  }
  PyObject *darr = dists_to_ndarray(n, dists_buf);
  if (darr == NULL) {
    Py_DECREF(parr);
    return NULL;
  }
  return Py_BuildValue("NN", parr, darr);
}

static PyObject *filtered_iter_close(PyObject *self, PyObject *args) {
  PyObject *cap;
  if (!PyArg_ParseTuple(args, "O", &cap)) {
    return NULL;
  }
  FilteredIterCapsule *wrap =
      (FilteredIterCapsule *)PyCapsule_GetPointer(cap, filtered_iter_name);
  if (wrap == NULL) {
    PyErr_Clear();
    Py_RETURN_NONE;
  }
  PairChunkIter_free(wrap->it);
  wrap->it = NULL;
  Py_XDECREF(wrap->keep0);
  wrap->keep0 = NULL;
  Py_XDECREF(wrap->keep1);
  wrap->keep1 = NULL;
  Py_XDECREF(wrap->keep_cell);
  wrap->keep_cell = NULL;
  Py_RETURN_NONE;
}
