#include <Python.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>
#include <limits.h>
#include <stdlib.h>
#include "pairlist.h"

static PyObject *pairs(PyObject *self, PyObject *args);
static PyObject *pairs2(PyObject *self, PyObject *args);

static PyMethodDef module_methods[] = {
  {"pairs", pairs, METH_VARARGS,
   "Rough neighbor list from given fractional coordinates of the particles and a grid."},
  {"pairs2", pairs2, METH_VARARGS,
   "Rough neighbor list from given two fractional coordinates of the different particles and a grid."},
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
  int *pairs_buf;
  int npairs = Pairs(n, a, ngrid, &pairs_buf);
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
  int *pairs_buf;
  int npairs = Pairs2(n0, a0, n1, a1, ngrid, &pairs_buf);
  Py_DECREF(rpos0);
  Py_DECREF(rpos1);
  if (npairs < 0) {
    PyErr_SetString(PyExc_MemoryError,
                    "failed to allocate neighbor pair list");
    return NULL;
  }
  return pairs_to_ndarray(npairs, pairs_buf);
}
