//http://acooke.org/cute/ExampleCod0.html
//https://github.com/numpy/numpy/blob/master/numpy/core/src/dummymodule.c
//https://qiita.com/junkoda/items/17df11d7a20dc9d50e7d




#include <Python.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>
//#include "C_arraytest.h"
#include <math.h>
#include "pairlist.h"

static PyObject *pairs(PyObject *self, PyObject* args);
static PyObject *pairs2(PyObject *self, PyObject* args);

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


//my initializer
PyMODINIT_FUNC PyInit_cpairlist(void) {
  PyObject *m;
  import_array();
  m = PyModule_Create(&moduledef);
  if (!m)
    return NULL;
  return m;
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


static PyObject *pairs(PyObject *self, PyObject *args) {
  PyObject *rpos_obj;
  int ngrid[3];

  if (!PyArg_ParseTuple(args, "Oiii", &rpos_obj, &ngrid[0], &ngrid[1],
                        &ngrid[2])) {
    return NULL;
  }

  PyArrayObject *rpos = as_rpos_array(rpos_obj, "rpos");
  if (rpos == NULL) {
    return NULL;
  }

  int n = (int)PyArray_DIM(rpos, 0);
  double *a = (double *)PyArray_DATA(rpos);
  int *pairs;
  int npairs = Pairs(n, a, ngrid, &pairs);

  /* return the array as a numpy array (numpy will free it later) */
  npy_intp output_dims[2] = {npairs, 2};
  PyObject *narray = PyArray_SimpleNewFromData(2, output_dims, NPY_INT, pairs);
  /* this is the critical line - tell numpy it has to free the data */
  PyArray_ENABLEFLAGS((PyArrayObject *)narray, NPY_ARRAY_OWNDATA);
  Py_DECREF(rpos);
  return narray;
}


static PyObject *pairs2(PyObject *self, PyObject *args) {
  PyObject *rpos0_obj, *rpos1_obj;
  int ngrid[3];

  if (!PyArg_ParseTuple(args, "OOiii", &rpos0_obj, &rpos1_obj, &ngrid[0],
                        &ngrid[1], &ngrid[2])) {
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

  int n0 = (int)PyArray_DIM(rpos0, 0);
  int n1 = (int)PyArray_DIM(rpos1, 0);
  double *a0 = (double *)PyArray_DATA(rpos0);
  double *a1 = (double *)PyArray_DATA(rpos1);
  int *pairs;
  int npairs = Pairs2(n0, a0, n1, a1, ngrid, &pairs);

  /* return the array as a numpy array (numpy will free it later) */
  npy_intp output_dims[2] = {npairs, 2};
  PyObject *narray = PyArray_SimpleNewFromData(2, output_dims, NPY_INT, pairs);
  /* this is the critical line - tell numpy it has to free the data */
  PyArray_ENABLEFLAGS((PyArrayObject *)narray, NPY_ARRAY_OWNDATA);
  Py_DECREF(rpos0);
  Py_DECREF(rpos1);
  return narray;
}
