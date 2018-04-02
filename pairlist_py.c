//http://acooke.org/cute/ExampleCod0.html
//https://github.com/numpy/numpy/blob/master/numpy/core/src/dummymodule.c
//https://qiita.com/junkoda/items/17df11d7a20dc9d50e7d




#include <Python.h>
// i think this means we're using the 1.9 api
// http://docs.scipy.org/doc/numpy/reference/c-api.deprecations.html
#define NPY_NO_DEPRECATED_API NPY_1_9_API_VERSION
#include <numpy/arrayobject.h>


static PyObject *new_vec(PyObject *self, PyObject *args);

static PyMethodDef module_methods[] = {
  {"new_vec", new_vec, METH_VARARGS,
   "create a new vector with n entries, numbered 0 to n-1"},
  {NULL, NULL, 0, NULL}};

static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "dummy",
        NULL,
        -1,
        module_methods,
        NULL,
        NULL,
        NULL,
        NULL
};


PyMODINIT_FUNC PyInit_pairlist_py(void) {
  PyObject *m;
  import_array();
  m = PyModule_Create(&moduledef);
  if (!m)
    return NULL;
  return m;
}


static PyObject *new_vec(PyObject *self, PyObject *args) {
  // expect a single integer argument
  int i, n;
  if (!(PyArg_ParseTuple(args, "i", &n))) return NULL;
  // create the array in C on the heap
  int *array = NULL;
  if (!(array = malloc(n * sizeof(int)))) return NULL;
  for (i = 0; i < n; ++i) array[i] = i;
  // return the array as a numpy array (numpy will free it later)
  npy_intp dims[1] = {n};
  PyObject *narray = PyArray_SimpleNewFromData(1, dims, NPY_INT, array);
  // this is the critical line - tell numpy it has to free the data
  PyArray_ENABLEFLAGS((PyArrayObject*)narray, NPY_ARRAY_OWNDATA);
  return narray;
}
