//http://acooke.org/cute/ExampleCod0.html
//https://github.com/numpy/numpy/blob/master/numpy/core/src/dummymodule.c
//https://qiita.com/junkoda/items/17df11d7a20dc9d50e7d




#include <Python.h>
//#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
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



//taken from C_arraytest.c in Scipy.

/* ==== Check that PyArrayObject is a double (Float) type and a matrix ==============
    return 1 if an error and raise exception */
int  not_doublematrix(PyArrayObject *mat)  {
	if (mat->descr->type_num != NPY_DOUBLE || mat->nd != 2)  {
		PyErr_SetString(PyExc_ValueError,
			"In not_doublematrix: array must be of type Float and 2 dimensional (n x m).");
		return 1;  }
	return 0;
}
/* ==== Allocate a double *vector (vec of pointers) ======================
    Memory is Allocated!  See void free_Carray(double ** )                  */
double **ptrvector(long n)  {
	double **v;
	v=(double **)malloc((size_t) (n*sizeof(double)));
	if (!v)   {
		printf("In **ptrvector. Allocation of memory for double array failed.");
		exit(0);  }
	return v;
}


/* ==== Create Carray from PyArray ======================
    Assumes PyArray is contiguous in memory.
    Memory is allocated!                                    */
double **pymatrix_to_Carrayptrs(PyArrayObject *arrayin)  {
	double **c, *a;
	int i,n,m;

	n=arrayin->dimensions[0];
	m=arrayin->dimensions[1];
	c=ptrvector(n);
	a=(double *) arrayin->data;  // pointer to arrayin data as double
	for ( i=0; i<n; i++)  {
		c[i]=a+i*m;  }
	return c;
}


/* ==== Free a double *vector (vec of pointers) ========================== */
void free_Carrayptrs(double **v)  {
	free((char*) v);
}



static PyObject *pairs(PyObject *self, PyObject* args) {
  // expect two arguments.
  PyArrayObject *rpos;
  int dimss[2], ngrid[3];

  /* Parse tuples separately since args will differ between C fcns */
  if (!PyArg_ParseTuple(args, "O!iii", &PyArray_Type,
			&rpos, &ngrid[0], &ngrid[1], &ngrid[2])) return NULL;
  //if (!PyArg_ParseTuple(args, "I", &n, &m, &ngrid[0], &ngrid[1], &ngrid[2])) return NULL;
  if (NULL == rpos) return NULL;
  if (not_doublematrix(rpos)) return NULL;

  /* Get the dimensions of the input */
  int n=dimss[0]=rpos->dimensions[0];
  //int m=dimss[1]=rpos->dimensions[1];

  double* a = (double*)rpos->data;
  int* pairs;
  int npairs = Pairs(n, a, ngrid, &pairs);
  //pairs is allocated.

  // return the array as a numpy array (numpy will free it later)
  npy_intp dims[2] = {npairs,2};
  PyObject *narray = PyArray_SimpleNewFromData(2, dims, NPY_INT, pairs);
  // this is the critical line - tell numpy it has to free the data
  PyArray_ENABLEFLAGS((PyArrayObject*)narray, NPY_ARRAY_OWNDATA);
  return narray;
}



static PyObject *pairs2(PyObject *self, PyObject* args) {
  // expect two arguments.
  PyArrayObject *rpos0, *rpos1;
  int ngrid[3];

  /* Parse tuples separately since args will differ between C fcns */
  if (!PyArg_ParseTuple(args, "O!O!iii", &PyArray_Type, &rpos0,
			&PyArray_Type, &rpos1, &ngrid[0], &ngrid[1], &ngrid[2])) return NULL;
  if (NULL == rpos0) return NULL;
  if (NULL == rpos1) return NULL;
  if (not_doublematrix(rpos0)) return NULL;
  if (not_doublematrix(rpos1)) return NULL;

  /* Get the dimensions of the input */
  int n0=rpos0->dimensions[0];
  //int m0=rpos0->dimensions[1];
  int n1=rpos1->dimensions[0];
  //int m1=rpos1->dimensions[1];

  double* a0 = (double*)rpos0->data;
  double* a1 = (double*)rpos1->data;
  int* pairs;
  int npairs = Pairs2(n0, a0, n1, a1, ngrid, &pairs);
  //pairs is allocated.

  // return the array as a numpy array (numpy will free it later)
  npy_intp dims[2] = {npairs,2};
  PyObject *narray = PyArray_SimpleNewFromData(2, dims, NPY_INT, pairs);
  // this is the critical line - tell numpy it has to free the data
  PyArray_ENABLEFLAGS((PyArrayObject*)narray, NPY_ARRAY_OWNDATA);
  return narray;
}
