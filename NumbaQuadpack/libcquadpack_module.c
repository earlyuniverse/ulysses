#define PY_SSIZE_T_CLEAN
#include <Python.h>

/* Minimal Python extension stub.
 * All QUADPACK functions are accessed via ctypes in driver.py.
 * This file exists only to let setuptools compile and package the
 * C sources as a standard Python extension, which is then loaded
 * as a shared library at runtime.
 */
static PyMethodDef methods[] = { {NULL, NULL, 0, NULL} };
static struct PyModuleDef module = {
    PyModuleDef_HEAD_INIT, "libcquadpack", NULL, -1, methods
};
PyMODINIT_FUNC PyInit_libcquadpack(void) {
    return PyModule_Create(&module);
}
