#include "YODA/ROOTCnv.h"
#include "Python.h"
#include "YODA/Profile1D.h"
#include "TPython.h"

#ifndef ROOT_VERSION_CODE
#define ROOT_VERSION_CODE 397313
#endif

// #define XSTR(x) STR(x)
// #define STR(x) #x
// #pragma message "ROOT = " XSTR(ROOT_VERSION_CODE)


/// Get a PyROOT object from a ROOT one
inline PyObject* root_to_py_owned(TObject* root_obj) {
  // Different signatures in different ROOT versions... *sigh*
  #if ROOT_VERSION_CODE >= ROOT_VERSION(6,22,0)
  return TPython::CPPInstance_FromVoidPtr(root_obj, root_obj->ClassName());
  #elif ROOT_VERSION_CODE >= ROOT_VERSION(6,0,0)
  return TPython::ObjectProxy_FromVoidPtr(root_obj, root_obj->ClassName());
  #else
  return TPython::ObjectProxy_FromVoidPtr(root_obj, root_obj->ClassName(), kFALSE);
  #endif
}


/// Get the native ROOT object from a PyROOT wrapper
inline TObject* py_owned_to_root(PyObject* pyroot_obj) {
  // Different signatures in different ROOT versions... *sigh*
  #if ROOT_VERSION_CODE >= ROOT_VERSION(6,22,0)
  return (TObject*) TPython::CPPInstance_AsVoidPtr(pyroot_obj);
  #else
  return (TObject*) TPython::ObjectProxy_AsVoidPtr(pyroot_obj);
  #endif
}
