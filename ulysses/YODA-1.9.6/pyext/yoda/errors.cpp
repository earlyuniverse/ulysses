#include "Python.h"

#include <exception>
#include <stdexcept>
#include <iostream>
#include <typeinfo>

#include "core.h"
#include "YODA/Exceptions.h"
#include "pyerrors.h"

void translate_yoda_error() {
  try {
    if (PyErr_Occurred())
      ; // let the latest Python exn pass through and ignore the current one
    else
      throw;
  } catch (const std::out_of_range& exn) {
    PyErr_SetString(PyExc_IndexError, exn.what());
  } catch (const std::bad_cast& exn) {
    PyErr_SetString(PyExc_TypeError, exn.what());
  } catch (const std::domain_error& exn) {
    PyErr_SetString(PyExc_ValueError, exn.what());
  } catch (const std::invalid_argument& exn) {
    PyErr_SetString(PyExc_ValueError, exn.what());
  } catch (const std::ios_base::failure& exn) {
    PyErr_SetString(PyExc_IOError, exn.what());
  } catch (const std::overflow_error& exn) {
    PyErr_SetString(PyExc_OverflowError, exn.what());
  } catch (const std::range_error& exn) {
    PyErr_SetString(PyExc_ArithmeticError, exn.what());
  } catch (const std::underflow_error& exn) {
    PyErr_SetString(PyExc_ArithmeticError, exn.what());
  } catch (const std::bad_alloc& exn) {
    PyErr_SetString(PyExc_MemoryError, exn.what());
  } catch (const YODA::BinningError& exn) {
    PyErr_SetString(YodaExc_BinningError, exn.what());
  } catch (const YODA::RangeError& exn) {
    PyErr_SetString(YodaExc_RangeError, exn.what());
  } catch (const YODA::LockError& exn) {
    PyErr_SetString(YodaExc_LockError, exn.what());
  } catch (const YODA::GridError& exn) {
    PyErr_SetString(YodaExc_GridError, exn.what());
  } catch (const YODA::LogicError& exn) {
    PyErr_SetString(YodaExc_LogicError, exn.what());
  } catch (const YODA::WeightError& exn) {
    PyErr_SetString(YodaExc_WeightError, exn.what());
  } catch (const YODA::LowStatsError& exn) {
    PyErr_SetString(YodaExc_LowStatsError, exn.what());
  } catch (const YODA::AnnotationError& exn) {
    PyErr_SetString(YodaExc_AnnotationError, exn.what());
  } catch (const YODA::ReadError& exn) {
    PyErr_SetString(YodaExc_ReadError, exn.what());
  } catch (const YODA::UserError& exn) {
    PyErr_SetString(YodaExc_UserError, exn.what());
  } catch (const YODA::Exception& exn) {
    PyErr_SetString(YodaExc_Exception, exn.what());
  } catch (...) {
    PyErr_SetString(PyExc_RuntimeError, "Unknown exception");
  }
}
