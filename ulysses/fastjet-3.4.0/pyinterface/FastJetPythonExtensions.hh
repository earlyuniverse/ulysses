#ifndef __FASTJET_PYTHONUSERINFO_HH__
#include "fastjet/PseudoJet.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/Selector.hh"
#include "fastjet/Error.hh"
#include "Python.h"

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

//----------------------------------------------------------------------
/// \class UserInfoPython
/// Internal helper class for making user info usable within python
///
/// This is an internal class that makes possible the calls to
///   pseudojet->set_python_info(object)
///   object = pseudojet->python_info()
/// through the PseudoJet::UserInfoBase interface and for any object
/// class in python
class UserInfoPython : public fastjet::PseudoJet::UserInfoBase {
public:
  UserInfoPython(PyObject * pyobj) : _pyobj(pyobj) {
    Py_XINCREF(_pyobj);
  }

  PyObject * get_pyobj() const {
    // since there's going to be an extra reference to this object
    // one must increase the reference count; it seems that this
    // is _our_ responsibility
    Py_XINCREF(_pyobj);
    return _pyobj;
  }
  //const PyObject * get_pyobj() const {return _pyobj;}
  
  ~UserInfoPython() {
    Py_XDECREF(_pyobj);
  }
private:
  PyObject * _pyobj;
};

//----------------------------------------------------------------------
/// get a C++ string from a Python object that is assumed to be
/// of string type (unicode for Py3).
///
/// Later we might imagine moving to using something like
/// SWIG_AsPtr_std_string.
///
/// Implementation with Python calls is inspired from discussions
/// at
///   https://stackoverflow.com/questions/22487780/what-do-i-use-instead-of-pystring-asstring-when-loading-a-python-module-in-3-3
///   https://mail.python.org/pipermail/python-list/2009-March/527813.html
///
inline std::string cpp_string_from_py_str(PyObject *py_str) {
  const char *char_result;
#if PY_VERSION_HEX >= 0x03000000
  char_result = PyUnicode_AsUTF8(py_str);
#else
  char_result = PyString_AsString(py_str);
#endif
  return std::string(char_result);
}

//----------------------------------------------------------------------
/// Invokes the str call on a python object and returns the corresponding
/// C++ string
inline std::string cpp_string_from_str_py_obj(PyObject *py_obj) {
  PyObject* py_str = PyObject_Str(py_obj);
  std::string cpp_str = cpp_string_from_py_str(py_str);
  Py_XDECREF(py_str);
  return cpp_str;
}

//----------------------------------------------------------------------
/// Invokes the name call on a python object and returns the corresponding
/// C++ string
inline std::string cpp_string_from_name_py_obj(PyObject *py_obj) {
  PyObject* py_str = PyObject_GetAttrString(py_obj, "__name__");
  std::string cpp_str = cpp_string_from_py_str(py_str);
  Py_XDECREF(py_str);
  return cpp_str;
}


//----------------------------------------------------------------------
/// \class SelectorWorkerPython
/// Internal class for making python classes/functions usable as selectors
///
/// This is an internal class that makes possible the calls to
///   selector = SelectorPython(pyton_function)
/// where python_function will take a PseudoJet as argument and return
/// a bool
class SelectorWorkerPython : public SelectorWorker{
public:
  // ctor based on a PYObject which should be callable (typically a
  // class or a function)
  SelectorWorkerPython(PyObject *py_class_or_function) : _py_class_or_function(py_class_or_function){
    // increment ref count on the py object
    Py_XINCREF(_py_class_or_function);

    // we directly make sure that the function is callable
    if (!PyCallable_Check(_py_class_or_function)){
      PyErr_SetString(PyExc_TypeError,
          "SelectorWorkerPython::SelectorWorkerPython: the argument should be callable");
      // do we also throw a fastjet error?
    }
  }

  // dtor
  ~SelectorWorkerPython(){
    // decrement ref count on the py object
    Py_XDECREF(_py_class_or_function);
  }

  // description of the Selector
  virtual std::string description() const{
    // Functions define __name__ which gives a more readable output
    // that __str__. So we'll use it for the description
    if (PyObject_HasAttrString(_py_class_or_function, "__name__")){
      //Py_XINCREF(_py_class_or_function); // GPS: not needed?
      std::string cpp_str = cpp_string_from_name_py_obj(_py_class_or_function);
      //Py_XDECREF(_py_class_or_function); // GPS: not needed?
      return std::string("Selector based on python function ")+cpp_str;
    }
    
    // reuse the python string if it is available
    //
    // Note: this will take the __str__ method in classes but it would
    // also be available for functions, producing a rather inelegant
    // output of the form
    //   Selector based on python condition <function is_pileup at 0x...>
    // Not sure how to avoid this?
    if (PyObject_HasAttrString(_py_class_or_function, "__str__")){
      //Py_XINCREF(_py_class_or_function); // GPS: not needed?
      std::string cpp_str = cpp_string_from_str_py_obj(_py_class_or_function);
      //Py_XDECREF(_py_class_or_function); // GPS: not needed?
      return std::string("Selector based on python condition ")+cpp_str;
    }
    return "Selector based on python function";
  }

  // implement the check whether the given PseudoJet passes the
  // selection or not
  virtual bool pass(const PseudoJet &jet) const{
    // first make a copy of the jet in a PyObject* managed by swig
    PyObject *py_jet = 0;
    py_jet = SWIG_NewPointerObj((new fastjet::PseudoJet(static_cast< const fastjet::PseudoJet& >(jet))), SWIGTYPE_p_fastjet__PseudoJet, SWIG_POINTER_OWN |  0 );

    // now call the user-defined selection function (in python)
    // GPS: are the INCREF and DECREF really needed here?
    // GS:  I don't think so but I am unsure so I preferred to play it safe.
    Py_XINCREF(_py_class_or_function);
    PyObject * args = Py_BuildValue("(O)", py_jet);
    PyObject *py_result = PyObject_CallObject(_py_class_or_function, args);
    Py_XDECREF(_py_class_or_function);

    // and interpret the result as a bool
    //
    // Note that somehow the conversion from bool via SWIG_AsVal_bool
    // is not available at this stage so we cannot simply do:
    //   bool result;
    //   int conversion_result = SWIG_AsVal_bool(py_result, &result);
    //   if (!SWIG_IsOK(conversion_result)){
    //     throw Error("SelectorWorkerPython::pass(): the value returned by the python function could not be casted to a bool");
    //   }
    if (py_result == NULL)
      throw Error("SelectorWorkerPython::pass(): call to python function returned a NULL result.");

    if (!PyBool_Check(py_result))
      throw Error("SelectorWorkerPython::pass(): the value returned by the python function could not be cast to a bool");
    int result = PyObject_IsTrue(py_result);
    if (result == -1)
      throw Error("SelectorWorkerPython::pass(): the value returned by the python function could not be cast to a bool");
    Py_XDECREF(py_result); ///???

    return result;
  }
  
private:
  PyObject *_py_class_or_function;
};

// effectively create a Selector for python
Selector SelectorPython(PyObject *py_function_or_class) {
  return Selector(new SelectorWorkerPython(py_function_or_class));
}

//----------------------------------------------------------------------
/// \class RecombinerPython
/// Class allowing user-defined Recombiners in python
///
/// If a (python) user implements a (python) class providing the 
///   __str__()
///   PseudoJet recombine(PseudoJet pa, PseudoJe pb)
///   PseudoJet preprocess(PseudoJet pa)
/// methods
///
/// Note that compared to the C++ implementation, the result is
/// returned by the method rather than being passed as a reference.
class RecombinerPython : public JetDefinition::Recombiner{
public:
  /// ctor with the python recombier class as an argument
  RecombinerPython(PyObject *py_class) : _py_class(py_class){
    Py_XINCREF(_py_class);
    
    // here we could add some tests that the class has the required
    // methods
  }
  
  /// dtor
  virtual ~RecombinerPython(){
    Py_XDECREF(_py_class);
  }

  /// return a textual description of the recombiner
  virtual std::string description() const{
    if (! PyObject_HasAttrString(_py_class, "__str__")){
      throw Error("RecombinerPython: the provided class should implement the __str__ method (for description");
    }
    
    //Py_XINCREF(_py_class); // GPS not needed
    std::string cpp_str = cpp_string_from_str_py_obj(_py_class);
    //Py_XDECREF(_py_class); // GPS not needed
    return std::string("User-defined recombiner based on python recombiner ")+cpp_str;
  }
  
  /// recombine pa and pb and put result into pab
  virtual void recombine(const PseudoJet & pa, const PseudoJet & pb, 
                         PseudoJet & pab) const{
    // first make a copy of the "input" arguments as PyObject* managed by swig
    //PseudoJet pa_copy = pa;  // not sure this is needed
    PyObject *py_pa = 0;
    py_pa = SWIG_NewPointerObj((new fastjet::PseudoJet(static_cast< const fastjet::PseudoJet& >(pa))), SWIGTYPE_p_fastjet__PseudoJet, SWIG_POINTER_OWN |  0 );

    //PseudoJet pb_copy = pb;  // not sure this is needed
    PyObject *py_pb = 0;
    py_pb = SWIG_NewPointerObj((new fastjet::PseudoJet(static_cast< const fastjet::PseudoJet& >(pb))), SWIGTYPE_p_fastjet__PseudoJet, SWIG_POINTER_OWN |  0 );

    // now call the recombiner
    Py_XINCREF(_py_class);
    PyObject *py_result = PyObject_CallMethod(_py_class, (char *) "recombine",
                                              (char *) "(OO)", py_pa, py_pb);
    Py_XDECREF(_py_class);
    if (py_result == NULL)
      throw Error("RecombinerPython::recombine(): call to python function returned a NULL result.");

    // put the result in pab
    void *pab_void_ptr = 0;
    PseudoJet *pab_ptr = 0;
    int res1 = SWIG_ConvertPtr(py_result, &pab_void_ptr, SWIGTYPE_p_fastjet__PseudoJet, 0 );
    if (!SWIG_IsOK(res1)) {
      throw Error("RecombinerPython::recombine(): cannot reinterpret the last argument as a fastjet::PseudoJet.");
    }
    pab_ptr = reinterpret_cast< fastjet::PseudoJet * >(pab_void_ptr);
    pab = *pab_ptr;
    Py_XDECREF(py_result);  ///needed???
  }

  /// routine called to preprocess each input jet (to make all input
  /// jets compatible with the scheme requirements (e.g. massless).
  virtual void preprocess(PseudoJet & pa) const {
    // first make a copy of the arguments as PyObject* managed by swig
    //PseudoJet pa_copy = pa;  // not sure this is needed
    PyObject *py_pa = 0;
    py_pa = SWIG_NewPointerObj((new fastjet::PseudoJet(static_cast< fastjet::PseudoJet& >(pa))), SWIGTYPE_p_fastjet__PseudoJet, SWIG_POINTER_OWN |  0 );

    // then call the user-defined python function
    Py_XINCREF(_py_class);
    PyObject *py_result = PyObject_CallMethod(_py_class, (char *) "preprocess",
                                              (char *) "(O)", py_pa);
    Py_XDECREF(_py_class);
    if (py_result == NULL)
      throw Error("RecombinerPython::preprocess(): call to python function returned a NULL result.");

    // copy the result back in pa
    void *pa_void_ptr = 0;
    PseudoJet *pa_ptr = 0;
    int res1 = SWIG_ConvertPtr(py_pa, &pa_void_ptr, SWIGTYPE_p_fastjet__PseudoJet, 0);
    if (!SWIG_IsOK(res1)) {
      throw Error("RecombinerPython::preprocess(): cannot reinterpret the last argument as a fastjet::PseudoJet.");
    }
    pa_ptr = reinterpret_cast< fastjet::PseudoJet * >(pa_void_ptr);
    pa = *pa_ptr;
    Py_XDECREF(py_result);  ///needed???
  }

private:
   PyObject *_py_class;  
};

//----------------------------------------------------------------------
// Since Python handles enum types as int, there can be some confusion
// between different JetDefinition ctors, where a int param (intended
// as a double, like using R=1 or p=-1 for the genkt algorithm) is
// actually interpreted an te enum (for the recombination scheme).
//
// We therefore provide a few helpers to force the construction of a
// Jet Definition with a fied number of parameters (+recombiner+strategy)
//
// JetDefinition0Param(algorithm, recomb_scheme, strategy)
JetDefinition JetDefinition0Param(JetAlgorithm jet_algorithm, 
                                  RecombinationScheme recomb_scheme = E_scheme,
                                  Strategy strategy = Best){
  return JetDefinition(jet_algorithm, recomb_scheme, strategy);
}

// JetDefinition1Param(algorithm, R, recomb_scheme, strategy)
JetDefinition JetDefinition1Param(JetAlgorithm jet_algorithm, 
                                  double R_in, 
                                  RecombinationScheme recomb_scheme = E_scheme,
                                  Strategy strategy = Best){
  return JetDefinition(jet_algorithm, R_in, recomb_scheme, strategy);
}

// JetDefinition2Param(algorithm, R, extrarecomb_scheme, strategy)
JetDefinition JetDefinition2Param(JetAlgorithm jet_algorithm, 
                                  double R_in, 
                                  double xtra_param,
                                  RecombinationScheme recomb_scheme = E_scheme,
                                  Strategy strategy = Best){
  return JetDefinition(jet_algorithm, R_in, xtra_param, recomb_scheme, strategy);
}


FASTJET_END_NAMESPACE      // defined in fastjet/internal/base.hh

#endif // __FASTJET_PYTHONUSERINFO_HH__
