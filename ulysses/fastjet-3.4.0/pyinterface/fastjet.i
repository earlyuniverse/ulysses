// -*-c++-*-
// for info on getting documentation, see https://github.com/m7thon/doxy2swig
//
%include "std_string.i"
%include "std_vector.i"
// read in documentation generated using doxygen and https://github.com/m7thon/doxy2swig
%include "fastjet-doc.i"


%define DOCSTRING
"Python interface to the FastJet jet clustering package.  

Usage is similar to the C++ case, with a few small changes noted
below.

Notes
-----

- You can pass a python list such as [PseudoJet0, PseudoJet1, ...]
  to any FastJet call that expects a vector of PseudoJets

- Any FastJet call that in C++ returns a vector of PseudoJets will in
  python return a tuple (or list) of PseudoJets

- for many objects that provide definitions of some kind, __str__
  call maps to description(). So, for example, you can just do

       jet_def = fastjet.JetDefinition(fastjet.antikt_algorithm, 0.4)
       print jet_def

- for combinations of selectors, (&&, || and !) in C++ map to 
  (&, | and ~) in python

- Selector::pass is remapped to Selector._pass

- remember that python uses reference, e.g. a = b means that a is a
  reference to b. If you need to copy a PseudoJet (pj), with a view to
  altering it, do 'pjcopy = PseudoJet(pj)'

- the python documentation has been automatically generated from the
  C++ doxygen documentation: python/C++ differences are not indicated,
  and certain methods and classes may be documented that were not
  included in the python conversion and/or configured for this
  particular installation.

Example
-------

  from fastjet import *
  particles = []
  particles.append(PseudoJet(100.0, 0.0, 0.0, 100.0)) # px, py, pz, E
  particles.append(PseudoJet(150.0, 0.0, 0.0, 150.0))

  R = 0.4
  jet_def = JetDefinition(antikt_algorithm, R)

  jets = jet_def(particles)

  print jet_def
  for jet in jets: print jet

"
%enddef

%module(docstring=DOCSTRING) fastjet

%{
#include "fastjet/config_auto.h"
#include "fastjet/config.h"
#include "fastjet/internal/base.hh"
#include "fastjet/internal/numconsts.hh"
#include "fastjet/internal/IsBase.hh"
#include "fastjet/internal/deprecated.hh"
#include "fastjet/internal/BasicRandom.hh"
#include "fastjet/SharedPtr.hh"
#include "fastjet/LimitedWarning.hh"
  //#include "fastjet/Error.hh"
#include "fastjet/PseudoJetStructureBase.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/FunctionOfPseudoJet.hh"
#include "fastjet/RangeDefinition.hh"
#include "fastjet/Selector.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/CompositeJetStructure.hh"
#include "fastjet/WrappedStructure.hh"
#include "fastjet/ClusterSequenceStructure.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/RectangularGrid.hh"
#include "fastjet/NNBase.hh"
#include "fastjet/NNH.hh"
#include "fastjet/NNFJN2Plain.hh"
#include "fastjet/NNFJN2Tiled.hh"
#include "fastjet/GhostedAreaSpec.hh"
#include "fastjet/AreaDefinition.hh"
#include "fastjet/ClusterSequenceAreaBase.hh"
#include "fastjet/ClusterSequenceActiveAreaExplicitGhosts.hh"
#include "fastjet/ClusterSequenceActiveArea.hh"
#include "fastjet/ClusterSequence1GhostPassiveArea.hh"
#include "fastjet/ClusterSequencePassiveArea.hh"
#include "fastjet/ClusterSequenceVoronoiArea.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "FastJetPythonExtensions.hh"
  %}


%template(vectorPJ) std::vector<fastjet::PseudoJet>;

//----------------------------------------------------------------------
// handle exceptoins by introducing a FastJetError python class
// derived from python's Exception class
//
// Thanks to Patrick Komiske for the fix
//----------------------------------------------------------------------

// define as macro for use in contrib files
%define FASTJET_ERRORS_AS_PYTHON_EXCEPTIONS(module)
%{
// Python class for representing errors from FastJet
static PyObject * FastJetError_;
%}

// this gets placed in the SWIG_init function
%init %{
  fastjet::Error::set_print_errors(false);
  unsigned int mlen = strlen(`module`);
  char * msg = (char*) calloc(mlen+15, sizeof(char));
  strcpy(msg, `module`);
  strcat(msg, ".FastJetError");
  FastJetError_ = PyErr_NewException(msg, NULL, NULL);
  Py_INCREF(FastJetError_);
  if (PyModule_AddObject(m, "FastJetError", FastJetError_) < 0) {
    Py_DECREF(m);
    Py_DECREF(FastJetError_);
    //return NULL;
  }
%}
%enddef

// include FastJetError in python module
%pythoncode {
  from _fastjet import FastJetError
}

FASTJET_ERRORS_AS_PYTHON_EXCEPTIONS(fastjet)

// check for fastjet errors after each action
%exception {
  try { $action }
  catch (fastjet::Error & e) {
    PyErr_SetString(FastJetError_, e.message().c_str());
    SWIG_fail;
  }
}


// // old code
// //
// // code to ensure that FJ C++ exceptions are passed on to Python
// // cf https://stackoverflow.com/questions/15006048/dynamically-rethrowing-self-defined-c-exceptions-as-python-exceptions-using-sw
// %exception {
//   try { $action }
//   catch (fastjet::Error &_e) {
//     SWIG_Python_Raise(SWIG_NewPointerObj(
//         (new fastjet::Error(static_cast<const fastjet::Error& >(_e))),  
//             SWIGTYPE_p_fastjet__Error,SWIG_POINTER_OWN),
//         "fastjet::Error", SWIGTYPE_p_fastjet__Error); 
//     SWIG_fail;
//   } 
// }


//----------------------------------------------------------------------
// a list of FastJet includes
//----------------------------------------------------------------------

%include "fastjet/config_auto.h"
%include "fastjet/config.h"
%include "fastjet/internal/base.hh"
%include "fastjet/internal/numconsts.hh"
 //%include "fastjet/internal/IsBase.hh"
%include "fastjet/internal/deprecated.hh"
%include "fastjet/internal/BasicRandom.hh"
%include "fastjet/SharedPtr.hh"
%include "fastjet/LimitedWarning.hh"
%include "fastjet/Error.hh"
%include "fastjet/PseudoJetStructureBase.hh"
%include "fastjet/PseudoJet.hh"
%include "fastjet/FunctionOfPseudoJet.hh"
%include "fastjet/RangeDefinition.hh"
%include "fastjet/Selector.hh"
%include "fastjet/JetDefinition.hh"
%include "fastjet/CompositeJetStructure.hh"
%include "fastjet/ClusterSequenceStructure.hh"
%include "fastjet/ClusterSequence.hh"
%include "fastjet/RectangularGrid.hh"
%include "fastjet/NNBase.hh"
%include "fastjet/NNH.hh"
%include "fastjet/NNFJN2Plain.hh"
%include "fastjet/NNFJN2Tiled.hh"
%include "fastjet/GhostedAreaSpec.hh"
%include "fastjet/AreaDefinition.hh"
%include "fastjet/ClusterSequenceAreaBase.hh"
%include "fastjet/ClusterSequenceActiveAreaExplicitGhosts.hh"
%include "fastjet/ClusterSequenceActiveArea.hh"
%include "fastjet/ClusterSequence1GhostPassiveArea.hh"
%include "fastjet/ClusterSequencePassiveArea.hh"
%include "fastjet/ClusterSequenceVoronoiArea.hh"
%include "fastjet/ClusterSequenceArea.hh"
%include "FastJetPythonExtensions.hh"

//----------------------------------------------------------------------
// extra bits and pieces to make a series of specialised things
// directly available in Python
//----------------------------------------------------------------------

namespace fastjet {

// the templated ctors must be specialised for PseudoJet
%define FASTJET_TEMPLATED_CTOR_FOR_PSEUDOJET(Class)
%extend Class {
   %template(Class) Class<PseudoJet>;
}
%enddef

// a macro to get support for description through __str__ method
%define FASTJET_SWIG_ADD_STR(Class)
%extend Class {
  std::string __str__() const {return $self->description();}
}
%enddef


// class descriptions (__str__)
FASTJET_SWIG_ADD_STR(JetDefinition)  
FASTJET_SWIG_ADD_STR(AreaDefinition)  
FASTJET_SWIG_ADD_STR(Selector)
FASTJET_SWIG_ADD_STR(GhostedAreaSpec)
FASTJET_SWIG_ADD_STR(FunctionOfPseudoJet)
FASTJET_SWIG_ADD_STR(RectangularGrid)

// templated ctors
FASTJET_TEMPLATED_CTOR_FOR_PSEUDOJET(ClusterSequence)
FASTJET_TEMPLATED_CTOR_FOR_PSEUDOJET(ClusterSequenceActiveAreaExplicitGhosts)
FASTJET_TEMPLATED_CTOR_FOR_PSEUDOJET(ClusterSequenceActiveArea)
FASTJET_TEMPLATED_CTOR_FOR_PSEUDOJET(ClusterSequence1GhostPassiveArea)
FASTJET_TEMPLATED_CTOR_FOR_PSEUDOJET(ClusterSequencePassiveArea)
FASTJET_TEMPLATED_CTOR_FOR_PSEUDOJET(ClusterSequenceVoronoiArea)
FASTJET_TEMPLATED_CTOR_FOR_PSEUDOJET(ClusterSequenceArea)

// These make JetDefinition, Selector and PseudoJet all printable
%extend JetDefinition {
  std::vector<PseudoJet> __call__(const std::vector<PseudoJet> & particles) {
    return (*self)(particles);
  }

  void set_python_recombiner(PyObject * pyobj){
    fastjet::RecombinerPython *new_python_recombiner = new fastjet::RecombinerPython(pyobj);
    $self->set_recombiner(new_python_recombiner);
    $self->delete_recombiner_when_unused();
  }
}

%extend Selector {
  std::string __str__() {return $self->description();}

  // The C++ operators [* && || !] map to [* & | ~] in python
  Selector __mul__   (const Selector & other) {return *($self) *  other;}
  Selector __and__   (const Selector & other) {return *($self) && other;}
  Selector __or__    (const Selector & other) {return *($self) || other;}
  Selector __invert__()                       {return !(*($self));}
 }

%extend PseudoJet {
  //PseudoJet(const PseudoJet & p) {return new fastjet::PseudoJet(p);}

  std::string __str__() {
    const unsigned int len_max=4096;
    char temp[len_max];
    snprintf(temp,len_max, "[%f, %f, %f, %f]",$self->px(), $self->py(), $self->pz(), $self->E());
    return std::string(temp);
  }

  void set_python_info(PyObject * pyobj) {
    fastjet::UserInfoPython * new_python_info = new fastjet::UserInfoPython(pyobj);
    $self->set_user_info(new_python_info);
  }

  PyObject * python_info() const {
    return $self->user_info<fastjet::UserInfoPython>().get_pyobj();
  }
  
  // these C++ operators are not automatically handled by SWIG (would only
  // be handled if there were part of the class)
  PseudoJet __add__ (const PseudoJet & p) {return *($self) + p;}
  PseudoJet __sub__ (const PseudoJet & p) {return *($self) - p;}
  bool      __eq__  (const PseudoJet & p) {return *($self) == p;}
  bool      __ne__  (const PseudoJet & p) {return *($self) != p;}
  PseudoJet __mul__ (double x) {return *($self) * x;}
  PseudoJet __rmul__(double x) {return *($self) * x;}
  PseudoJet __div__ (double x) {return *($self) / x;}
  bool      __eq__  (double x) {return *($self) == x;}
  bool      __ne__  (double x) {return *($self) != x;}
}

%template(FunctionOfPseudoJetDouble) FunctionOfPseudoJet<double>;
%template(FunctionOfPseudoJetPseudoJet) FunctionOfPseudoJet<fastjet::PseudoJet>;
  
}

%{
#include "fastjet/tools/Transformer.hh"
#include "fastjet/tools/Boost.hh"
#include "fastjet/tools/Recluster.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/Pruner.hh"
#include "fastjet/tools/CASubJetTagger.hh"
#include "fastjet/tools/MassDropTagger.hh"
#include "fastjet/tools/RestFrameNSubjettinessTagger.hh"
#include "fastjet/tools/TopTaggerBase.hh"
#include "fastjet/tools/JHTopTagger.hh"
#include "fastjet/tools/BackgroundEstimatorBase.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"
#include "fastjet/tools/GridMedianBackgroundEstimator.hh"
#include "fastjet/tools/Subtractor.hh"
%}


%include "fastjet/tools/Transformer.hh"
%include "fastjet/tools/Boost.hh"
%include "fastjet/tools/Recluster.hh"
%include "fastjet/tools/Filter.hh"
%include "fastjet/tools/Pruner.hh"
%include "fastjet/tools/CASubJetTagger.hh"
%include "fastjet/tools/MassDropTagger.hh"
%include "fastjet/tools/RestFrameNSubjettinessTagger.hh"
%include "fastjet/tools/TopTaggerBase.hh"
%include "fastjet/tools/JHTopTagger.hh"
%include "fastjet/tools/BackgroundEstimatorBase.hh"
%include "fastjet/tools/JetMedianBackgroundEstimator.hh"
%include "fastjet/tools/GridMedianBackgroundEstimator.hh"
%include "fastjet/tools/Subtractor.hh"


namespace fastjet{
// the access to the tools operator() is automatically handled by SWIG
//
// Here, we add the __str__ suppport so the description can be printed
// easily
%define FASTJET_OPERATOR_PARENTHESIS_CALLABLE(Class)
%extend Class {
  std::string  __str__() {
    return $self->description();
  }
}
%enddef
  
// class descriptions (__str__) for the tools
FASTJET_SWIG_ADD_STR(Boost)
FASTJET_SWIG_ADD_STR(Unboost)
FASTJET_SWIG_ADD_STR(Recluster)
FASTJET_SWIG_ADD_STR(Filter)
FASTJET_SWIG_ADD_STR(Pruner)
FASTJET_SWIG_ADD_STR(CASubJetTagger)
FASTJET_SWIG_ADD_STR(MassDropTagger)
FASTJET_SWIG_ADD_STR(RestFrameNSubjettinessTagger)
FASTJET_SWIG_ADD_STR(JHTopTagger)
FASTJET_SWIG_ADD_STR(Subtractor)
FASTJET_SWIG_ADD_STR(Error)

    
}
