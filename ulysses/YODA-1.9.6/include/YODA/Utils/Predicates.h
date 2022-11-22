// -*- C++ -*-
//
// This file is part of YODA -- Yet more Objects for Data Analysis
// Copyright (C) 2008-2017 The YODA collaboration (see AUTHORS for details)
//
#ifndef YODA_PREDICATES_H
#define YODA_PREDICATES_H

#include <iostream>

namespace YODA {


  /// @brief Functor to compare two floating point numbers and return whether they are fuzzily equivalent
  ///
  /// The equivalence comparison is of the form dev = |b-a|/refscale. The @a
  /// refscale argument may used to force the reference scale in this, otherwise
  /// (or if explicitly set to zero) it will default to (|a| + |b|)/2. The
  /// returned boolean is determined by the comparison |dev| < tol.
  struct CmpFloats {
    CmpFloats(double tol=1e-3, double refscale=0.0) : _tol(tol), _refscale(refscale) {}
    bool operator () (const double& a, const double& b) {
      const double div = (_refscale == 0) ? 0.5*(std::abs(a)+std::abs(b)) : _refscale;
      const double dev = (b-a)/div;
      // std::cout << "CmpFloats: " << a << " vs. " << b << " -> dev = " << dev << std::endl;
      return std::abs(dev) < _tol;
    }
    double _tol, _refscale;
  };


}

#endif
