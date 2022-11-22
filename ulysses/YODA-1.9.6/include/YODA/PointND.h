// -*- C++ -*-
//
// This file is part of YODA -- Yet more Objects for Data Analysis
// Copyright (C) 2008-2021 The YODA collaboration (see AUTHORS for details)
//
#ifndef YODA_POINTND_H
#define YODA_POINTND_H

#include "YODA/Exceptions.h"
#include "YODA/ErrorND.h"
#include "YODA/Utils/MathUtils.h"
#include "YODA/Utils/sortedvector.h"
#include "YODA/Utils/ndarray.h"
#include <utility>
#include <algorithm>

namespace YODA {


  /// An N-dimensional data point to be contained in a Scatter<N>
  template<int N>
  class Point {
  public:

    // Typedefs
    typedef Utils::ndarray<double, N> NdVal;
    typedef Utils::ndarray<std::pair<double,double>, N> NdValPair;
    typedef Utils::sortedvector< Error<N> > Errors;


    /// @name Constructors
    /// @{

    // Default constructor
    Point() {
      clear();
    }


    /// Constructor from position values without errors
    Point(const NdVal& pos)
      : _pos(pos)
    {    }


    // /// Constructor from values with a single set of symmetric errors
    // /// @todo Unnecessary since Error can be implicitly constructed this way
    // Point(const NdVal& pos, const NdVal& errs)
    //   : _pos(pos)
    // {
    //   _errs.insert(Error<N>(errs));
    // }


    /// Constructor from values with a single set of asymmetric errors
    Point(const NdVal& pos, const NdValPair& errs)
      : _pos(pos)
    {
      _errs.insert(Error<N>(errs));
    }


    /// Constructor from values and a single Error object
    Point(const NdVal& pos, const Error<N>& err)
      : _pos(pos)
    {
      _errs.insert(err);
    }


    /// Constructor from values and a collection of Error objects
    Point(const std::vector<double>& pos, const std::vector< Error<N> >& errs)
      : _pos(pos), _errs(errs)
    {    }

    /// @}


    /// @name Modifiers
    /// @{

    /// Clear the point values and errors
    void clear() {
      for (size_t i = 0; i < N; ++i) _pos[i] = 0;
      _errs.clear();
    }

    /// @todo addError, addErrors, setErrors

    /// @}


  public:

    /// @name Coordinate accessors
    /// @{

    /// Get the coordinate vector
    NdVal& pos() { return _pos; }

    /// Get the coordinate vector (const version)
    const NdVal& pos() const { return _pos; }

    /// Set the coordinate vector
    void setPos(const NdVal& pos) {
      _pos = pos;
    }

    /// @}


    /// @name Error accessors
    /// @{

    /// Get error values
    Errors& errs() {
      return _errs;
    }

    /// Get error values (const version)
    const Errors& errs() const {
      return _errs;
    }

    /// Set the error values
    void setErrs(const Errors& errs) {
      _errs = errs;
    }

    /// @}


    /// @name Scaling and transformations
    /// @{

    /// Uniform scaling
    void scale(const NdVal& scales) {
      for (size_t i = 0; i < N; ++i) _pos[i] *= scales[i];
      for (Error<N>& e : errs()) e.scale(scales);
    }


    // /// Generalised transformations with functors
    // void scale(const Trf<N>& trf) {
    //   for (size_t i = 0; i < N; ++i)
    //     _pos = trf.transform(_pos);
    //   for (Error e : errs())
    //     rf.transformErrs(_pos, e);
    // }

    /// @}


  protected:

    /// @name Value and error variables
    /// @{

    NdVal _pos;
    Errors _errs;

    /// @}

  };



  /// @name Comparison operators
  /// @{

  /// Equality test
  template <int N>
  inline bool operator==(const Point<N>& a, const Point<N>& b) {
    // Compare positions
    for (size_t i = 0; i < N; ++i) {
      if ( !fuzzyEquals(a.pos()[i], b.pos()[i]) ) return false;
    }
    // Compare number of errors and then (sorted) error equality
    if (a.errs().size() != b.errs().size()) return false;
    for (size_t i = 0; i < a.errs().size(); ++i) {
      if (a.errs()[i] != b.errs()[i]) return false;
    }
    return true;
  }

  /// Inequality test
  template <int N>
  inline bool operator!=(const Point<N>& a, const Point<N>& b) {
    return !(a == b);
  }


  /// Less-than operator used to sort points
  template <int N>
  inline bool operator<(const Point<N>& a, const Point<N>& b) {
    #define LT_IF_NOT_EQ(a,b) { if (!fuzzyEquals(a, b)) return a < b; }
    for (size_t i = 0; i < N; ++i) LT_IF_NOT_EQ(a.pos()[i], b.pos()[i]);
    if (a.errs().size() != b.errs().size()) return a.errs().size() < b.errs().size();
    for (size_t i = 0; i < a.errs().size(); ++i) {
      if (a.errs()[i] != b.errs()[i]) return a.errs()[i] < b.errs()[i];
    }
    #undef LT_IF_NOT_EQ
    return false;
  }

  /// Less-than-or-equals operator used to sort points
  template <int N>
  inline bool operator<=(const Point<N>& a, const Point<N>& b) {
    if (a == b) return true;
    return a < b;
  }

  /// Greater-than operator used to sort points
  template <int N>
  inline bool operator>(const Point<N>& a, const Point<N>& b) {
    return !(a <= b);
  }

  /// Greater-than-or-equals operator used to sort points
  template <int N>
  inline bool operator>=(const Point<N>& a, const Point<N>& b) {
    return !(a < b);
  }

  /// @}


}

#endif
