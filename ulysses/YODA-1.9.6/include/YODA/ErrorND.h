// -*- C++ -*-
//
// This file is part of YODA -- Yet more Objects for Data Analysis
// Copyright (C) 2008-2021 The YODA collaboration (see AUTHORS for details)
//
#ifndef YODA_ERRORND_H
#define YODA_ERRORND_H

#include "YODA/Exceptions.h"
#include "YODA/Utils/MathUtils.h"
#include "YODA/Utils/sortedvector.h"
#include "YODA/Utils/ndarray.h"
#include <utility>
#include <algorithm>

namespace YODA {


  /// An N-dimensional error to be contained in a Point<N>
  template<int N>
  class Error {
  public:

    typedef Utils::ndarray<double, N> NdVal;
    typedef Utils::ndarray<std::pair<double,double>, N> NdValPair;


    /// @name Constructors
    /// @{

    // Default constructor
    Error(const std::string& name="")
      : _name(name)
    {
      clear();
    }


    /// Constructor from ND array of error pairs
    Error(const NdValPair& err, const std::string& name="")
      : _name(name), _val(err)
    {    }


    /// Constructor of a symmetric error from one ND array of values
    Error(const NdVal& errsymm, const std::string& name="")
      : _name(name)
    {
      for (size_t i = 0; i < N; ++i) {
        _val[i] = std::make_pair(errsymm[i], errsymm[i]);
      }
    }


    /// Constructor of an asymmetric error from two ND arrays of values
    Error(const NdVal& errminus, const NdVal& errplus, const std::string& name="")
      : _name(name)
    {
      for (size_t i = 0; i < N; ++i) {
        _val[i] = std::make_pair(errminus[i], errplus[i]);
      }
    }

    /// @}


    /// @name Modifiers and accessors
    /// @{

    /// Clear the point values and errors
    const std::string& name() const {
      return _name;
    }

    /// Clear the point values and errors
    void setName(const std::string& name) {
      _name = name;
    }

    /// Clear the point values and errors
    void clear() {
      _val.clear();
    }

    // /// Get the error pair array
    // NdValPair& errs() {
    //   return _val;
    // }

    // /// Get the error pair array (const version)
    // const NdValPair& errs() const {
    //   return _val;
    // }

    /// Access the error pair in dimension @a dim
    std::pair<double,double>& err(size_t dim) {
      return _val[dim];
    }

    /// Get the error pair in dimension @a dim (const version)
    const std::pair<double,double>& err(size_t dim) const {
      return _val[dim];
    }

    /// Access the error pair in dimension @a dim
    std::pair<double,double>& operator[](size_t dim) {
      return _val[dim];
    }

    /// Get the error pair in dimension @a dim (const version)
    const std::pair<double,double>& operator[](size_t dim) const {
      return _val[dim];
    }

    /// Get the minus error in dimension @a dim
    double errMinus(size_t dim) const {
      return _val[dim].first;
    }

    /// Get the plus error in dimension @a dim
    double errPlus(size_t dim) const {
      return _val[dim].second;
    }

    /// Get the mean error in dimension @a dim
    double errAvg(size_t dim) const {
      return (_val[dim].first + _val[dim].second)/2.0;
    }

    /// @}



    /// @name Scaling and transformations
    /// @{

    /// Uniform scaling
    void scale(const NdVal& scales) {
      for (size_t i = 0; i < N; ++i) {
        _val[i].first *= scales[i];
        _val[i].second *= scales[i];
      }
    }

    /// @todo Generic trf functor support -- need abs position of error bar

    /// @}


  protected:

    /// @name Value and error variables
    /// @{

    std::string _name;
    NdValPair _val;

    /// @}

  };



  /// @name Comparison operators
  /// @{

  /// Equality test
  template <int N>
  inline bool operator==(const Error<N>& a, const Error<N>& b) {
    if (a.name() != b.name()) return false;
    for (size_t i = 0; i < N; ++i) {
      if (!fuzzyEquals(a.errMinus(i), b.errMinus(i))) return false;
      if (!fuzzyEquals(a.errPlus(i), b.errPlus(i))) return false;
    }
    return true;
  }

  /// Inequality test
  template <int N>
  inline bool operator!=(const Error<N>& a, const Error<N>& b) {
    return !(a == b);
  }


  /// Less-than operator used to sort errors
  template <int N>
  inline bool operator<(const Error<N>& a, const Error<N>& b) {
    #define LT_IF_NOT_EQ(a,b) { if (!fuzzyEquals(a, b)) return a < b; }
    for (size_t i = 0; i < N; ++i) {
      LT_IF_NOT_EQ(a.err(i).first, b.err(i).first);
      LT_IF_NOT_EQ(a.err(i).second, b.err(i).second);
    }
    #undef LT_IF_NOT_EQ
    return a.name() < b.name();
  }

  /// Less-than-or-equals operator used to sort errors
  template <int N>
  inline bool operator<=(const Error<N>& a, const Error<N>& b) {
    if (a == b) return true;
    return a < b;
  }

  /// Greater-than operator used to sort errors
  template <int N>
  inline bool operator>(const Error<N>& a, const Error<N>& b) {
    return !(a <= b);
  }

  /// Greater-than-or-equals operator used to sort errors
  template <int N>
  inline bool operator>=(const Error<N>& a, const Error<N>& b) {
    return !(a < b);
  }

  /// @}


}

#endif
