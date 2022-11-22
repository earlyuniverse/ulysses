// -*- C++ -*-
//
// This file is part of YODA -- Yet more Objects for Data Analysis
// Copyright (C) 2008-2021 The YODA collaboration (see AUTHORS for details)
//
#ifndef YODA_Counter_h
#define YODA_Counter_h

#include "YODA/AnalysisObject.h"
#include "YODA/Fillable.h"
#include "YODA/Dbn0D.h"
#include "YODA/Scatter1D.h"
#include "YODA/Exceptions.h"
#include <vector>
#include <string>
#include <map>
#include <tuple>
#include <memory>

namespace YODA {


  /// A weighted counter.
  class Counter : public AnalysisObject, public Fillable {
  public:


    typedef std::tuple<> FillType;
    typedef std::shared_ptr<Counter> Ptr;

    /// @name Constructors
    /// @{

    /// Default constructor
    Counter(const std::string& path="", const std::string& title="")
      : AnalysisObject("Counter", path, title)
    { }


    /// @brief Constructor accepting an explicit Dbn0D.
    ///
    /// Intended both for internal persistency and user use.
    Counter(const Dbn0D& dbn,
            const std::string& path="", const std::string& title="")
            : AnalysisObject("Counter", path, title),
            _dbn(dbn)
    { }


    /// @brief Constructor accepting a double (treated as the weight of a single fill).
    ///
    /// Intended for user convenience only, so Counter can be treated as a number.
    Counter(double w,
            const std::string& path="", const std::string& title="")
            : AnalysisObject("Counter", path, title)
    {
      _dbn.fill(w);
    }


    /// Copy constructor with optional new path
    /// @todo Don't copy the path?
    Counter(const Counter& c, const std::string& path="");


    /// Assignment operator
    Counter& operator = (const Counter& c) {
      setPath(c.path());
      setTitle(c.title());
      _dbn = c._dbn;
      return *this;
    }

    /// Make a copy on the stack
    Counter clone() const {
      return Counter(*this);
    }

    /// Make a copy on the heap, via 'new'
    Counter* newclone() const {
      return new Counter(*this);
    }

    /// @}


    /// @name Dimensions
    /// @{

    /// @brief Fill dimension of this data object
    ///
    /// @todo Change to return total dimension
    size_t dim() const { return 0; }

    /// Fill dimension of this data object
    size_t fillDim() const { return 0; }

    /// @}


    /// @name Modifiers
    /// @{

    /// Fill histo by value and weight
    virtual void fill(double weight=1.0, double fraction=1.0) {
      _dbn.fill(weight, fraction);
    }

    virtual void fill(FillType, double weight=1.0, double fraction=1.0) {
      fill(weight, fraction);
    }

    /// @brief Reset the histogram.
    ///
    /// Keep the binning but set all bin contents and related quantities to zero
    virtual void reset() {
      _dbn.reset();
    }


    /// Rescale as if all fill weights had been different by factor @a scalefactor.
    void scaleW(double scalefactor) {
      setAnnotation("ScaledBy", annotation<double>("ScaledBy", 1.0) * scalefactor);
      _dbn.scaleW(scalefactor);
    }

    /// @}


    /// @name Data access
    /// @{

    /// Get the number of fills
    double numEntries(bool=false) const { return _dbn.numEntries(); }

    /// Get the effective number of fills
    double effNumEntries(bool=false) const { return _dbn.effNumEntries(); }

    /// Get the sum of weights
    double sumW(bool=false) const { return _dbn.sumW(); }

    /// Get the sum of squared weights
    double sumW2(bool=false) const { return _dbn.sumW2(); }

    /// Get the value
    double val(bool=false) const { return sumW(); }

    /// Get the uncertainty on the value
    /// @todo Implement on Dbn0D and feed through to this and Dbn1D, 2D, etc.
    // double err() const { return _dbn.err(); }
    double err() const { return sqrt(sumW2()); }

    /// Get the relative uncertainty on the value
    /// @todo Implement on Dbn0D and feed through to this and Dbn1D, 2D, etc.
    // double err() const { return _dbn.err(); }
    double relErr() const {
      /// @todo Throw excp if sumW2 is 0?
      return sumW2() != 0 ? err()/sumW() : 0;
    }

    /// @}


    /// @name Internal state access and modification (mainly for persistency use)
    /// @{

    /// Get the internal distribution object
    const Dbn0D& dbn() const {
      return _dbn;
    }

    /// Set the internal distribution object: CAREFUL!
    void setDbn(const Dbn0D& dbn) {
      _dbn = dbn;
    }

    // /// Set the whole object state
    // void setState(const Dbn0D& dbn, const AnalysisObject::Annotations& anns=AnalysisObject::Annotations()) {
    //   setDbn(dbn);
    //   setAnnotations(anns);
    // }

    /// @}


    /// @name Adding and subtracting counters
    /// @{

    /// Add another counter to this
    Counter& operator += (const Counter& toAdd) {
      _dbn += toAdd._dbn;
      return *this;
    }

    /// Subtract another counter from this
    Counter& operator -= (const Counter& toSubtract) {
      _dbn -= toSubtract._dbn;
      return *this;
    }

    /// Increment as if by a fill of weight = 1
    /// @note This is post-increment only, i.e. cn++ not ++cn
    Counter& operator ++ () {
      *this += 1;
      return *this;
    }

    /// Increment as if by a fill of weight = -1
    /// @note This is post-decrement only, i.e. cn-- not --cn
    Counter& operator -- () {
      *this -= 1;
      return *this;
    }

    /// Scale by a double (syntactic sugar for @c scaleW(s))
    Counter& operator *= (double s) {
      scaleW(s);
      return *this;
    }

    /// Inverse-scale by a double (syntactic sugar for @c scaleW(1/s))
    Counter& operator /= (double s) {
      scaleW(1/s);
      return *this;
    }

    /// @}


  private:

    /// @name Data
    /// @{

    /// Contained 0D distribution
    Dbn0D _dbn;

    /// @}

  };


  /// @name Combining counters: global operators
  /// @{

  /// Add two counters
  inline Counter add(const Counter& first, const Counter& second) {
    Counter tmp = first;
    tmp += second;
    return tmp;
  }

  /// Add two counters
  inline Counter operator + (const Counter& first, const Counter& second) {
    return add(first, second);
  }


  /// Subtract two counters
  inline Counter subtract(const Counter& first, const Counter& second) {
    Counter tmp = first;
    tmp -= second;
    return tmp;
  }

  /// Subtract two counters
  inline Counter operator - (const Counter& first, const Counter& second) {
    return subtract(first, second);
  }


  /// Divide two counters, with an uncorrelated error treatment
  /// @todo Or just return a Point1D?
  Scatter1D divide(const Counter& numer, const Counter& denom);

  /// Divide two counters, with an uncorrelated error treatment
  /// @todo Or just return a Point1D?
  inline Scatter1D operator / (const Counter& numer, const Counter& denom) {
    return divide(numer, denom);
  }

  /// @todo Add divide functions/operators on pointers



  /// @brief Calculate an efficiency ratio of two counters
  ///
  /// Note that an efficiency is not the same thing as a standard division of two
  /// histograms: the errors must be treated as correlated
  ///
  /// @todo Or just return a Point1D?
  Scatter1D efficiency(const Counter& accepted, const Counter& total);

  /// @}


}

#endif
