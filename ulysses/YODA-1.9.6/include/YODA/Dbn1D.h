// -*- C++ -*-
//
// This file is part of YODA -- Yet more Objects for Data Analysis
// Copyright (C) 2008-2021 The YODA collaboration (see AUTHORS for details)
//
#ifndef YODA_Dbn1D_h
#define YODA_Dbn1D_h

#include "YODA/Dbn0D.h"
#include "YODA/Exceptions.h"
#include "YODA/Utils/MathUtils.h"
#include <cmath>

namespace YODA {


  /// @brief A 1D distribution
  ///
  /// This class is used internally by YODA to centralise the calculation of
  /// statistics of unbounded, unbinned sampled distributions. Each distribution
  /// fill contributes a weight, \f$ w \f$, and a value, \f$ x \f$. By storing
  /// the total number of fills (ignoring weights), \f$ \sum w \f$, \f$ \sum w^2
  /// \f$, \f$ \sum wx \f$, and \f$ \sum wx^2 \f$, the Dbn1D can calculate the
  /// mean and spread (\f$ \sigma^2 \f$, \f$ \sigma \f$ and \f$ \hat{\sigma}
  /// \f$) of the sampled distribution. It is used to provide this information
  /// in bins and for the "hidden" \f$ y \f$ distribution in profile histogram
  /// bins.
  class Dbn1D {
  public:

    /// @name Constructors
    /// @{

    /// Default constructor of a new distribution.
    Dbn1D() {
      reset();
    }


    /// @brief Constructor to set a distribution with a pre-filled state.
    ///
    /// Principally designed for internal persistency use.
    Dbn1D(double numEntries, double sumW, double sumW2, double sumWX, double sumWX2)
      : _dbnW(numEntries, sumW, sumW2),
        _sumWX(sumWX),
        _sumWX2(sumWX2)
    { }


    /// Copy constructor
    ///
    /// Sets all the parameters using the ones provided from an existing Dbn1D.
    Dbn1D(const Dbn1D& toCopy) {
      _dbnW = toCopy._dbnW;
      _sumWX = toCopy._sumWX;
      _sumWX2 = toCopy._sumWX2;
    }


    /// Copy assignment
    ///
    /// Sets all the parameters using the ones provided from an existing Dbn1D.
    Dbn1D& operator=(const Dbn1D& toCopy) {
      _dbnW = toCopy._dbnW;
      _sumWX = toCopy._sumWX;
      _sumWX2 = toCopy._sumWX2;
      return *this;
    }

    /// @}


    /// @name Modifiers
    /// @{

    /// @brief Contribute a sample at @a val with weight @a weight.
    void fill(double val, double weight=1.0, double fraction=1.0) {
      _dbnW.fill(weight, fraction);
      _sumWX += fraction*weight*val;
      _sumWX2 += fraction*weight*val*val;
    }


    /// Reset the internal counters.
    void reset() {
      _dbnW.reset();
      _sumWX = 0;
      _sumWX2 = 0;
    }


    /// Rescale as if all fill weights had been different by factor @a scalefactor.
    void scaleW(double scalefactor) {
      _dbnW.scaleW(scalefactor);
      _sumWX *= scalefactor;
      _sumWX2 *= scalefactor;
    }


    /// Rescale x: needed if histo bin edges are rescaled.
    void scaleX(double factor) {
      _sumWX *= factor;
      _sumWX2 *= factor*factor;
    }

    /// @}

  public:


    /// @name Distribution statistics
    /// @{

    /// The absolute error on sumW
    double errW() const { return _dbnW.errW(); }

    /// The relative error on sumW
    double relErrW() const { return _dbnW.relErrW(); }

    /// Weighted mean, \f$ \bar{x} \f$, of distribution.
    double xMean() const;

    /// Weighted variance, \f$ \sigma^2 \f$, of distribution.
    double xVariance() const;

    /// Weighted standard deviation, \f$ \sigma \f$, of distribution.
    double xStdDev() const { return std::sqrt(xVariance()); }

    /// Weighted standard error on the mean, \f$ \sim \sigma/\sqrt{N-1} \f$, of distribution.
    double xStdErr() const;

    /// Weighted RMS, \f$ \sqrt{ \sum{w x^2}/\sum{w} } \f$, of distribution.
    double xRMS() const;

    /// @}


    /// @name Raw distribution running sums
    /// @{

    /// Number of entries (number of times @c fill was called, ignoring weights)
    double numEntries() const {
      return _dbnW.numEntries();
    }

    /// Effective number of entries \f$ = (\sum w)^2 / \sum w^2 \f$
    double effNumEntries() const {
      return _dbnW.effNumEntries();
    }

    /// The sum of weights
    double sumW() const {
      return _dbnW.sumW();
    }

    /// The sum of weights squared
    double sumW2() const {
      return _dbnW.sumW2();
    }

    /// The sum of x*weight
    double sumWX() const {
      return _sumWX;
    }

    /// The sum of x^2*weight
    double sumWX2() const {
      return _sumWX2;
    }

    /// @}


    /// @name Operators
    /// @{

    /// Add two dbns
    Dbn1D& operator += (const Dbn1D& d) {
      return add(d);
    }

    /// Subtract one dbn from another
    Dbn1D& operator -= (const Dbn1D& d) {
      return subtract(d);
    }

    /// @}


  protected:

    /// Add two dbns (internal, explicitly named version)
    Dbn1D& add(const Dbn1D& d);

    /// Subtract one dbn from another (internal, explicitly named version)
    Dbn1D& subtract(const Dbn1D& d);


  private:

    /// @name Storage
    /// @{

    /// The pure weight moments are stored in a 0D distribution
    Dbn0D _dbnW;

    /// The 1st order weighted x moment
    double _sumWX;

    /// The 2nd order weighted x moment
    double _sumWX2;

    /// @}

  };


  /// Add two dbns
  inline Dbn1D operator + (const Dbn1D& a, const Dbn1D& b) {
    Dbn1D rtn = a;
    rtn += b;
    return rtn;
  }

  /// Subtract one dbn from another
  inline Dbn1D operator - (const Dbn1D& a, const Dbn1D& b) {
    Dbn1D rtn = a;
    rtn -= b;
    return rtn;
  }


}

#endif
