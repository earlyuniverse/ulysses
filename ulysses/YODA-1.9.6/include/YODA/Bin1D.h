// -*- C++ -*-
//
// This file is part of YODA -- Yet more Objects for Data Analysis
// Copyright (C) 2008-2021 The YODA collaboration (see AUTHORS for details)
//
#ifndef YODA_Bin1D_h
#define YODA_Bin1D_h

#include "YODA/Utils/MathUtils.h"
#include "YODA/Bin.h"
#include <utility>

namespace YODA {


  /// @brief A generic 1D bin type
  ///
  /// This is a generic 1D bin type which supplies the accessors for the two "x"
  /// and "y" axis directions in which it is defined. Bin1D is not intended to be
  /// directly instantiated: it is inherited from to make specific histogram and
  /// profile bin types as HistoBin1D and ProfileBin1D.
  /// The lower bin edge is inclusive. This base class provides no fill
  /// method, since the signatures for standard and profile histos differ.
  ///
  /// @todo It would also be nice to have an *untemplated* generic Bin1D interface
  template <class DBN>
  class Bin1D : public Bin {
  public:

    /// @name Constructors
    /// @{

    // /// Make a new, empty bin with a pair of edges.
    // Bin1D(double lowedge, double highedge)
    //   : _edges( std::make_pair(lowedge, highedge) )
    // {
    //   if (_edges.second < _edges.first) {
    //     throw RangeError("The bin edges are wrongly defined!");
    //   }
    // }


    /// Make a new, empty bin with a pair of edges.
    Bin1D(const std::pair<double,double>& edges)
      : _edges(edges)
    {
      if (_edges.second < _edges.first) {
        throw RangeError("The bin edges are wrongly defined!");
      }
    }


    /// @brief Make a bin with all the components of a fill history.
    ///
    /// Mainly intended for internal persistency use.
    Bin1D(const std::pair<double,double>& edges, const DBN& dbn)
      : _edges(edges), _dbn(dbn)
    {
      if (_edges.second < _edges.first) {
        throw RangeError("The bin edges are wrongly defined!");
      }
    }


    /// Copy constructor
    Bin1D(const Bin1D<DBN>& b)
      : _edges(b._edges),
        _dbn(b._dbn)
    { }


    /// Copy assignment
    Bin1D& operator = (const Bin1D<DBN>& b) {
      _edges = b._edges;
      _dbn = b._dbn;
      return *this;
    }

    /// @}


    /// @name Dimensions
    /// @{

    /// Dimension of the fill space
    ///
    /// @todo Convert to total dimension
    size_t dim() { return 1; }

    /// Dimension of the fill space
    size_t fillDim() { return 1; }

    /// @}


    /// @name Modifiers
    /// @{

    /// Reset this bin
    virtual void reset() {
      _dbn.reset();
    }

    /// Rescale as if all fill weights had been different by factor @a scalefactor
    ///
    /// @note This should not be used, since it breaks histogram consistency. It will be removed in a future version.
    void scaleW(double scalefactor) {
      _dbn.scaleW(scalefactor);
    }

    /// Scale the x dimension
    ///
    /// @note This should not be used, since it breaks histogram consistency. It will be removed in a future version.
    void scaleX(double factor) {
      _edges.first *= factor;
      _edges.second *= factor;
      _dbn.scaleX(factor);
    }

    /// @}


  public:

    /// @name X-axis info
    /// @{

    /// Get the {low,high} edges as an STL @c pair.
    std::pair<double,double> xEdges() const {
      return _edges;
    }

    /// Lower limit of the bin (inclusive).
    double xMin() const {
      return _edges.first;
    }

    /// Upper limit of the bin (exclusive).
    double xMax() const {
      return _edges.second;
    }

    /// Geometric centre of the bin, i.e. high+low/2.0
    double xMid() const {
      return ( _edges.second + _edges.first ) / 2;
    }
    /// Alias for xMid
    /// @deprecated Only retained for temporary backward compatibility: use xMid
    double midpoint() const { return xMid(); }

    /// Separation of low and high edges, i.e. high-low.
    double xWidth() const {
      return _edges.second - _edges.first;
    }
    /// Alias for xWidth
    /// @deprecated Only retained for temporary backward compatibility: use xWidth
    double width() const { return xWidth(); }


    /// The mean position in the bin, or the midpoint if that is not available.
    double xFocus() const {
      return (!isZero(sumW())) ? xMean() : xMid();
    }

    /// @}


    /// @name X distribution statistics
    /// @{

    /// Mean value of x-values in the bin.
    double xMean() const {
      return _dbn.xMean();
    }

    /// The variance of x-values in the bin.
    double xVariance() const {
      return _dbn.xVariance();
    }

    /// The standard deviation (spread) of x-values in the bin.
    double xStdDev() const {
      return _dbn.xStdDev();
    }

    /// The standard error on the bin focus.
    double xStdErr() const {
      return _dbn.xStdErr();
    }

    /// The x RMS in the bin.
    double xRMS() const {
      return _dbn.xRMS();
    }

    /// @}


  public:

    /// @name Raw distribution statistics
    /// @{

    /// Statistical distribution in this bin (non-const)
    DBN& dbn() {
      return _dbn;
    }

    /// Statistical distribution in this bin (const)
    const DBN& dbn() const {
      return _dbn;
    }


    /// The number of entries
    double numEntries() const {
      return _dbn.numEntries();
    }

    /// The effective number of entries
    double effNumEntries() const {
      return _dbn.effNumEntries();
    }

    /// The sum of weights
    double sumW() const {
      return _dbn.sumW();
    }

    /// The sum of weights squared
    double sumW2() const {
      return _dbn.sumW2();
    }

    /// The sum of x*weight
    double sumWX() const {
      return _dbn.sumWX();
    }

    /// The sum of x^2 * weight
    double sumWX2() const {
      return _dbn.sumWX2();
    }

    /// @}


  public:

    /// @name Operators
    /// @{

    /// Add two bins
    Bin1D<DBN>& operator += (const Bin1D<DBN>& b) {
      return add(b);
    }

    /// Subtract one bin from another
    Bin1D<DBN>& operator -= (const Bin1D<DBN>& b) {
      return subtract(b);
    }

    /// @}


    /// @name Named operators
    /// @{

    /// Merge two adjacent bins
    Bin1D<DBN>& merge(const Bin1D<DBN>& b) {
      if (fuzzyEquals(_edges.second, b._edges.first)) {
        _edges.second = b._edges.second;
      } else if (fuzzyEquals(_edges.second, b._edges.first)) {
        _edges.first = b._edges.first;
      } else {
        throw LogicError("Attempted to merge two non-adjacent bins");
      }
      // std::cout << "a " << _dbn.sumW() << std::endl;
      _dbn += b._dbn;
      // std::cout << "b " << _dbn.sumW() << std::endl;
      return *this;
    }


    /// Add two bins (internal, explicitly named version)
    ///
    /// This operator is defined for adding two bins with equivalent binning.
    /// It cannot be used to merge two bins into one larger bin.
    Bin1D<DBN>& add(const Bin1D<DBN>& b) {
      if (!fuzzyEquals(_edges.first, b._edges.first) ||
          !fuzzyEquals(_edges.second, b._edges.second)) {
        throw LogicError("Attempted to add two bins with different edges");
      }
      _dbn += b._dbn;
      return *this;
    }


    /// Subtract one bin from another (internal, explicitly named version)
    ///
    /// This operator is defined for subtracting two bins with equivalent binning.
    /// It cannot be used to merge two bins into one larger bin.
    Bin1D<DBN>& subtract(const Bin1D<DBN>& b) {
      if (!fuzzyEquals(_edges.first, b._edges.first) ||
          !fuzzyEquals(_edges.second, b._edges.second)) {
        throw LogicError("Attempted to subtract two bins with different edges");
      }
      _dbn -= b._dbn;
      return *this;
    }

    /// @}


  protected:

    /// The bin limits
    std::pair<double,double> _edges;

    // Distribution of weighted x (and perhaps y) values
    DBN _dbn;

    /// @todo Make Axis1D<BIN, DBN> -> Axis1D<DBN> and hold a pointer to the one that contains this bin

  };



  /// Add two bins
  ///
  /// This "add" operator is defined for adding two bins with equivalent binning.
  /// It cannot be used to merge two bins into one larger bin.
  template <class DBN>
  inline Bin1D<DBN> operator + (const Bin1D<DBN>& a, const Bin1D<DBN>& b) {
    Bin1D<DBN> rtn = a;
    rtn += b;
    return rtn;
  }


  /// Subtract one bin from another
  ///
  /// This "subtraction" operator is defined for subtracting two bins with equivalent binning.
  /// It cannot be used to merge two bins into one larger bin.
  template <class DBN>
  inline Bin1D<DBN> operator - (const Bin1D<DBN>& a, const Bin1D<DBN>& b) {
    Bin1D<DBN> rtn = a;
    rtn -= b;
    return rtn;
  }


  /// Bin1Ds are compared for axis sorting by lower edge position
  template <class DBN>
  inline bool operator<(const Bin1D<DBN>& a, const Bin1D<DBN>& b) {
    return b.xEdges().first > a.xEdges().first;
  }


}



#endif
