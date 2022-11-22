// -*- C++ -*-
//
// This file is part of YODA -- Yet more Objects for Data Analysis
// Copyright (C) 2008-2021 The YODA collaboration (see AUTHORS for details)
//
#ifndef YODA_Bin2D_h
#define YODA_Bin2D_h

#include "YODA/Utils/MathUtils.h"
#include "YODA/Bin.h"
#include <utility>

namespace YODA {


  /// @brief A generic 2D bin type
  ///
  /// This is a generic 2D bin type which supplies the accessors for the two "x"
  /// and "y" axis directions in which it is defined. Bin2D is not intended to be
  /// directly instantiated: it is inherited from to make specific histogram and
  /// profile bin types as HistoBin2D and ProfileBin2D.
  /// The lower bin edges in x and y are inclusive. This base class provides no fill
  /// method, since the signatures for standard and profile histos differ.
  template <class DBN>
  class Bin2D : public Bin {
  public:

    /// @name Constructors
    /// @{

    // /// Make a new, empty bin with two pairs of edges.
    // Bin2D(double xmin, double ymin, double xmax, double ymax)
    //   : _xedges( std::make_pair(xmin, xmax) ),
    //     _yedges( std::make_pair(ymin, ymax) )
    // {
    //   if (_xedges.second < _xedges.first) {
    //     throw RangeError("The bin x-edges are wrongly defined!");
    //   }
    //   if (_yedges.second < _yedges.first) {
    //     throw RangeError("The bin y-edges are wrongly defined!");
    //   }
    // }


    /// Make a new, empty bin with two pairs of edges
    Bin2D(const std::pair<double, double>& xedges, const std::pair<double, double>& yedges)
      : _xedges(xedges), _yedges(yedges)
    {
      if (_xedges.second < _xedges.first) {
        throw RangeError("The bin x-edges are wrongly defined!");
      }
      if (_yedges.second < _yedges.first) {
        throw RangeError("The bin y-edges are wrongly defined!");
      }
    }


    /// @brief Make a bin with all the components of a fill history
    ///
    /// Mainly intended for internal persistency use.
    Bin2D(const std::pair<double, double>& xedges,
          const std::pair<double, double>& yedges, const DBN& dbn)
      : _xedges(xedges), _yedges(yedges), _dbn(dbn)
    {
      if (_xedges.second < _xedges.first) {
        throw RangeError("The bin x-edges are wrongly defined!");
      }
      if (_yedges.second < _yedges.first) {
        throw RangeError("The bin y-edges are wrongly defined!");
      }
    }


    /// Copy constructor
    Bin2D(const Bin2D<DBN>& b)
      : _xedges(b._xedges),
        _yedges(b._yedges),
        _dbn(b._dbn)
    { }


    /// Copy assignment
    Bin2D<DBN>& operator = (const Bin2D<DBN>& b) {
      _xedges = b._xedges;
      _yedges = b._yedges;
      _dbn = b._dbn;
      return *this;
    }

    /// @}


    /// @name Dimensions
    /// @{

    /// Dimension of the fill space
    ///
    /// @todo Convert to total dimension
    size_t dim() { return 2; }

    /// Dimension of the fill space
    size_t fillDim() { return 2; }

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

    /// Scale the x and y coordinates and distributions.
    ///
    /// @note This should not be used, since it breaks histogram consistency. It will be removed in a future version.
    void scaleXY(double scaleX, double scaleY) {
      _xedges.first *= scaleX;
      _xedges.second *= scaleX;

      _yedges.first *= scaleY;
      _yedges.second *= scaleY;

      _dbn.scaleX(scaleX);
      _dbn.scaleY(scaleY);
    }

    /// @}


  public:

    /// @name X-axis info
    /// @{

    /// Get the {low,high} edges as an STL @c pair.
    std::pair<double,double> xEdges() const {
      return _xedges;
    }

    /// Lower x limit of the bin (inclusive).
    double xMin() const {
      return _xedges.first;
    }

    /// Upper x limit of the bin (exclusive).
    double xMax() const {
      return _xedges.second;
    }


    /// Get the {low,high} edges as an STL @c pair.
    std::pair<double,double> yEdges() const {
      return _yedges;
    }

    /// Lower y limit of the bin (inclusive).
    double yMin() const {
      return _yedges.first;
    }

    /// Upper y limit of the bin (exclusive).
    double yMax() const {
      return _yedges.second;
    }


    /// Middle of the bin in x
    double xMid() const {
      return (xMax() + xMin())/2.0;
    }

    /// Middle of the bin in y
    double yMid() const {
      return (yMax() + yMin())/2.0;
    }

    /// The geometric centre of the bin
    std::pair<double, double> xyMid() const {
      return std::make_pair(xMid(), yMid());
    }


    /// Width of the bin in x
    double xWidth() const {
      return xMax() - xMin();
    }

    /// Width of the bin in y
    double yWidth() const {
      return yMax() - yMin();
    }

    /// Widths of the bin in x and y
    std::pair<double, double> xyWidths() const {
      return std::make_pair(xWidth(), yWidth());
    }


    /// Area of the bin in x-y
    double area() const {
      return xWidth() * yWidth();
    }


    /// The mean x position in the bin, or the x midpoint if that is not available.
    double xFocus() const {
      return (!isZero(sumW())) ? xMean() : xMid();
    }

    /// The mean y position in the bin, or the y midpoint if that is not available.
    double yFocus() const {
      return (!isZero(sumW())) ? yMean() : yMid();
    }

    /// The mean position in the bin, or the midpoint if that is not available.
    std::pair<double, double> xyFocus() const {
      return std::make_pair(xFocus(), yFocus());
    }

    /// @}


  public:

    /// @name Distribution statistics
    /// @{

    /// Mean value of x-values in the bin.
    double xMean() const {
      return _dbn.xMean();
    }

    /// Mean value of y-values in the bin.
    double yMean() const {
      return _dbn.yMean();
    }

    /// The variance of x-values in the bin.
    double xVariance() const {
      return _dbn.xVariance();
    }

    /// The variance of y-values in the bin.
    double yVariance() const {
      return _dbn.yVariance();
    }

    /// The standard deviation (spread) of x-values in the bin.
    double xStdDev() const {
      return _dbn.xStdDev();
    }

    /// The standard deviation (spread) of y-values in the bin.
    double yStdDev() const {
      return _dbn.yStdDev();
    }

    /// The standard error on the bin x focus.
    double xStdErr() const {
      return _dbn.xStdErr();
    }

    /// The standard error on the bin y focus.
    double yStdErr() const {
      return _dbn.yStdErr();
    }

    /// The x RMS in the bin.
    double xRMS() const {
      return _dbn.xRMS();
    }

    /// The y RMS in the bin.
    double yRMS() const {
      return _dbn.yRMS();
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

    /// The sum of y*weight
    double sumWY() const {
      return _dbn.sumWY();
    }

    /// The sum of x*y*weight
    double sumWXY() const {
      return _dbn.sumWXY();
    }

    /// The sum of x^2 * weight
    double sumWX2() const {
      return _dbn.sumWX2();
    }

    /// The sum of y^2 * weight
    double sumWY2() const {
      return _dbn.sumWY2();
    }

    /// @}


  public:

    /// @name Operators
    /// @{

    /// Add two bins
    Bin2D<DBN>& operator += (const Bin2D<DBN>& b) {
      return add(b);
    }

    /// Subtract one bin from another
    Bin2D<DBN>& operator -= (const Bin2D<DBN>& b) {
      return subtract(b);
    }

    /// @}


    /// @name Named operators
    /// @{

    /// Merge two adjacent bins
    /*
    Bin2D<DBN>& merge(const Bin2D<DBN>& b) {

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
    */


    /// Add two bins (internal, explicitly named version)
    ///
    /// This operator is defined for adding two bins with equivalent binning.
    /// It cannot be used to merge two bins into one larger bin.
    Bin2D<DBN>& add(const Bin2D<DBN>& b) {
      if (!fuzzyEquals(_xedges.first, b._xedges.first) ||
          !fuzzyEquals(_xedges.second, b._xedges.second) ||
          !fuzzyEquals(_yedges.first, b._yedges.first) ||
          !fuzzyEquals(_yedges.second, b._yedges.second)) {
        throw LogicError("Attempted to add two bins with different edges");
      }
      _dbn += b._dbn;
      return *this;
    }


    /// Subtract one bin from another (internal, explicitly named version)
    ///
    /// This operator is defined for subtracting two bins with equivalent binning.
    /// It cannot be used to merge two bins into one larger bin.
    Bin2D<DBN>& subtract(const Bin2D<DBN>& b) {
      if (!fuzzyEquals(_xedges.first, b._xedges.first) ||
          !fuzzyEquals(_xedges.second, b._xedges.second) ||
          !fuzzyEquals(_yedges.first, b._yedges.first) ||
          !fuzzyEquals(_yedges.second, b._yedges.second)) {
        throw LogicError("Attempted to subtract two bins with different edges");
      }
      _dbn -= b._dbn;
      return *this;
    }

    /// Test whether this bin would fit inside the given area.
    bool fitsInside(std::pair<double, double> xrange,
                    std::pair<double, double> yrange) const {
      return (xMin() >= xrange.first &&
              yMin() >= yrange.first &&
              xMax() < xrange.second &&
              yMax() < yrange.second);
    }

    /// Test whether a point lies within the current bin
    bool bounds(double x, double y) const {
      return (x >= xMin() && x < xMax() && y >= yMin() && y < yMax());
    }


    /// Test whether two bins are adjacent and, if so, return how as an integer.
    int adjacentTo(const Bin2D<DBN> &b) const {
      for (int i = 0; i < 4; i++) {
        if (_edges_equal(b, i, (i+2) % 4))
          return i;
      }
      return -1;
    }

    /// @}


  protected:

    /// @todo Remove?
    std::pair<double, double> _edge_par(int i) const {
      if (i % 2)
        return xEdges();
      else
        return yEdges();
    }

    /// @todo Remove?
    double _edge_perp(size_t i) const {
      double output = 0.0;
      switch (i) {
        case 0: output = xMax(); break;
        case 1: output = yMax(); break;
        case 2: output = xMin(); break;
        case 3: output = yMin(); break;
      }
      return output;
    }

    // Check if common edge.
    /// @todo Remove?
    bool _edges_equal(const Bin2D<DBN>& other, const int i, const int j) const {
      return other._edges_equal(_edge_perp(i), _edge_par(i), j);
    }

    /// @todo Remove?
    bool _edges_equal(const double perp, const std::pair<double, double> par, int j) const {
      return (fuzzyEquals(perp, _edge_perp(j)) &&
              fuzzyEquals(par.first, _edge_par(j).first) &&
              fuzzyEquals(par.second, _edge_par(j).second));
    }


  protected:

    /// The bin limits
    std::pair<double,double> _xedges;
    std::pair<double,double> _yedges;

    // Distribution of weighted x (and perhaps y) values
    DBN _dbn;

    /// @todo Make Axis2D<BIN, DBN> -> Axis2D<DBN> and hold a pointer to the one that contains this bin

  };





  /// Add two bins
  ///
  /// This "add" operator is defined for adding two bins with equivalent binning.
  /// It cannot be used to merge two bins into one larger bin.
  template <class DBN>
  inline Bin2D<DBN> operator + (const Bin2D<DBN>& a, const Bin2D<DBN>& b) {
    Bin2D<DBN> rtn = a;
    rtn += b;
    return rtn;
  }


  /// Subtract one bin from another
  ///
  /// This "subtraction" operator is defined for subtracting two bins with equivalent binning.
  /// It cannot be used to merge two bins into one larger bin.
  template <class DBN>
  inline Bin2D<DBN> operator - (const Bin2D<DBN>& a, const Bin2D<DBN>& b) {
    Bin2D<DBN> rtn = a;
    rtn -= b;
    return rtn;
  }


  /// Bin2Ds are compared for axis sorting by lower edge position in first x and then y directions
  template <class DBN>
  inline bool operator<(const Bin2D<DBN>& a, const Bin2D<DBN>& b) {
    if (!fuzzyEquals(a.xMin(), b.xMin())) return b.xMin() > a.xMin();
    return b.yMin() > a.yMin();
  }



}



#endif
