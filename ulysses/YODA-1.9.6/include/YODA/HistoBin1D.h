// -*- C++ -*-
//
// This file is part of YODA -- Yet more Objects for Data Analysis
// Copyright (C) 2008-2021 The YODA collaboration (see AUTHORS for details)
//
#ifndef YODA_HistoBin1D_h
#define YODA_HistoBin1D_h

#include "YODA/Bin1D.h"
#include "YODA/Dbn1D.h"
#include "YODA/Exceptions.h"

namespace YODA {


  /// @brief A Bin1D specialised for handling histogram-type information
  ///
  /// This is a 1D bin type, which supports all the operations defined for
  /// a generic Bin1D object, but also supplies the specific member functions
  /// for histogram-type data, as opposed to profile-type.
  class HistoBin1D : public Bin1D<Dbn1D> {
  public:

    /// @name Constructor giving bin low and high edges.
    /// @{

    /// Make a new, empty bin with a pair of edges.
    HistoBin1D(double lowedge, double highedge)
      : Bin1D<Dbn1D>(std::make_pair(lowedge, highedge))
    { }


    /// Make a new, empty bin with a pair of edges.
    HistoBin1D(const std::pair<double,double>& edges)
      : Bin1D<Dbn1D>(edges)
    { }


    /// @brief Make a bin with all the components of a fill history.
    ///
    /// Mainly intended for internal persistency use.
    HistoBin1D(std::pair<double, double> edges, const Dbn1D& dbnx)
      : Bin1D<Dbn1D>(edges, dbnx)
    { }


    /// Copy constructor
    HistoBin1D(const HistoBin1D& hb)
      : Bin1D<Dbn1D>(hb)
    { }


    /// Copy assignment
    HistoBin1D& operator = (const HistoBin1D& hb) {
      Bin1D<Dbn1D>::operator=(hb);
      return *this;
    }

    /// @}


  public:

    /// @name Modifiers
    /// @{

    /// Fill this bin with weight @a weight at position @a x.
    ///
    /// @note This should not be used, since it breaks histogram consistency. It will be removed in a future version.
    void fill(double x, double weight=1.0, double fraction=1.0) {
      _dbn.fill(x, weight, fraction);
    }

    /// Fill this bin with weight @a weight.
    ///
    /// @note This should not be used, since it breaks histogram consistency. It will be removed in a future version.
    void fillBin(double weight=1.0, double fraction=1.0) {
      fill(xMid(), weight, fraction);
    }

    /// @}


  public:

    /// @name Bin content info
    /// @{

    /// The area is the sum of weights in the bin, i.e. the
    /// width of the bin has no influence on this figure.
    double area() const {
      return sumW();
    }

    /// The height is defined as area/width.
    double height() const {
      return area() / xWidth();
    }

    /// @}


    /// @name Error info
    /// @{

    /// Error computed using binomial statistics on the sum of bin weights,
    /// i.e. err_area = sqrt{sum{weights}}
    double areaErr() const {
      return sqrt(sumW2());
    }

    /// As for the height vs. area, the height error includes a scaling factor
    /// of the bin width, i.e. err_height = sqrt{sum{weights}} / width.
    double heightErr() const {
      return areaErr() / xWidth();
    }

    /// The relative size of the error (same for either area or height errors)
    double relErr() const {
      /// @todo Throw excp if sumW2 is 0?
      return sumW2() != 0 ? sqrt(sumW2()) / sumW() : 0;
    }

    /// @}


  public:

    /// Add two bins (for use by Histo1D).
    HistoBin1D& operator += (const HistoBin1D& toAdd) {
      return add(toAdd);
    }

    /// Subtract two bins
    HistoBin1D& operator -= (const HistoBin1D& toSubtract) {
      return subtract(toSubtract);
    }


  protected:

    /// Add two bins (internal, explicitly named version)
    HistoBin1D& add(const HistoBin1D& hb) {
      Bin1D<Dbn1D>::add(hb);
      return *this;
    }

    /// Subtract one bin from another (internal, explicitly named version)
    HistoBin1D& subtract(const HistoBin1D& hb) {
      Bin1D<Dbn1D>::subtract(hb);
      return *this;
    }

  };


  /// Add two bins
  inline HistoBin1D operator + (const HistoBin1D& a, const HistoBin1D& b) {
    HistoBin1D rtn(a);
    rtn += b;
    return rtn;
  }

  /// Subtract two bins
  inline HistoBin1D operator - (const HistoBin1D& a, const HistoBin1D& b) {
    HistoBin1D rtn(a);
    rtn -= b;
    return rtn;
  }


}

#endif
