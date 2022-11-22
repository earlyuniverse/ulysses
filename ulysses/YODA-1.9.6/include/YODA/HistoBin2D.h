// -*- C++ -*-
//
// This file is part of YODA -- Yet more Objects for Data Analysis
// Copyright (C) 2008-2021 The YODA collaboration (see AUTHORS for details)
//
#ifndef YODA_HistoBin2D_h
#define YODA_HistoBin2D_h

#include "YODA/Bin2D.h"
#include "YODA/Dbn2D.h"
#include "YODA/Exceptions.h"

namespace YODA {


  /// @brief A Bin2D specialised for handling histogram-type information
  ///
  /// This is a 2D bin type, which supports all the operations defined for
  /// a generic Bin2D object, but also supplies the specific member functions
  /// for histogram-type data, as opposed to profile-type.
  class HistoBin2D : public Bin2D<Dbn2D> {
  public:

    /// @name Constructors
    /// @{

    /// Make a new, empty bin with two pairs of edges.
    HistoBin2D(double xmin, double xmax, double ymin, double ymax)
      : Bin2D<Dbn2D>(std::make_pair(xmin, xmax), std::make_pair(ymin, ymax))
    { }

    /// Constructor accepting a set of all edges of a bin
    HistoBin2D(const std::pair<double,double>& xedges,
               const std::pair<double,double>& yedges)
      : Bin2D<Dbn2D>(xedges, yedges)
    { }

    /// @brief Make a bin with all the components of a fill history.
    ///
    /// Mainly intended for internal persistency use.
    HistoBin2D(const std::pair<double, double>& xedges,
               const std::pair<double, double>& yedges, const Dbn2D& dbn)
      : Bin2D<Dbn2D>(xedges, yedges, dbn)
    { }

    /// Copy constructor
    HistoBin2D(const HistoBin2D& pb)
      : Bin2D<Dbn2D>(pb)
    { }

    /// Copy assignment
    HistoBin2D& operator = (const HistoBin2D& hb) {
      Bin2D<Dbn2D>::operator=(hb);
      return *this;
    }

    /// @}


    /// @name Modifiers
    /// @{

    /// A fill() function accepting coordinates as separate numbers
    ///
    /// @note This should not be used, since it breaks histogram consistency. It will be removed in a future version.
    void fill(double x, double y, double weight=1.0, double fraction=1.0) {
      _dbn.fill(x, y, weight, fraction);
    }

    /// A fill() function accepting the coordinates as std::pair
    ///
    /// @note This should not be used, since it breaks histogram consistency. It will be removed in a future version.
    void fill(std::pair<double,double> coords, double weight=1.0, double fraction=1.0) {
      fill(coords.first, coords.second, weight, fraction);
    }

    /// A function that fills this particular bin.
    ///
    /// @note This should not be used, since it breaks histogram consistency. It will be removed in a future version.
    void fillBin(double weight=1.0, double fraction=1.0) {
      fill(xyMid(), weight, fraction);
    }

    /// A reset function
    void reset() {
      Bin2D<Dbn2D>::reset();
    }

    /// @}


    /// @name Accessors
    /// @{

    /// The volume of a bin
    double volume() const {
      return sumW();
    }

    /// Error on volume
    double volumeErr() const {
      return sqrt(sumW2());
    }

    /// The height of a bin
    double height() const {
      return volume()/(xWidth()*yWidth());
    }

    /// Error on height
    double heightErr() const {
      return volumeErr()/(xWidth()*yWidth());
    }

    /// The relative size of the error (same for either volume or height errors)
    double relErr() const {
      return sumW2() != 0 ? sqrt(sumW2()) / sumW() : 0;
    }

    /// @}


    /// @name Transformers
    /// @{

    // /// @brief Transformer taking x as the primary axis of ProfileBin1D
    // ///
    // /// @todo Need to think about the name, and clarify what "primary axis" means
    // ProfileBin1D transformX() {
    //   ProfileBin1D ret(std::make_pair(xMin(), xMax()), Dbn2D(_dbn));
    //   return ret;
    // }

    // /// @brief Transformer taking y as the primary axis of ProfileBin1D
    // ///
    // /// @todo Need to think about the name, and clarify what "primary axis" means
    // ProfileBin1D transformY() {
    //   Dbn2D dbn = _dbn; dbn.flipXY();
    //   ProfileBin1D ret(std::make_pair(yMin(), yMax()), Dbn2D(dbn));
    //   return ret;
    // }

    /// @}
  };


  /// Bin addition operator
  inline HistoBin2D operator + (const HistoBin2D& a, const HistoBin2D& b) {
    HistoBin2D rtn(a);
    rtn += b;
    return rtn;
  }


  /// Bin subtraction operator
  inline HistoBin2D operator - (const HistoBin2D& a, const HistoBin2D& b) {
    HistoBin2D rtn(a);
    rtn -= b;
    return rtn;
  }


}

#endif
