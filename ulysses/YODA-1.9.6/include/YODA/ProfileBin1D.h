// -*- C++ -*-
//
// This file is part of YODA -- Yet more Objects for Data Analysis
// Copyright (C) 2008-2021 The YODA collaboration (see AUTHORS for details)
//
#ifndef YODA_ProfileBin1D_h
#define YODA_ProfileBin1D_h

#include "YODA/Bin1D.h"
#include "YODA/Dbn2D.h"
#include "YODA/Exceptions.h"

namespace YODA {


  /// @brief A Bin1D specialised for handling profile-type information
  ///
  /// This is a 1D bin type, which supports all the operations defined for
  /// a generic Bin1D object, but also supplies the specific member functions
  /// for profile-type data, as opposed to histogram-type. This means that
  /// extra internal distribution statistics are stored for the extra
  /// "y-direction" specified in the profile fill operation.
  class ProfileBin1D : public Bin1D<Dbn2D> {
  public:

    /// @name Constructors
    /// @{

    /// Constructor giving bin low and high edges.
    ProfileBin1D(double lowedge, double highedge)
      : Bin1D<Dbn2D>(std::make_pair(lowedge, highedge))
    { }


    /// Constructor giving bin low and high edges as a pair.
    ProfileBin1D(const std::pair<double,double>& edges)
      : Bin1D<Dbn2D>(edges)
    { }


    /// @brief Make a profile bin with all the components of a fill history.
    ///
    /// Mainly intended for internal persistency use.
    ProfileBin1D(std::pair<double, double> edges, const Dbn2D& dbnxy)
      : Bin1D<Dbn2D>(edges, dbnxy)
    {  }


    /// Copy constructor
    ProfileBin1D(const ProfileBin1D& pb)
      : Bin1D<Dbn2D>(pb)
    { }


    /// Copy assignment
    ProfileBin1D& operator = (const ProfileBin1D& pb) {
      Bin1D<Dbn2D>::operator=(pb);
      return *this;
    }

    /// @}


    /// @name Modifiers
    /// @{

    /// Fill histo by x and y values and weight.
    ///
    /// @note This should not be used, since it breaks histogram consistency. It will be removed in a future version.
    void fill(double x, double y, double weight=1.0, double fraction=1.0) {
      _dbn.fill(x, y, weight, fraction);
    }

    /// Fill histo with @a weight and y-value @c y at x = bin midpoint.
    ///
    /// @note This should not be used, since it breaks histogram consistency. It will be removed in a future version.
    void fillBin(double y, double weight=1.0, double fraction=1.0) {
      fill(xMid(), y, weight, fraction);
    }

    /// @}


    /// @name Bin scaling (x scaling is inherited)
    /// @{

    /// Scale the y (profiled) dimension
    ///
    /// @note This should not be used, since it breaks histogram consistency. It will be removed in a future version.
    inline void scaleY(double ay) {
      _dbn.scaleY(ay);
    }

    /// Scale the x and y dimensions
    ///
    /// @note This should not be used, since it breaks histogram consistency. It will be removed in a future version.
    inline void scaleXY(double ax, double ay) {
      scaleX(ax);
      scaleY(ay);
    }

    /// @}


  public:

    /// @name Bin content info
    /// @{

    /// The mean of the y distribution
    double mean() const {
      return _dbn.yMean();
    }

    /// The std deviation of the y distribution about the mean
    double stdDev() const {
      return _dbn.yStdDev();
    }

    /// The variance of the y distribution about the mean
    double variance() const {
      return _dbn.yVariance();
    }

    /// The standard error on the mean
    double stdErr() const {
      return _dbn.yStdErr();
    }

    /// The relative size of the error on the mean
    double relErr() const {
      return stdErr() != 0 ? stdErr() / mean() : 0;
    }

    /// The RMS of the y distribution
    double rms() const {
      return _dbn.yRMS();
    }

    /// @}


    /// @name Raw y distribution statistics
    /// @{

    /// The sum of y*weight
    double sumWY() const {
      return _dbn.sumWY();
    }

    /// The sum of y^2 * weight
    double sumWY2() const {
      return _dbn.sumWY2();
    }

    /// @}


  public:

    /// Add two bins (for use by Profile1D).
    ProfileBin1D& operator += (const ProfileBin1D& toAdd) {
      return add(toAdd);
    }

    /// Subtract two bins
    ProfileBin1D& operator -= (const ProfileBin1D& toSubtract) {
      return subtract(toSubtract);
    }


  protected:

    /// Add two bins (internal, explicitly named version)
    ProfileBin1D& add(const ProfileBin1D& pb) {
      Bin1D<Dbn2D>::add(pb);
      return *this;
    }

    /// Subtract one bin from another (internal, explicitly named version)
    ProfileBin1D& subtract(const ProfileBin1D& pb) {
      Bin1D<Dbn2D>::subtract(pb);
      return *this;
    }

  };


  inline ProfileBin1D operator + (const ProfileBin1D& a, const ProfileBin1D& b) {
    ProfileBin1D rtn(a);
    rtn += b;
    return rtn;
  }

  inline ProfileBin1D operator - (const ProfileBin1D& a, const ProfileBin1D& b) {
    ProfileBin1D rtn(a);
    rtn -= b;
    return rtn;
  }


}

#endif
