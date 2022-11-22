// -*- C++ -*-
//
// This file is part of YODA -- Yet more Objects for Data Analysis
// Copyright (C) 2008-2021 The YODA collaboration (see AUTHORS for details)
//
#ifndef YODA_Dbn3D_h
#define YODA_Dbn3D_h

#include "YODA/Exceptions.h"
#include "YODA/Dbn1D.h"

namespace YODA {


  /// A 2D distribution
  class Dbn3D {
  public:

    /// @name Constructors
    /// @{

    /// Default constructor of a new distribution.
    Dbn3D() {
      reset();
    }


    /// Constructor to set a distribution with a pre-filled state.
    ///
    /// Principally designed for internal persistency use.
    Dbn3D(double numEntries,
          double sumW, double sumW2,
          double sumWX, double sumWX2,
          double sumWY, double sumWY2,
          double sumWZ, double sumWZ2,
          double sumWXY, double sumWXZ, double sumWYZ)
      : _dbnX(numEntries, sumW, sumW2, sumWX, sumWX2),
        _dbnY(numEntries, sumW, sumW2, sumWY, sumWY2),
        _dbnZ(numEntries, sumW, sumW2, sumWZ, sumWZ2),
        _sumWXY(sumWXY), _sumWXZ(sumWXZ), _sumWYZ(sumWYZ)
    { }


    /// Copy constructor
    ///
    /// Sets all the parameters using the ones provided from an existing Dbn3D.
    Dbn3D(const Dbn3D& toCopy) {
      _dbnX = toCopy._dbnX;
      _dbnY = toCopy._dbnY;
      _dbnZ = toCopy._dbnZ;
      _sumWXY = toCopy._sumWXY;
      _sumWXZ = toCopy._sumWXZ;
      _sumWYZ = toCopy._sumWYZ;
    }


    /// Copy assignment
    ///
    /// Sets all the parameters using the ones provided from an existing Dbn3D.
    Dbn3D& operator=(const Dbn3D& toCopy) {
      _dbnX = toCopy._dbnX;
      _dbnY = toCopy._dbnY;
      _dbnZ = toCopy._dbnZ;
      _sumWXY = toCopy._sumWXY;
      _sumWXZ = toCopy._sumWXZ;
      _sumWYZ = toCopy._sumWYZ;
      return *this;
    }

    /// @}


    /// @name Modifiers
    /// @{

    /// Fill, providing the fill coordinates as three different numbers.
    void fill(double valX, double valY, double valZ, double weight=1.0, double fraction=1.0) {
      _dbnX.fill(valX, weight, fraction);
      _dbnY.fill(valY, weight, fraction);
      _dbnZ.fill(valZ, weight, fraction);
      _sumWXY += fraction*weight*valX*valY;
      _sumWXZ += fraction*weight*valX*valZ;
      _sumWYZ += fraction*weight*valY*valZ;
    }


    /// Fill, providing the fill coordinates as a vector.
    void fill(std::vector<double> val, double weight=1.0, double fraction=1.0) {
      assert (val.size() == 3);
      fill(val[0], val[1], val[2], weight, fraction);
    }


    /// Reset the distribution to an unfilled state.
    void reset() {
      _dbnX.reset();
      _dbnY.reset();
      _dbnZ.reset();
      _sumWXY = 0;
      _sumWXZ = 0;
      _sumWYZ = 0;
    }


    /// Rescale as if all fill weights had been different by factor @a scalefactor.
    void scaleW(double scalefactor) {
      _dbnX.scaleW(scalefactor);
      _dbnY.scaleW(scalefactor);
      _dbnZ.scaleW(scalefactor);
      _sumWXY *= scalefactor;
      _sumWXZ *= scalefactor;
      _sumWYZ *= scalefactor;
    }


    /// Rescale x: needed if x histo bin edges are rescaled.
    void scaleX(double xscale) {
      _dbnX.scaleX(xscale);
      _sumWXY *= xscale;
      _sumWXZ *= xscale;
    }


    /// Rescale y: needed if y histo bin edges are rescaled.
    void scaleY(double yscale) {
      _dbnY.scaleX(yscale);
      _sumWXY *= yscale;
      _sumWYZ *= yscale;
    }


    /// Rescale z: needed if z histo bin edges are rescaled.
    void scaleZ(double zscale) {
      _dbnZ.scaleX(zscale);
      _sumWXZ *= zscale;
      _sumWYZ *= zscale;
    }


    // /// Rescale x and y: needed if histo bin edges are rescaled.
    // void scaleXY(double xscale, double yscale) {
    //   scaleX(xscale);
    //   scaleY(yscale);
    // }


    // /// Rescale x and z: needed if histo bin edges are rescaled.
    // void scaleXZ(double xscale, double zscale) {
    //   scaleX(xscale);
    //   scaleZ(zscale);
    // }


    // /// Rescale y and z: needed if histo bin edges are rescaled.
    // void scaleYZ(double yscale, double zscale) {
    //   scaleY(yscale);
    //   scaleZ(zscale);
    // }


    /// Rescale x, y and z: needed if histo bin edges are rescaled.
    void scaleXYZ(double xscale, double yscale, double zscale) {
      scaleX(xscale);
      scaleY(yscale);
      scaleZ(zscale);
    }

    /// @}

  public:


    /// @name Distribution statistics
    /// @{

    /// The absolute error on sumW
    double errW() const { return _dbnX.errW(); }

    /// The relative error on sumW
    double relErrW() const { return _dbnX.relErrW(); }

    /// Weighted mean, \f$ \bar{x} \f$, of distribution.
    double xMean() const { return _dbnX.xMean(); }

    /// Weighted mean, \f$ \bar{y} \f$, of distribution.
    double yMean() const { return _dbnY.xMean(); }

    /// Weighted mean, \f$ \bar{z} \f$, of distribution.
    double zMean() const { return _dbnZ.xMean(); }

    /// Weighted \f$ x \f$ variance, \f$ \sigma_x^2 \f$, of distribution.
    double xVariance() const { return _dbnX.xVariance(); }

    /// Weighted \f$ y \f$ variance, \f$ \sigma_y^2 \f$, of distribution.
    double yVariance() const { return _dbnY.xVariance(); }

    /// Weighted \f$ z \f$ variance, \f$ \sigma_z^2 \f$, of distribution.
    double zVariance() const { return _dbnZ.xVariance(); }

    /// Weighted \f$ x \f$ standard deviation, \f$ \sigma_x \f$, of distribution.
    double xStdDev() const { return _dbnX.xStdDev(); }

    /// Weighted \f$ y \f$ standard deviation, \f$ \sigma_y \f$, of distribution.
    double yStdDev() const { return _dbnY.xStdDev(); }

    /// Weighted \f$ z \f$ standard deviation, \f$ \sigma_z \f$, of distribution.
    double zStdDev() const { return _dbnZ.xStdDev(); }

    /// Weighted standard error on the \f$ x \f$ mean, \f$ \sim \sigma_x/\sqrt{N-1} \f$, of distribution.
    double xStdErr() const { return _dbnX.xStdErr(); }

    /// Weighted standard error on the \f$ y \f$ mean, \f$ \sim \sigma_y/\sqrt{N-1} \f$, of distribution.
    double yStdErr() const { return _dbnY.xStdErr(); }

    /// Weighted standard error on the \f$ z \f$ mean, \f$ \sim \sigma_z/\sqrt{N-1} \f$, of distribution.
    double zStdErr() const { return _dbnZ.xStdErr(); }

    /// Weighted RMS, \f$ \sqrt{ \sum{w x^2}/\sum{w} } \f$, of distribution.
    double xRMS() const { return _dbnX.xRMS(); }

    /// Weighted RMS, \f$ \sqrt{ \sum{w y^2}/\sum{w} } \f$, of distribution.
    double yRMS() const { return _dbnY.xRMS(); }

    /// Weighted RMS, \f$ \sqrt{ \sum{w z^2}/\sum{w} } \f$, of distribution.
    double zRMS() const { return _dbnZ.xRMS(); }

    /// @}


    /// @name Raw distribution running sums
    /// @{

    /// Number of entries (number of times @c fill was called, ignoring weights)
    double numEntries() const {
      return _dbnX.numEntries();
    }

    /// Effective number of entries \f$ = (\sum w)^2 / \sum w^2 \f$
    double effNumEntries() const {
      return _dbnX.effNumEntries();
    }

    /// The sum of weights
    double sumW() const {
      return _dbnX.sumW();
    }

    /// The sum of weights squared
    double sumW2() const {
      return _dbnX.sumW2();
    }

    /// The sum of x*weight
    double sumWX() const {
      return _dbnX.sumWX();
    }

    /// The sum of x^2*weight
    double sumWX2() const {
      return _dbnX.sumWX2();
    }

    /// The sum of y*weight
    double sumWY() const {
      return _dbnY.sumWX();
    }

    /// The sum of y^2*weight
    double sumWY2() const {
      return _dbnY.sumWX2();
    }

    /// The sum of z*weight
    double sumWZ() const {
      return _dbnZ.sumWX();
    }

    /// The sum of z^2*weight
    double sumWZ2() const {
      return _dbnZ.sumWX2();
    }

    /// The sum of x*y*weight
    double sumWXY() const {
      return _sumWXY;
    }

    /// The sum of x*z*weight
    double sumWXZ() const {
      return _sumWXZ;
    }

    /// The sum of y*z*weight
    double sumWYZ() const {
      return _sumWXZ;
    }

    /// @}


    /// @name Operators
    /// @{

    /// Add two dbns
    Dbn3D& operator += (const Dbn3D& d) {
      return add(d);
    }

    /// Subtract one dbn from another
    Dbn3D& operator -= (const Dbn3D& d) {
      return subtract(d);
    }

    /// @brief Interchange X and Y subdistributions
    ///
    /// Mostly used for operations on total distribution of an Axis
    void flipXY() {
      Dbn1D temp(_dbnX);
      _dbnX = _dbnY;
      _dbnY = temp;
    }

    /// @brief Interchange X and Z subdistributions
    ///
    /// Mostly used for operations on total distribution of an Axis
    void flipXZ() {
      Dbn1D temp(_dbnX);
      _dbnX = _dbnZ;
      _dbnZ = temp;
    }

    /// @brief Interchange Y and Z subdistributions
    ///
    /// Mostly used for operations on total distribution of an Axis
    void flipYZ() {
      Dbn1D temp(_dbnY);
      _dbnY = _dbnZ;
      _dbnZ = temp;
    }

    /// @brief Transform into a Dbn1D parallel to X axis (dropping Y and Z terms)
    ///
    /// @todo Rename
    Dbn1D transformX() {
      Dbn1D ret(_dbnX);
      return ret;
    }

    /// @brief Transform into a Dbn1D parallel to Y axis (dropping X and Z terms)
    ///
    /// @todo Rename
    Dbn1D transformY() {
      Dbn1D ret(_dbnY);
      return ret;
    }

    /// @brief Transform into a Dbn1D parallel to Z axis (dropping X and Y terms)
    ///
    /// @todo Rename
    Dbn1D transformZ() {
      Dbn1D ret(_dbnZ);
      return ret;
    }
    /// @}


  protected:

    /// Add two dbns (internal, explicitly named version)
    Dbn3D& add(const Dbn3D& d) {
      _dbnX += d._dbnX;
      _dbnY += d._dbnY;
      _dbnZ += d._dbnZ;
      _sumWXY += d._sumWXY;
      _sumWXZ += d._sumWXZ;
      _sumWYZ += d._sumWYZ;
      return *this;
    }

    /// Subtract one dbn from another (internal, explicitly named version)
    Dbn3D& subtract(const Dbn3D& d) {
      _dbnX -= d._dbnX;
      _dbnY -= d._dbnY;
      _dbnZ -= d._dbnZ;
      _sumWXY -= d._sumWXY;
      _sumWXZ -= d._sumWXZ;
      _sumWYZ -= d._sumWYZ;
      return *this;
    }


  private:

    /// @name Storage
    /// @{

    /// The x moments and the pure-weight quantities are stored in a 1D "x" distribution
    Dbn1D _dbnX;

    /// The y moments are stored in a 1D "y" distribution
    Dbn1D _dbnY;

    /// The z moments are stored in a 1D "z" distribution
    Dbn1D _dbnZ;

    /// The higher-order "cross-term" that can't be handled using the 1D distributions
    double _sumWXY;
    double _sumWXZ;
    double _sumWYZ;
    /// @}

  };


  /// Add two dbns
  inline Dbn3D operator + (const Dbn3D& a, const Dbn3D& b) {
    Dbn3D rtn = a;
    rtn += b;
    return rtn;
  }

  /// Subtract one dbn from another
  inline Dbn3D operator - (const Dbn3D& a, const Dbn3D& b) {
    Dbn3D rtn = a;
    rtn -= b;
    return rtn;
  }


}

#endif
