// -*- C++ -*-
//
// This file is part of YODA -- Yet more Objects for Data Analysis
// Copyright (C) 2008-2020 The YODA collaboration (see AUTHORS for details)
//
#ifndef YODA_Scatter_h
#define YODA_Scatter_h

#include "YODA/AnalysisObject.h"

namespace YODA {


  /// A base class for common operations on scatter types (Scatter1D, etc.)
  class Scatter {
  public:

    /// @todo Add a generic Scatter base class, providing reset(), rmPoint(), etc.

    /// Virtual destructor for inheritance
    virtual ~Scatter() {}

    //@}


    /// Dimension of this data object
    virtual size_t dim() const = 0;


    /// @name Modifiers
    //@{

    /// Clear all points
    virtual void reset() = 0;

    /// Scaling along direction @a i
    virtual void scale(size_t i, double scale) = 0;

    //@}


    ///////////////////////////////////////////////////


    // virtual void parseVariations() const = 0;

    /// Get the list of variations stored in the points
    virtual std::vector<std::string> variations() const = 0;

    /// Clear the variations
    virtual void rmVariations() = 0;

    /// Number of points in the scatter
    virtual size_t numPoints() const = 0;


    /// @name Point removers
    /// @{

    /// Remove the point with index @a index
    virtual void rmPoint(size_t index) = 0;

    /// Safely remove the points with indices @a indices
    virtual void rmPoints(std::vector<size_t> indices) {
      // remove points in decreasing order, so the numbering isn't invalidated mid-loop:
      std::sort(indices.begin(), indices.end(), std::greater<size_t>()); //< reverse-sort
      for (size_t i : indices) rmPoint(i);
    }

    /// @}


    /// @name Point inserters
    ///
    /// These don't currently make sense for a generic base class, but passing
    /// generic tuples of vals/errs would be quite nice.
    ///
    /// @{

    // /// Insert a new point, defined as the x value and no errors
    // void addPoint(double x) {
    //   Point1D thisPoint=Point1D(x);
    //   thisPoint.setParentAO(this);
    //   _points.insert(thisPoint);
    // }

    // /// Insert a new point, defined as the x value and symmetric errors
    // void addPoint(double x, double ex) {
    //   Point1D thisPoint=Point1D(x, ex);
    //   thisPoint.setParentAO(this);
    //   _points.insert(thisPoint);
    // }

    // /// Insert a new point, defined as the x value and an asymmetric error pair
    // void addPoint(double x, const std::pair<double,double>& ex) {
    //   Point1D thisPoint=Point1D(x, ex);
    //   thisPoint.setParentAO(this);
    //   _points.insert(thisPoint);
    // }

    // /// Insert a new point, defined as the x value and explicit asymmetric errors
    // void addPoint(double x, double exminus, double explus) {
    //   Point1D thisPoint=Point1D(x, exminus, explus);
    //   thisPoint.setParentAO(this);
    //   _points.insert(thisPoint);
    // }

    // /// Insert a collection of new points
    // void addPoints(const Points& pts) {
    //   for (const Point1D& pt : pts) addPoint(pt);
    // }

    //@}


    // /// Equality operator
    // bool operator == (const Scatter1D& other) {
    //   return _points == other._points;
    // }

    // /// Non-equality operator
    // bool operator != (const Scatter1D& other) {
    //   return ! operator == (other);
    // }

  };


}

#endif
