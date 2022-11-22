// -*- C++ -*-
//
// This file is part of YODA -- Yet more Objects for Data Analysis
// Copyright (C) 2008-2021 The YODA collaboration (see AUTHORS for details)
//
#ifndef YODA_SCATTERND_H
#define YODA_SCATTERND_H

#include "YODA/AnalysisObject.h"
#include "YODA/PointND.h"
#include "YODA/Utils/sortedvector.h"
#include "YODA/Utils/ndarray.h"
#include <vector>
#include <set>
#include <string>
#include <utility>
#include <memory>

namespace YODA {


  /// A very generic data type which is just a collection of ND data points with errors
  template <int N>
  class Scatter : public AnalysisObject {
  public:

    // Typedefs
    typedef Utils::ndarray<double, N> NdVal;
    typedef Utils::ndarray<std::pair<double,double>, N> NdValPair;
    typedef Utils::sortedvector< Point<N> > Points;
    typedef std::shared_ptr<Scatter> Ptr;


    /// @name Constructors
    /// @{

    /// Empty constructor
    Scatter(const std::string& path="", const std::string& title="")
      : AnalysisObject("Scatter", path, title)
    {  }


    /// Constructor from a set of points
    Scatter(const Points& points,
            const std::string& path="", const std::string& title="")
      : AnalysisObject("Scatter", path, title),
        _points(points)
    {  }


    /// Constructor from a vector of position values with no errors
    Scatter(const std::vector<NdVal>& positions,
            const std::string& path="", const std::string& title="")
      : AnalysisObject("Scatter", path, title)
    {
      for (size_t i = 0; i < positions.size(); ++i) {
        addPoint(Point<N>(positions[i]));
      }
    }


    /// Constructor from vectors of values for positions and a single set of symmetric errors
    Scatter(const std::vector<NdVal>& positions,
            const std::vector<NdVal>& errors,
            const std::string& path="", const std::string& title="")
      : AnalysisObject("Scatter", path, title)
    {
      assert(positions.size() == errors.size());
      for (size_t i = 0; i < positions.size(); ++i) {
        addPoint(Point<N>(positions[i], errors[i]));
      }
    }


    // /// Constructor from values with completely explicit asymmetric errors
    // Scatter(const std::vector<double>& x, const std::vector<double>& y,
    //         const std::vector<double>& exminus,
    //         const std::vector<double>& explus,
    //         const std::vector<double>& eyminus,
    //         const std::vector<double>& eyplus,
    //         const std::string& path="", const std::string& title="")
    //   : AnalysisObject("Scatter", path, title)
    // {
    //   assert(x.size() == y.size() &&
    //          x.size() == exminus.size() && x.size() == explus.size() &&
    //          x.size() == eyminus.size() && x.size() == eyplus.size());
    //   for (size_t i = 0; i < x.size(); ++i) {
    //     addPoint(Point<N>2D(x[i], exminus[i], explus[i], y[i], eyminus[i], eyplus[i]));
    //   }
    // }


    /// Copy constructor with optional new path
    Scatter(const Scatter<N>& s, const std::string& path="")
      : AnalysisObject("Scatter", (path.size() == 0) ? s.path() : path, s, s.title()),
        _points(s._points)
    {  }


    /// Assignment operator
    Scatter<N>& operator = (const Scatter<N>& s) {
      setPath(s.path());
      setTitle(s.title());
      _points = s._points;
      return *this;
    }

    /// Make a copy on the stack
    Scatter<N> clone() const {
      return Scatter<N>(*this);
    }

    /// Make a copy on the heap, via 'new'
    Scatter<N>* newclone() const {
      return new Scatter<N>(*this);
    }

    /// @}


    /// @name Modifiers
    /// @{

    /// Clear all points
    void reset() {
      _points.clear();
    }

    /// Scaling
    void scale(const NdVal& scales) {
      for (Point<N>& p : _points) p.scale(scales);
    }

    /// @}


    ///////////////////////////////////////////////////


    /// @name Point accessors
    /// @{

    /// Number of points in the scatter
    size_t numPoints() const {
      return _points.size();
    }


    /// Get the collection of points
    Points& points() {
      return _points;
    }


    /// Get the collection of points (const version)
    const Points& points() const {
      return _points;
    }


    /// Get a reference to the point with index @a index
    Point<N>& point(size_t index) {
      return _points.at(index);
    }


    /// Get the point with index @a index (const version)
    const Point<N>& point(size_t index) const {
      return _points.at(index);
    }

    /// @}


    /// @name Point inserters
    /// @{

    /// Insert a new point
    Scatter<N>& addPoint(const Point<N>& pt) {
      _points.insert(pt);
      return *this;
    }

    /// Insert a new point, from a position array
    Scatter<N>& addPoint(const NdVal& pos) {
      _points.insert(Point<N>(pos));
      return *this;
    }

    /// @todo More addPoint combinations with arrays for errors

    /// Insert a collection of new points
    Scatter<N>& addPoints(Points pts) {
      for (const Point<N>& pt : pts) addPoint(pt);
      return *this;
    }

    /// @}


    /// @name Combining sets of scatter points
    /// @{

    /// @todo Better name?
    Scatter<N>& combineWith(const Scatter<N>& other) {
      addPoints(other.points());
      return *this;
    }

    /// @todo Better name?
    Scatter<N>& combineWith(const std::vector< Scatter<N> >& others) {
      for (const Scatter<N>& s : others) combineWith(s);
      return *this;
    }

    /// @}


  private:

    Points _points;

  };


  /// @name Combining scatters by merging sets of points
  /// @{

  template <int N>
  inline Scatter<N> combine(const Scatter<N>& a, const Scatter<N>& b) {
    Scatter<N> rtn = a;
    rtn.combineWith(b);
    return rtn;
  }

  template <int N>
  inline Scatter<N> combine(const std::vector< Scatter<N> >& scatters) {
    Scatter<N> rtn;
    rtn.combineWith(scatters);
    return rtn;
  }

  /// @}


  //////////////////////////////////


  // /// @name Combining scatters: global operators, assuming aligned points
  // /// @{

  // /// Add two scatters
  // template <int N>
  // Scatter add(const Scatter& first, const Scatter& second);


  // /// Add two scatters
  // template <int N>
  // inline Scatter operator + (const Scatter& first, const Scatter& second) {
  //   return add(first, second);
  // }


  // /// Subtract two scatters
  // template <int N>
  // Scatter subtract(const Scatter& first, const Scatter& second);


  // /// Subtract two scatters
  // template <int N>
  // inline Scatter operator - (const Scatter& first, const Scatter& second) {
  //   return subtract(first, second);
  // }


  // /// Divide two scatters
  // template <int N>
  // Scatter divide(const Scatter& numer, const Scatter& denom);


  // /// Divide two scatters
  // template <int N>
  // inline Scatter operator / (const Scatter& numer, const Scatter& denom) {
  //   return divide(numer, denom);
  // }

  // /// @}


}

#endif
