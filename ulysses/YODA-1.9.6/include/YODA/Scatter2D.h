// -*- C++ -*-
//
// This file is part of YODA -- Yet more Objects for Data Analysis
// Copyright (C) 2008-2021 The YODA collaboration (see AUTHORS for details)
//
#ifndef YODA_Scatter2D_h
#define YODA_Scatter2D_h

#include "YODA/AnalysisObject.h"
#include "YODA/Scatter.h"
#include "YODA/Point2D.h"
#include "YODA/Utils/sortedvector.h"
#include <utility>
#include <memory>

namespace YODA {


  // Forward declarations
  class Histo1D;
  class Profile1D;


  /// A very generic data type which is just a collection of 2D data points with errors
  class Scatter2D : public AnalysisObject, public Scatter {
  public:

    /// Type of the native Point2D collection
    typedef Point2D Point;
    typedef Utils::sortedvector<Point2D> Points;
    typedef std::shared_ptr<Scatter2D> Ptr;


    /// @name Constructors
    /// @{

    /// Empty constructor
    Scatter2D(const std::string& path="", const std::string& title="")
      : AnalysisObject("Scatter2D", path, title)
    {  }


    /// Constructor from a set of points
    Scatter2D(const Points& points,
              const std::string& path="", const std::string& title="")
      : AnalysisObject("Scatter2D", path, title),
        _points(points)
    {  }


    /// Constructor from a vector of values with no errors
    Scatter2D(const std::vector<double>& x, const std::vector<double>& y,
              const std::string& path="", const std::string& title="")
      : AnalysisObject("Scatter2D", path, title)
    {
      if (x.size() != y.size()) throw UserError("x and y vectors must have same length");
      for (size_t i = 0; i < x.size(); ++i) addPoint(x[i], y[i]);
    }


    /// Constructor from vectors of values with symmetric errors on x and y
    Scatter2D(const std::vector<double>& x, const std::vector<double>& y,
              const std::vector<double>& ex, const std::vector<double>& ey,
              const std::string& path="", const std::string& title="")
      : AnalysisObject("Scatter2D", path, title)
    {
      if (x.size() != y.size()) throw UserError("x and y vectors must have same length");
      if (x.size() != ex.size()) throw UserError("x and ex vectors must have same length");
      if (y.size() != ey.size()) throw UserError("y and ey vectors must have same length");
      for (size_t i = 0; i < x.size(); ++i) addPoint(x[i], y[i], ex[i], ey[i]);
    }


    /// Constructor from values with asymmetric errors on both x and y
    Scatter2D(const std::vector<double>& x, const std::vector<double>& y,
              const std::vector<std::pair<double,double> >& ex, const std::vector<std::pair<double,double> >& ey,
              const std::string& path="", const std::string& title="")
      : AnalysisObject("Scatter2D", path, title)
    {
      if (x.size() != y.size()) throw UserError("x and y vectors must have same length");
      if (x.size() != ex.size()) throw UserError("x and ex vectors must have same length");
      if (y.size() != ey.size()) throw UserError("y and ey vectors must have same length");
      for (size_t i = 0; i < x.size(); ++i) addPoint(Point2D(x[i], y[i], ex[i], ey[i]));
    }


    /// Constructor from values with completely explicit asymmetric errors
    Scatter2D(const std::vector<double>& x, const std::vector<double>& y,
              const std::vector<double>& exminus, const std::vector<double>& explus,
              const std::vector<double>& eyminus, const std::vector<double>& eyplus,
              const std::string& path="", const std::string& title="")
      : AnalysisObject("Scatter2D", path, title)
    {
      if (x.size() != y.size()) throw UserError("x and y vectors must have same length");
      if (x.size() != exminus.size()) throw UserError("x and ex vectors must have same length");
      if (y.size() != eyminus.size()) throw UserError("y and ey vectors must have same length");
      if (exminus.size() != explus.size()) throw UserError("ex plus and minus vectors must have same length");
      if (eyminus.size() != eyplus.size()) throw UserError("ey plus and minus vectors must have same length");
      for (size_t i = 0; i < x.size(); ++i) addPoint(Point2D(x[i], y[i], exminus[i], explus[i], eyminus[i], eyplus[i]));
    }


    /// Copy constructor with optional new path
    /// @todo Also allow title setting from the constructor?
    Scatter2D(const Scatter2D& s2, const std::string& path="")
      : AnalysisObject("Scatter2D", (path.size() == 0) ? s2.path() : path, s2, s2.title()),
        _points(s2._points)
    {
      for ( auto &ann : annotations()) setAnnotation(ann, annotation(ann));
      for ( auto &pt : _points) pt.setParent(this);
    }


    /// Assignment operator
    Scatter2D& operator = (const Scatter2D& s2) {
      AnalysisObject::operator = (s2); //< AO treatment of paths etc.
      _points = s2._points;
      return *this;
    }

    /// Make a copy on the stack
    Scatter2D clone() const {
      return Scatter2D(*this);
    }

    /// Make a copy on the heap, via 'new'
    Scatter2D* newclone() const {
      return new Scatter2D(*this);
    }

    /// @}


    /// Dimension of this data object
    size_t dim() const { return 2; }


    /// @name Modifiers
    /// @{

    /// Clear all points
    void reset() {
      _points.clear();
    }

    /// Scaling of x axis
    void scaleX(double scalex) {
      for (Point2D& p : _points) p.scaleX(scalex);
    }

    /// Scaling of y axis
    void scaleY(double scaley) {
      for (Point2D& p : _points) p.scaleY(scaley);
    }

    /// Scaling of both axes
    void scaleXY(double scalex, double scaley) {
      for (Point2D& p : _points) p.scaleXY(scalex, scaley);
    }

    /// Scaling along direction @a i
    void scale(size_t i, double scale) {
      switch (i) {
      case 1: scaleX(scale); break;
      case 2: scaleY(scale); break;
      default: throw RangeError("Invalid axis int, must be in range 1..dim");
      }
    }

    /// @}


    ///////////////////////////////////////////////////

    /// read/write the variations stored in the points as annotations
    void parseVariations();
    void writeVariationsToAnnotations();
    void updateTotalUncertainty();

    /// Get the list of variations stored in the points
    std::vector<std::string> variations() const;

    /// Remove the variations
    void rmVariations();

    // Construct a covariance matrix from the error breakdown
    std::vector<std::vector<double> > covarianceMatrix(bool ignoreOffDiagonalTerms=false);

    /// @name Point accessors
    /// @{

    /// Number of points in the scatter
    size_t numPoints() const {
      return _points.size();
    }


    /// Get the collection of points (non-const)
    Points& points() {
      return _points;
    }


    /// Get the collection of points (const)
    const Points& points() const {
      return _points;
    }


    /// Get a reference to the point with index @a index (non-const)
    Point2D& point(size_t index) {
      if (index >= numPoints()) throw RangeError("There is no point with this index");
      return _points.at(index);
    }


    /// Get a reference to the point with index @a index (const)
    const Point2D& point(size_t index) const {
      if (index >= numPoints()) throw RangeError("There is no point with this index");
      return _points.at(index);
    }

    /// @}


    /// @name Point inserters
    /// @{

    // Insert a new point and assign
    /// this scatter as its parent
    void addPoint(Point2D pt) {
      pt.setParent(this);
      _points.insert(pt);
    }

    /// Insert a new point, defined as the x/y value pair and no errors
    void addPoint(double x, double y) {
      Point2D thisPoint = Point2D(x, y);
      thisPoint.setParent(this);
      _points.insert(thisPoint);
    }

    /// Insert a new point, defined as the x/y value pair and symmetric errors
    void addPoint(double x, double y,
                  double ex, double ey) {
      Point2D thisPoint = Point2D(x, y, ex, ey);
      thisPoint.setParent(this);
      _points.insert(thisPoint);
    }

    /// Insert a new point, defined as the x/y value pair and asymmetric error pairs
    void addPoint(double x, double y,
                  const std::pair<double,double>& ex, const std::pair<double,double>& ey) {
      Point2D thisPoint = Point2D(x, y, ex, ey);
      thisPoint.setParent(this);
      _points.insert(thisPoint);
    }

    /// Insert a new point, defined as the x/y value pair and asymmetric errors
    void addPoint(double x, double y,
                  double exminus, double explus,
                  double eyminus, double eyplus) {
      Point2D thisPoint = Point2D(x, y, exminus, explus, eyminus, eyplus);
      thisPoint.setParent(this);
      _points.insert(thisPoint);
    }

    /// Insert a collection of new points
    void addPoints(const Points& pts) {
      for (const Point2D& pt : pts) addPoint(pt);
    }

    /// @}


    /// @name Point removers
    /// @{

    /// Remove the point with index @a index
    void rmPoint(size_t index) {
      _points.erase(_points.begin()+index);
    }

    // /// Remove the points with indices @a indices
    // void rmPoints(std::vector<size_t> indices) {
    //   // reverse-sort so the erasure-loop doesn't invalidate the indices
    //   std::sort(indices.begin(), indices.end(), std::greater<size_t>());
    //   for (size_t i : indices) rmPoint(i);
    // }

    /// @}



    /// @name Combining sets of scatter points
    /// @{

    /// @todo Better name? Make this the add operation?
    void combineWith(const Scatter2D& other) {
      addPoints(other.points());
      //return *this;
    }

    /// @todo Better name?
    /// @todo Convert/extend to accept a Range or generic
    void combineWith(const std::vector<Scatter2D>& others) {
      for (const Scatter2D& s : others) combineWith(s);
      //return *this;
    }

    /// @}


    /// Equality operator
    bool operator == (const Scatter2D& other) {
      return _points == other._points;
    }

    /// Non-equality operator
    bool operator != (const Scatter2D& other) {
      return ! operator == (other);
    }



    //////////////////////////////////


  private:

    Points _points;

    bool _variationsParsed =false ;

  };


  /// Convenience typedef
  typedef Scatter2D S2D;


  /// @name Combining scatters by merging sets of points
  /// @{

  inline Scatter2D combine(const Scatter2D& a, const Scatter2D& b) {
    Scatter2D rtn = a;
    rtn.combineWith(b);
    return rtn;
  }

  inline Scatter2D combine(const std::vector<Scatter2D>& scatters) {
    Scatter2D rtn;
    rtn.combineWith(scatters);
    return rtn;
  }

  /// @}


  //////////////////////////////////


  /// @name Conversion functions from other data types
  /// @{

  /// Make a Scatter2D representation of a Histo1D
  ///
  /// Optional @c usefocus argument can be used to position the point at the bin
  /// focus rather than geometric midpoint. Optional @c binwidthdiv argument can be
  /// used to disable the default (physical, differential!) scaling of y values and
  /// errors by 1/bin-width.
  Scatter2D mkScatter(const Histo1D& h, bool usefocus=false, bool binwidthdiv=true,
                      double uflow_binwidth=-1, double oflow_binwidth=-1);

  /// Make a Scatter2D representation of a Profile1D
  ///
  /// Optional @c usefocus argument can be used to position the point at the bin
  /// focus rather than geometric midpoint. Optional @c usestddev argument can
  /// be used to draw the y-distribution sigma rather than the standard error on
  /// the mean as the y-error bar size.
  Scatter2D mkScatter(const Profile1D& p, bool usefocus=false, bool usestddev=false,
                      double uflow_binwidth=-1, double oflow_binwidth=-1);

  /// Make a Scatter2D representation of... erm, a Scatter2D!
  ///
  /// @note Mainly exists to allow (no-opt-args) mkScatter to be called on any AnalysisObject type
  inline Scatter2D mkScatter(const Scatter2D& s) { return Scatter2D(s); }
  // /// @note The usefocus arg is just for consistency and has no effect for Scatter -> Scatter
  // inline Scatter2D mkScatter(const Scatter2D& s, bool) { return mkScatter(s); }

  /// @}


  //////////////////////////////////


  /// @name Transforming operations on Scatter2D
  /// @{

  /// @brief Apply transformation fx(x) to all values and error positions (operates in-place on @a s)
  ///
  /// fx should be a function which takes double x -> double newx
  template<typename FNX>
  inline void transformX(Scatter2D& s, FNX fx) {
    for (size_t i = 0; i < s.numPoints(); ++i) {
      Point2D& p = s.point(i);
      const double newx = fx(p.x());
      const double fx_xmin = fx(p.xMin());
      const double fx_xmax = fx(p.xMax());
      // Deal with possible inversions of min/max ordering under the transformation
      const double newxmin = std::min(fx_xmin, fx_xmax);
      const double newxmax = std::max(fx_xmin, fx_xmax);
      // Set new point x values
      p.setX(newx);
      /// @todo Be careful about transforms which could switch around min and max errors, or send both in the same direction!
      p.setXErrMinus(newx - newxmin);
      p.setXErrPlus(newxmax - newx);
    }
  }


  /// @brief Apply transformation fy(y) to all values and error positions (operates in-place on @a s)
  ///
  /// fy should be a function which takes double y -> double newy
  template<typename FNY>
  inline void transformY(Scatter2D& s, FNY fy) {
    for (size_t i = 0; i < s.numPoints(); ++i) {
      Point2D& p = s.point(i);
      const double newy = fy(p.y());
      const double fy_ymin = fy(p.yMin());
      const double fy_ymax = fy(p.yMax());
      // Deal with possible inversions of min/max ordering under the transformation
      const double newymin = std::min(fy_ymin, fy_ymax);
      const double newymax = std::max(fy_ymin, fy_ymax);
      // Set new point y values
      p.setY(newy);
      /// @todo Be careful about transforms which could switch around min and max errors, or send both in the same direction!
      p.setYErrMinus(newy - newymin);
      p.setYErrPlus(newymax - newy);
    }
  }


  /// @todo Add external scale, scaleX, scaleY functions


  /// Exchange the x and y axes (operates in-place on @a s)
  inline void flip(Scatter2D& s) {
    for (size_t i = 0; i < s.numPoints(); ++i) {
      Point2D& p = s.point(i);
      const double newx = p.y();
      const double newy = p.x();
      const double newxmin = p.yMin();
      const double newxmax = p.yMax();
      const double newymin = p.xMin();
      const double newymax = p.xMax();
      p.setX(newx);
      p.setY(newy);
      /// @todo Be careful about transforms which could switch around min and max errors, or send both in the same direction!
      p.setXErrMinus(newx - newxmin);
      p.setXErrPlus(newxmax - newx);
      p.setYErrMinus(newy - newymin);
      p.setYErrPlus(newymax - newy);
    }
  }

  /// @}


}

#endif
