// -*- C++ -*-
//
// This file is part of YODA -- Yet more Objects for Data Analysis
// Copyright (C) 2008-2021 The YODA collaboration (see AUTHORS for details)
//
#ifndef YODA_Scatter1D_h
#define YODA_Scatter1D_h

#include "YODA/AnalysisObject.h"
#include "YODA/Scatter.h"
#include "YODA/Point1D.h"
#include "YODA/Utils/sortedvector.h"
#include <utility>
#include <memory>

namespace YODA {


  // Forward declarations
  class Counter;


  /// A very generic data type which is just a collection of 1D data points with errors
  ///
  /// @todo Add a generic Scatter base class, providing reset(), rmPoint(), etc.
  class Scatter1D : public AnalysisObject, public Scatter {
  public:

    /// Type of the native Point1D collection
    typedef Point1D Point;
    typedef Utils::sortedvector<Point1D> Points;
    typedef std::shared_ptr<Scatter1D> Ptr;


    /// @name Constructors
    /// @{

    /// Empty constructor
    Scatter1D(const std::string& path="", const std::string& title="")
      : AnalysisObject("Scatter1D", path, title)
    {  }


    /// Constructor from a set of points
    Scatter1D(const Points& points,
              const std::string& path="", const std::string& title="")
      : AnalysisObject("Scatter1D", path, title),
        _points(points)
    {  }


    /// Constructor from a vector of x values with no errors
    Scatter1D(const std::vector<double>& x,
              const std::string& path="", const std::string& title="")
      : AnalysisObject("Scatter1D", path, title)
    {
      for (size_t i = 0; i < x.size(); ++i) addPoint(x[i]);
    }


    /// Constructor from vectors of x values with symmetric errors
    Scatter1D(const std::vector<double>& x, const std::vector<double>& ex,
              const std::string& path="", const std::string& title="")
      : AnalysisObject("Scatter1D", path, title)
    {
      if (x.size() != ex.size()) throw UserError("x and ex vectors must have same length");
      for (size_t i = 0; i < x.size(); ++i) addPoint(x[i], ex[i]);
    }

    /// Constructor from x values with asymmetric errors
    Scatter1D(const std::vector<double>& x, const std::vector<std::pair<double,double> >& ex,
              const std::string& path="", const std::string& title="")
      : AnalysisObject("Scatter1D", path, title)
    {
      if (x.size() != ex.size()) throw UserError("x and ex vectors must have same length");
      for (size_t i = 0; i < x.size(); ++i) addPoint(Point1D(x[i], ex[i]));
    }


    /// Constructor from values with completely explicit asymmetric errors
    Scatter1D(const std::vector<double>& x,
              const std::vector<double>& exminus,
              const std::vector<double>& explus,
              const std::string& path="", const std::string& title="")
      : AnalysisObject("Scatter1D", path, title)
    {
      if (x.size() != exminus.size()) throw UserError("x and ex vectors must have same length");
      if (exminus.size() != explus.size()) throw UserError("ex plus and minus vectors must have same length");
      for (size_t i = 0; i < x.size(); ++i) addPoint(Point1D(x[i], exminus[i], explus[i]));
    }


    /// Copy constructor with optional new path
    /// @todo Also allow title setting from the constructor?
    Scatter1D(const Scatter1D& s1, const std::string& path="")
      : AnalysisObject("Scatter1D", (path.size() == 0) ? s1.path() : path, s1, s1.title()),
        _points(s1._points)
    {
      for ( auto &ann : annotations()) setAnnotation(ann, annotation(ann));
      for ( auto &pt : _points) pt.setParent(this);;
      
    }


    /// Assignment operator
    Scatter1D& operator = (const Scatter1D& s1) {
      AnalysisObject::operator = (s1); //< AO treatment of paths etc.
      _points = s1._points;
      return *this;
    }

    /// Make a copy on the stack
    Scatter1D clone() const {
      return Scatter1D(*this);
    }

    /// Make a copy on the heap, via 'new'
    Scatter1D* newclone() const {
      return new Scatter1D(*this);
    }

    /// @}


    /// Dimension of this data object
    size_t dim() const { return 1; }


    /// @name Modifiers
    /// @{

    /// Clear all points
    void reset() {
      _points.clear();
    }

    /// Scaling of x axis
    void scaleX(double scalex) {
      for (Point1D& p : _points) p.scaleX(scalex);
    }

    /// Scaling along direction @a i
    void scale(size_t i, double scale) {
      switch (i) {
      case 1: scaleX(scale); break;
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
    Point1D& point(size_t index) {
      if (index >= numPoints()) throw RangeError("There is no point with this index");
      return _points.at(index);
    }


    /// Get a reference to the point with index @a index (const)
    const Point1D& point(size_t index) const {
      if (index >= numPoints()) throw RangeError("There is no point with this index");
      return _points.at(index);
    }

    /// @}


    /// @name Point inserters
    /// @{

    // Insert a new point and assign 
    /// this scatter as its parent
    void addPoint(Point1D pt) {
      pt.setParent(this);
      _points.insert(pt);
    }

    /// Insert a new point, defined as the x value and no errors
    void addPoint(double x) {
      Point1D thisPoint = Point1D(x);
      thisPoint.setParent(this);
      _points.insert(thisPoint);
    }

    /// Insert a new point, defined as the x value and symmetric errors
    void addPoint(double x, double ex) {
      Point1D thisPoint = Point1D(x, ex);
      thisPoint.setParent(this);
      _points.insert(thisPoint);
    }

    /// Insert a new point, defined as the x value and an asymmetric error pair
    void addPoint(double x, const std::pair<double,double>& ex) {
      Point1D thisPoint = Point1D(x, ex);
      thisPoint.setParent(this);
      _points.insert(thisPoint);
    }

    /// Insert a new point, defined as the x value and explicit asymmetric errors
    void addPoint(double x, double exminus, double explus) {
      Point1D thisPoint = Point1D(x, exminus, explus);
      thisPoint.setParent(this);
      _points.insert(thisPoint);
    }

    /// Insert a collection of new points
    void addPoints(const Points& pts) {
      for (const Point1D& pt : pts) addPoint(pt);
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

    /// @todo Better name?
    void combineWith(const Scatter1D& other) {
      addPoints(other.points());
    }

    /// @todo Better name? Make this the add operation?
    /// @todo Convert/extend to accept a Range or generic
    void combineWith(const std::vector<Scatter1D>& others) {
      for (const Scatter1D& s : others) combineWith(s);
    }

    /// @}


    /// Equality operator
    bool operator == (const Scatter1D& other) {
      return _points == other._points;
    }

    /// Non-equality operator
    bool operator != (const Scatter1D& other) {
      return ! operator == (other);
    }


    //////////////////////////////////



  private:

    Points _points;

    bool _variationsParsed = false;

  };


  /// Convenience typedef
  typedef Scatter1D S1D;


  /// @name Combining scatters by merging sets of points
  /// @{

  inline Scatter1D combine(const Scatter1D& a, const Scatter1D& b) {
    Scatter1D rtn = a;
    rtn.combineWith(b);
    return rtn;
  }

  inline Scatter1D combine(const std::vector<Scatter1D>& scatters) {
    Scatter1D rtn;
    rtn.combineWith(scatters);
    return rtn;
  }

  /// @}

  //////////////////////////////////


  /// @name Conversion functions from other data types
  /// @{

  /// Make a Scatter1D representation of a Histo1D
  Scatter1D mkScatter(const Counter& c);

  /// Make a Scatter1D representation of... erm, a Scatter1D!
  /// @note Mainly exists to allow mkScatter to be called on any AnalysisObject type
  inline Scatter1D mkScatter(const Scatter1D& s) {
    return Scatter1D(s);
  }

  /// @}


  /////////////////////////////////


  /// @name Transforming operations on Scatter1D
  /// @{

  /// @brief Apply transformation fx(x) to all values and error positions (operates in-place on @a s)
  ///
  /// fx should be a function which takes double x -> double newx
  template<typename FNX>
  inline void transformX(Scatter1D& s, FNX fx) {
    for (size_t i = 0; i < s.numPoints(); ++i) {
      Point1D& p = s.point(i);
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

  /// @}


}

#endif
