// -*- C++ -*-
//
// This file is part of YODA -- Yet more Objects for Data Analysis
// Copyright (C) 2008-2021 The YODA collaboration (see AUTHORS for details)
//
#ifndef YODA_Scatter3D_h
#define YODA_Scatter3D_h

#include "YODA/AnalysisObject.h"
#include "YODA/Scatter.h"
#include "YODA/Point3D.h"
#include "YODA/Utils/sortedvector.h"
#include <utility>
#include <memory>

namespace YODA {


  // Forward declarations
  class Histo2D;
  class Profile2D;


  /// A very generic data type which is just a collection of 3D data points with errors
  class Scatter3D : public AnalysisObject, public Scatter {
  public:

    /// Types of the native Point3D collection
    typedef Point3D Point;
    typedef Utils::sortedvector<Point3D> Points;
    typedef std::shared_ptr<Scatter3D> Ptr;


    /// @name Constructors
    /// @{

    /// Empty constructor
    Scatter3D(const std::string& path="", const std::string& title="")
      : AnalysisObject("Scatter3D", path, title)
    { }


    /// Constructor from a set of points
    Scatter3D(const Points& points,
              const std::string& path="", const std::string& title="")
      : AnalysisObject("Scatter3D", path, title),
        _points(points)
    {
      std::sort(_points.begin(), _points.end());
    }


    /// Constructor from vectors of values with no errors
    Scatter3D(const std::vector<double>& x,
        const std::vector<double>& y,
        const std::vector<double>& z,
              const std::string& path="",
        const std::string& title="")
      : AnalysisObject("Scatter3D", path, title)
    {
      if (x.size() != y.size() || y.size() != z.size()) {
        throw RangeError("There are different numbers of x, y, and z values in the provided vectors.");
      }
      const std::pair<double,double> nullerr = std::make_pair(0.0, 0.0);
      for (size_t i = 0; i < x.size(); ++i) {
        addPoint(Point3D(x[i], y[i], z[i], nullerr, nullerr, nullerr));
      }
      std::sort(_points.begin(), _points.end());
    }


    /// Constructor from vectors of values with asymmetric errors on both x and y
    Scatter3D(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& z,
              const std::vector<std::pair<double,double> >& ex, const std::vector<std::pair<double,double> >& ey, const std::vector<std::pair<double,double> >& ez,
              const std::string& path="", const std::string& title="")
      : AnalysisObject("Scatter3D", path, title)
    {
      if (x.size() != y.size() || y.size() != z.size()) {
        throw RangeError("There are different numbers of x, y, and z values in the provided vectors.");
      }
      if (x.size() != ex.size() || y.size() != ey.size() || z.size() != ez.size()) {
        throw RangeError("The sizes of the provided error vectors don't match the corresponding x, y, or z value vectors.");
      }
      for (size_t i = 0; i < x.size(); ++i) {
        addPoint(Point3D(x[i], y[i], z[i], ex[i], ey[i], ez[i]));
      }
      std::sort(_points.begin(), _points.end());
    }


    /// Constructor from vectors of values with completely explicit asymmetric errors
    Scatter3D(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double> z,
              const std::vector<double>& exminus,
              const std::vector<double>& explus,
              const std::vector<double>& eyminus,
              const std::vector<double>& eyplus,
              const std::vector<double>& ezminus,
              const std::vector<double>& ezplus,
              const std::string& path="", const std::string& title="")
      : AnalysisObject("Scatter3D", path, title)
    {
      if(x.size() != y.size() || y.size() != z.size() ||
         x.size() != exminus.size() || x.size() != explus.size() ||
         y.size() != eyminus.size() || y.size() != eyplus.size() ||
         z.size() != ezminus.size() || z.size() != ezplus.size())
        throw RangeError("There are either different amounts of points on x/y/z vectors or not every of these vectors has properly defined error vectors!");

      for (size_t i = 0; i < x.size(); ++i) {
        addPoint(Point3D(x[i], y[i], z[i], exminus[i], explus[i], eyminus[i], eyplus[i], ezminus[i], ezplus[i]));
      }

      std::sort(_points.begin(), _points.end());
    }


    /// Copy constructor with optional new path
    /// @todo Also allow title setting from the constructor?
    Scatter3D(const Scatter3D& s3, const std::string& path="")
      : AnalysisObject("Scatter3D", (path.size() == 0) ? s3.path() : path, s3, s3.title()),
        _points(s3._points)
    {
      for ( auto &ann : annotations()) setAnnotation(ann, annotation(ann));
      for ( auto &pt : _points) pt.setParent(this);
    }

    /// Assignment operator
    Scatter3D& operator = (const Scatter3D& s3) {
      AnalysisObject::operator = (s3); //< AO treatment of paths etc.
      _points = s3._points;
      return *this;
    }

    /// Make a copy on the stack
    Scatter3D clone() const {
      return Scatter3D(*this);
    }

    /// Make a copy on the heap, via 'new'
    Scatter3D* newclone() const {
      return new Scatter3D(*this);
    }

    /// @}


    /// Dimension of this data object
    size_t dim() const { return 3; }


    /// @name Modifiers
    /// @{

    /// Clear all points
    void reset() {
      _points.clear();
    }

    /// Scaling of x axis
    void scaleX(double scalex) {
      for (Point3D& p : _points) p.scaleX(scalex);
    }

    /// Scaling of y axis
    void scaleY(double scaley) {
      for (Point3D& p : _points) p.scaleY(scaley);
    }

    /// Scaling of z axis
    void scaleZ(double scalez) {
      for (Point3D& p : _points) p.scaleZ(scalez);
    }

    /// Scaling of all three axes
    void scaleXYZ(double scalex, double scaley, double scalez) {
      for (Point3D& p : _points) p.scaleXYZ(scalex, scaley, scalez);
    }

    /// Scaling along direction @a i
    void scale(size_t i, double scale) {
      switch (i) {
      case 1: scaleX(scale); break;
      case 2: scaleY(scale); break;
      case 3: scaleZ(scale); break;
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


    /// Get a reference to the point with index @a index
    Point3D& point(size_t index) {
      if (index >= numPoints()) throw RangeError("There is no point with this index");
      return _points.at(index);
    }


    /// Get the point with index @a index (const version)
    const Point3D& point(size_t index) const {
      if (index >= numPoints()) throw RangeError("There is no point with such index!");
      return _points.at(index);
    }

    /// @}


    /// @name Point inserters
    /// @{

    /// Insert a new point and assign
    /// this scatter as its parent
    void addPoint(Point3D pt) {
      pt.setParent(this);
      _points.insert(pt);
    }

    /// Insert a new point, defined as the x/y/z value triplet and no errors
    void addPoint(double x, double y, double z) {
      Point3D thisPoint = Point3D(x, y, z);
      thisPoint.setParent(this);
      _points.insert(thisPoint);
    }

    /// Insert a new point, defined as the x/y/z value triplet and symmetric errors
    void addPoint(double x, double y, double z,
                  double ex, double ey, double ez) {
      Point3D thisPoint = Point3D(x, y, z, ex, ey, ez);
      thisPoint.setParent(this);
      _points.insert(thisPoint);
    }

    /// Insert a new point, defined as the x/y/z value triplet and asymmetric error pairs
    void addPoint(double x, double y, double z,
                  const std::pair<double,double>& ex, const std::pair<double,double>& ey, const std::pair<double,double>& ez) {
      Point3D thisPoint = Point3D(x, y, z, ex, ey, ez);
      thisPoint.setParent(this);
      _points.insert(thisPoint);
    }

    /// Insert a new point, defined as the x/y/z value triplet and asymmetric errors
    void addPoint(double x, double y, double z,
                  double exminus, double explus,
                  double eyminus, double eyplus,
                  double ezminus, double ezplus) {
      Point3D thisPoint = Point3D(x, y, z, exminus, explus, eyminus, eyplus, ezminus, ezplus);
      thisPoint.setParent(this);
      _points.insert(thisPoint);
    }

    /// Insert a collection of new points
    void addPoints(const Points& pts) {
      for (const Point3D& pt : pts) addPoint(pt);
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


    /// @todo Better name?
    void combineWith(const Scatter3D& other) {
      addPoints(other.points());
      //return *this;
    }


    /// @todo Better name?
    /// @todo Convert to accept a Range or generic
    void combineWith(const std::vector<Scatter3D>& others) {
      for (const Scatter3D& s : others) combineWith(s);
    }


    /// Equality operator
    bool operator == (const Scatter3D& other) {
      return _points == other._points;
    }

    /// Non-equality operator
    bool operator != (const Scatter3D& other) {
      return ! operator == (other);
    }


    //////////////////////////////////



  private:

    Points _points;

    bool _variationsParsed =false ;

  };


  /// Convenience typedef
  typedef Scatter3D S3D;


  /// @name Combining scatters by merging sets of points
  /// @{

  inline Scatter3D combine(const Scatter3D& a, const Scatter3D& b) {
    Scatter3D rtn = a;
    rtn.combineWith(b);
    return rtn;
  }

  inline Scatter3D combine(const std::vector<Scatter3D>& scatters) {
    Scatter3D rtn;
    rtn.combineWith(scatters);
    return rtn;
  }

  /// @}


  //////////////////////////////////


  /// @name Conversion functions from other data types
  /// @{

  /// Make a Scatter3D representation of a Histo2D
  ///
  /// Optional @c usefocus argument can be used to position the point at the bin
  /// focus rather than geometric midpoint.
  Scatter3D mkScatter(const Histo2D& h, bool usefocus=false, bool binareadiv=true);

  /// Make a Scatter3D representation of a Profile2D
  ///
  /// Optional @c usefocus argument can be used to position the point at the bin
  /// focus rather than geometric midpoint. Optional @c usestddev argument can
  /// be used to draw the distribution sigma rather than the standard error on
  /// the mean as the z-error bar size.
  Scatter3D mkScatter(const Profile2D& p, bool usefocus=false, bool usestddev=false);

  /// Make a Scatter3D representation of... erm, a Scatter3D!
  /// @note Mainly exists to allow mkScatter to be called on any AnalysisObject type
  inline Scatter3D mkScatter(const Scatter3D& s) { return Scatter3D(s); }
  // /// @note The usefocus arg is just for consistency and has no effect for Scatter -> Scatter
  //inline Scatter3D mkScatter(const Scatter3D& s, bool) { return mkScatter(s); }

  /// @}


  /////////////////////////////////


  /// @name Transforming operations on Scatter3D
  /// @{

  /// @brief Apply transformation fx(x) to all values and error positions (operates in-place on @a s)
  ///
  /// fx should be a function which takes double x -> double newx
  template<typename FNX>
  inline void transformX(Scatter3D& s, FNX fx) {
    for (size_t i = 0; i < s.numPoints(); ++i) {
      Point3D& p = s.point(i);
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
  inline void transformY(Scatter3D& s, FNY fy) {
    for (size_t i = 0; i < s.numPoints(); ++i) {
      Point3D& p = s.point(i);
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


  /// @brief Apply transformation fz(z) to all values and error positions (operates in-place on @a s)
  ///
  /// fz should be a function which takes double z -> double newz
  template<typename FNZ>
  inline void transformZ(Scatter3D& s, FNZ fz) {
    for (size_t i = 0; i < s.numPoints(); ++i) {
      Point3D& p = s.point(i);
      const double newz = fz(p.z());
      const double fz_zmin = fz(p.zMin());
      const double fz_zmax = fz(p.zMax());
      // Deal with possible inversions of min/max ordering under the transformation
      const double newzmin = std::min(fz_zmin, fz_zmax);
      const double newzmax = std::max(fz_zmin, fz_zmax);
      // Set new point z values
      p.setZ(newz);
      /// @todo Be careful about transforms which could switch around min and max errors, or send both in the same direction!
      p.setZErrMinus(newz - newzmin);
      p.setZErrPlus(newzmax - newz);
    }
  }


  /// @todo Add external scale, scaleX, scaleY, scaleZ functions

  /// @}



}

#endif
