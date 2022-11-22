// -*- C++ -*-
//
// This file is part of YODA -- Yet more Objects for Data Analysis
// Copyright (C) 2008-2021 The YODA collaboration (see AUTHORS for details)
//
#ifndef YODA_POINT2D_H
#define YODA_POINT2D_H

#include "YODA/Point.h"
#include "YODA/Exceptions.h"
#include "YODA/Utils/MathUtils.h"
#include <utility>


namespace YODA {


  /// A 2D data point to be contained in a Scatter2D
  class Point2D : public Point {
  public:

    /// @name Constructors
    /// @{

    // Default constructor
    Point2D() {  }


    /// Constructor from values with optional symmetric errors
    Point2D(double x, double y, double ex=0.0, double ey=0.0, std::string source="")
      : _x(x), _y(y)
    {
      _ex = std::make_pair(ex, ex);
      _ey[source] = std::make_pair(ey, ey);
    }


    /// Constructor from values with explicit asymmetric errors
    Point2D(double x, double y,
            double exminus,
            double explus,
            double eyminus,
            double eyplus, std::string source="")
      : _x(x), _y(y)
    {
      _ex = std::make_pair(exminus, explus);
      _ey[source] = std::make_pair(eyminus, eyplus);
    }


    // /// Constructor from values with symmetric errors on x and asymmetric errors on y
    // Point2D(double x, double y, double ex, const std::pair<double,double>& ey)
    //   : _x(x), _y(y), _ey(ey)
    // {
    //   _ex = std::make_pair(ex, ex);
    // }


    // /// Constructor from values with asymmetric errors on x and symmetric errors on y
    // Point2D(double x, double y, const std::pair<double,double>& ex, double ey)
    //   : _x(x), _y(y), _ex(ex)
    // {
    //   _ey = std::make_pair(ey, ey);
    // }


    /// Constructor from values with asymmetric errors on both x and y
    Point2D(double x, double y, const std::pair<double,double>& ex, const std::pair<double,double>& ey, std::string source="")
      : _x(x), _y(y)
    {
      _ex = ex;
      _ey[source] = ey;
    }


    /// Copy constructor
    Point2D(const Point2D& p)
      : _x(p._x), _y(p._y)
    {
      _ex = p._ex;
      _ey = p._ey;
      this->setParent( p.getParent() );
    }


    /// Copy assignment
    Point2D& operator = (const Point2D& p) {
      _x = p._x;
      _y = p._y;
      _ex = p._ex;
      _ey = p._ey;
      this->setParent( p.getParent() );
      return *this;
    }

    /// @}


  public:

    /// Space dimension of the point
    size_t dim() { return 2; }
    

    /// @name Value accessors
    /// @{

    /// Get x value
    double x() const { return _x; }

    /// Set x value
    void setX(double x) { _x = x; }

    /// Get y value
    double y() const { return _y; }

    /// Set y value
    void setY(double y) { _y = y; }

    /// @todo Uniform "coords" accessor across all Scatters: returning fixed-size tuple?

    /// Get x,y value pair
    std::pair<double,double> xy() const { return std::make_pair(_x, _y); }

    /// Set x and y values
    void setXY(double x, double y) { setX(x); setY(y); }

    /// Set x and y values
    void setXY(const std::pair<double,double>& xy) { setX(xy.first); setY(xy.second); }

    /// @}


    /// @name x error accessors
    /// @{


    /// Get x-error values
    const std::pair<double,double>& xErrs() const {
      return _ex;
    }

    /// Get negative x-error value
    double xErrMinus() const {
      return _ex.first;
    }

    /// Get positive x-error value
    double xErrPlus() const {
      return _ex.second;

    }

    /// Get average x-error value
    double xErrAvg() const {
      return (_ex.first + _ex.second)/2.0;
    }

    /// Set negative x error
    void setXErrMinus(double exminus) {
      _ex.first = exminus;
    }

    /// Set positive x error
    void setXErrPlus(double explus) {
      _ex.second = explus;
    }

    /// Set symmetric x error
    void setXErr(double ex) {
      setXErrMinus(ex);
      setXErrPlus(ex);
    }

    /// Set symmetric x error (alias)
    void setXErrs(double ex) {
      setXErr(ex);
    }

    /// Set asymmetric x error
    void setXErrs(double exminus, double explus) {
      setXErrMinus(exminus);
      setXErrPlus(explus);
    }

    /// Set asymmetric x error
    void setXErrs(const std::pair<double,double>& ex) {
      _ex = ex;
    }

    /// Get value minus negative x-error
    /// @todo Remove (or extend) when multiple errors are supported
    /// No: doesn't need to change since (for now) we only store multiple
    /// errors for the highest dimentsion
    double xMin() const {
      return _x - _ex.first;
    }

    /// Get value plus positive x-error
    /// @todo Remove (or extend) when multiple errors are supported
    /// No: doesn't need to change since (for now) we only store multiple
    /// errors for the highest dimentsion
    double xMax() const {
      return _x + _ex.second;
    }

    /// @}


    /// @name y error accessors
    /// @{

    /// Get y-error values
    const std::pair<double,double>& yErrs(std::string source="") const {
      if (source!="") getVariationsFromParent();
      if (!_ey.count(source)) throw RangeError("yErrs has no such key: "+source);
      return _ey.at(source);
    }

    /// Get negative y-error value
    double yErrMinus(std::string source="") const {
      if (source!="") getVariationsFromParent();
      if (!_ey.count(source)) throw RangeError("yErrs has no such key: "+source);
      return _ey.at(source).first;
    }

    /// Get positive y-error value
    double yErrPlus(std::string source="") const {
      if (source!="") getVariationsFromParent();
      if (!_ey.count(source)) throw RangeError("yErrs has no such key: "+source);
      return _ey.at(source).second;
    }

    /// Get average y-error value
    double yErrAvg(std::string source="") const {
      if (source!="") getVariationsFromParent();
      if (!_ey.count(source)) throw RangeError("yErrs has no such key: "+source);
      double res=(fabs(_ey.at(source).first) + fabs(_ey.at(source).second))/2.;
      return res;
    }

    /// Set negative y error
    void setYErrMinus(double eyminus, std::string source="") {
      if (!_ey.count(source)) {
        _ey[source] = std::make_pair(0.,0.);
      }
      _ey.at(source).first = eyminus;
    }

    /// Set positive y error
    void setYErrPlus(double eyplus, std::string source="") {
      if (!_ey.count(source)) _ey[source] = std::make_pair(0.,0.);
      _ey.at(source).second = eyplus;
    }

    /// Set symmetric y error
    void setYErr(double ey, std::string source="") {
      setYErrMinus(ey, source );
      setYErrPlus(ey, source );
    }

    /// Set symmetric y error (alias)
    void setYErrs(double ey, std::string source="") {
      setYErr(ey, source);
    }

    /// Set asymmetric y error
    void setYErrs(double eyminus, double eyplus, std::string source="") {
      setYErrMinus(eyminus, source);
      setYErrPlus(eyplus, source );
    }

    /// Set asymmetric y error
    void setYErrs(const std::pair<double,double>& ey, std::string source="") {
      _ey[source] = ey;
    }

    /// Get value minus negative y-error
    double yMin(std::string source="") const {
      if (source!="") getVariationsFromParent();
      if (!_ey.count(source)) throw RangeError("yErrs has no such key: "+source);
      return _y - _ey.at(source).first;
    }

    /// Get value plus positive y-error
    double yMax(std::string source="") const {
      if (source!="") getVariationsFromParent();
      if (!_ey.count(source)) throw RangeError("yErrs has no such key: "+source);
      return _y + _ey.at(source).second;
    }

    /// @}




    /// @name Combined x/y value and error setters
    /// @{

    /// Set x value and symmetric error
    void setX(double x, double ex) {
      setX(x);
      setXErrs(ex);
    }

    /// Set x value and asymmetric error
    void setX(double x, double exminus, double explus) {
      setX(x);
      setXErrs(exminus, explus);
    }

    /// Set x value and asymmetric error
    void setX(double x, std::pair<double,double>& ex) {
      setX(x);
      setXErrs(ex);
    }


    /// Set y value and symmetric error
    void setY(double y, double ey, std::string source="") {
      setY(y);
      setYErrs(ey, source);
    }

    /// Set y value and asymmetric error
    void setY(double y, double eyminus, double eyplus, std::string source="") {
      setY(y);
      setYErrs(eyminus, eyplus, source);
    }

    /// Set y value and asymmetric error
    void setY(double y, std::pair<double,double>& ey, std::string source="") {
      setY(y);
      setYErrs(ey, source);
    }

    /// @}


    // @name Manipulations
    /// @{

    /// Scaling of x axis
    void scaleX(double scalex) {
      setX(x()*scalex);
      setXErrs(xErrMinus()*scalex, xErrPlus()*scalex);
    }

    /// Scaling of y axis
    void scaleY(double scaley) {
      setY(y()*scaley);
      for (const auto   &source : _ey){
        setYErrs(yErrMinus()*scaley, yErrPlus()*scaley, source.first);
      }
    }

    /// Scaling of both axes
    void scaleXY(double scalex, double scaley) {
      scaleX(scalex);
      scaleY(scaley);
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


    /// @name Integer axis accessor equivalents
    /// @{

    /// Get the point value for direction @a i
    double val(size_t i) const {
      switch (i) {
      case 1: return x();
      case 2: return y();
      default: throw RangeError("Invalid axis int, must be in range 1..dim");
      }
    }
    /// Set the point value for direction @a i
    void setVal(size_t i, double val) {
      switch (i) {
      case 1: setX(val); break;
      case 2: setY(val); break;
      default: throw RangeError("Invalid axis int, must be in range 1..dim");
      }
    }

    /// Get error map for direction @a i
    const std::map< std::string, std::pair<double,double>> & errMap() const;

    /// Parse the variations from the parent AO if it exists
    void getVariationsFromParent() const;

    /// Remove the parsed variations, but keep the total
    void rmVariations() { 
      std::map< std::string, std::pair<double,double> > tmp;
      const auto& it = _ey.find("");
      if (it != _ey.end())  tmp[""] = it->second;
      _ey.clear();  _ey = tmp;
    }
    
    /// set the "" error source to the sum in quad of the existing variations 
    void updateTotalUncertainty() {
        float totalUp = 0.;
        float totalDn = 0.;
        for (const auto& variation : getParent()->variations()) {
          if (variation=="") continue;
          float thisUp = yErrPlus(variation);
          float thisDn = yErrMinus(variation);
          totalUp += thisUp*thisUp;
          totalDn += thisDn*thisDn;
        }
    
        totalUp = sqrt(totalUp);
        totalDn = sqrt(totalDn);
        setErrPlus(2, totalUp);
        setErrMinus(2, totalDn);
    }

    /// Get error values for direction @a i
    const std::pair<double,double>& errs(size_t i, std::string source="") const {
      switch (i) {
      case 1: return xErrs();
      case 2: return yErrs(source);
      default: throw RangeError("Invalid axis int, must be in range 1..dim");
      }
    }
    /// Get negative error value for direction @a i
    double errMinus(size_t i, std::string source="") const {
      switch (i) {
      case 1: return xErrMinus();
      case 2: return yErrMinus(source);
      default: throw RangeError("Invalid axis int, must be in range 1..dim");
      }
    }
    /// Get positive error value for direction @a i
    double errPlus(size_t i, std::string source="") const {
      switch (i) {
      case 1: return xErrPlus();
      case 2: return yErrPlus(source);
      default: throw RangeError("Invalid axis int, must be in range 1..dim");
      }
    }
    /// Get average error value for direction @a i
    double errAvg(size_t i, std::string source="") const {
      switch (i) {
      case 1: return xErrAvg();
      case 2: return yErrAvg(source);
      default: throw RangeError("Invalid axis int, must be in range 1..dim");
      }
    }

    /// Set negative error for direction @a i
    void setErrMinus(size_t i, double eminus, std::string source="") {
      switch (i) {
      case 1: setXErrMinus(eminus); break;
      case 2: setYErrMinus(eminus, source); break;
      default: throw RangeError("Invalid axis int, must be in range 1..dim");
      }
    }
    /// Set positive error for direction @a i
    void setErrPlus(size_t i, double eplus, std::string source="") {
      switch (i) {
      case 1: setXErrPlus(eplus); break;
      case 2: setYErrPlus(eplus, source); break;
      default: throw RangeError("Invalid axis int, must be in range 1..dim");
      }
    }

    /// Set symmetric error for direction @a i
    void setErr(size_t i, double e, std::string source="") {
      switch (i) {
      case 1: setXErrs(e); break;
      case 2: setYErrs(e, source); break;
      default: throw RangeError("Invalid axis int, must be in range 1..dim");
      }
    }
    /// Set asymmetric error for direction @a i
    void setErrs(size_t i, double eminus, double eplus, std::string source="") {
      switch (i) {
      case 1: setXErrs(eminus, eplus); break;
      case 2: setYErrs(eminus, eplus, source); break;
      default: throw RangeError("Invalid axis int, must be in range 1..dim");
      }
    }
    /// Set asymmetric error for direction @a i
    void setErrs(size_t i, std::pair<double,double>& e, std::string source="") {
      switch (i) {
      case 1: setXErrs(e); break;
      case 2: setYErrs(e, source); break;
      default: throw RangeError("Invalid axis int, must be in range 1..dim");
      }
    }

    /// Set value and symmetric error for direction @a i
    void set(size_t i, double val, double e,  std::string source="") {
      switch (i) {
      case 1: setX(val, e); break;
      case 2: setY(val, e, source); break;
      default: throw RangeError("Invalid axis int, must be in range 1..dim");
      }
    }
    /// Set value and asymmetric error for direction @a i
    void set(size_t i, double val, double eminus, double eplus,  std::string source="") {
      switch (i) {
      case 1: setX(val, eminus, eplus); break;
      case 2: setY(val, eminus, eplus, source); break;
      default: throw RangeError("Invalid axis int, must be in range 1..dim");
      }
    }
    /// Set value and asymmetric error for direction @a i
    void set(size_t i, double val, std::pair<double,double>& e,  std::string source="") {
      switch (i) {
      case 1: setX(val, e); break;
      case 2: setY(val, e, source); break;
      default: throw RangeError("Invalid axis int, must be in range 1..dim");
      }
    }


    /// @}


  protected:

    /// @name Value and error variables
    /// @{


    double _x;
    double _y;
    std::pair<double,double> _ex;
    // a map of the errors for each source. Nominal stored under ""
    // to ensure backward compatibility
    std::map< std::string, std::pair<double,double> > _ey;

    /// @}

  };



  /// @name Comparison operators
  /// @{

  /// Equality test of x & y characteristics only
  /// @todo Base on a named fuzzyEquals(a,b,tol=1e-3) unbound function
  inline bool operator==(const YODA::Point2D& a, const YODA::Point2D& b) {
    if (!YODA::fuzzyEquals(a.x(), b.x()) ||
        !YODA::fuzzyEquals(a.xErrMinus(), b.xErrMinus()) ||
        !YODA::fuzzyEquals(a.xErrPlus(),  b.xErrPlus()) ) return false;
    if (!YODA::fuzzyEquals(a.y(), b.y()) ||
        !YODA::fuzzyEquals(a.yErrMinus(), b.yErrMinus()) ||
        !YODA::fuzzyEquals(a.yErrPlus(),  b.yErrPlus()) ) return false;
    return true;
  }

  /// Equality test of x characteristics only
  inline bool operator != (const YODA::Point2D& a, const YODA::Point2D& b) {
    return !(a == b);
  }

  /// Less-than operator used to sort bins by x-ordering
  inline bool operator < (const YODA::Point2D& a, const YODA::Point2D& b) {
    if (!YODA::fuzzyEquals(a.x(), b.x())) {
      return a.x() < b.x();
    }
    if (!YODA::fuzzyEquals(a.xErrMinus(), b.xErrMinus())) {
      return a.xErrMinus() < b.xErrMinus();
    }
    if (!YODA::fuzzyEquals(a.xErrPlus(), b.xErrPlus())) {
      return a.xErrPlus() < b.xErrPlus();
    }
    return false;
  }

  /// Less-than-or-equals operator used to sort bins by x-ordering
  inline bool operator <= (const YODA::Point2D& a, const YODA::Point2D& b) {
    if (a == b) return true;
    return a < b;
  }

  /// Greater-than operator used to sort bins by x-ordering
  inline bool operator > (const YODA::Point2D& a, const YODA::Point2D& b) {
    return !(a <= b);
  }

  /// Greater-than-or-equals operator used to sort bins by x-ordering
  inline bool operator >= (const YODA::Point2D& a, const YODA::Point2D& b) {
    return !(a < b);
  }

  /// @}


}

#endif
