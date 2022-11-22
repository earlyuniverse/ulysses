// -*- C++ -*-
//
// This file is part of YODA -- Yet more Objects for Data Analysis
// Copyright (C) 2008-2021 The YODA collaboration (see AUTHORS for details)
//
#ifndef YODA_POINT1D_H
#define YODA_POINT1D_H

#include "YODA/Point.h"
#include "YODA/Exceptions.h"
#include "YODA/Utils/MathUtils.h"
#include <utility>

namespace YODA {


  /// A 1D data point to be contained in a Scatter1D
  class Point1D : public Point {
  public:

    /// @name Constructors
    /// @{

    // Default constructor
    Point1D() {  }


    /// Constructor from values with optional symmetric errors
    Point1D(double x, double ex=0.0, std::string source="")
      : _x(x)
    {
      _ex[source] = std::make_pair(ex, ex);
    }


    /// Constructor from values with explicit asymmetric errors
    Point1D(double x, double exminus, double explus, std::string source="")
      : _x(x)
    {
      _ex[source] = std::make_pair(exminus, explus);
    }


    /// Constructor from values with asymmetric errors
    Point1D(double x, const std::pair<double,double>& ex,  std::string source="")
      : _x(x)
    {
      _ex[source] = ex;
    }


    /// Copy constructor
    Point1D(const Point1D& p)
      : _x(p._x), _ex(p._ex)
    {
      this->setParent( p.getParent() );
    }


    /// Copy assignment
    Point1D& operator = (const Point1D& p) {
      _x = p._x;
      _ex = p._ex;
      this->setParent( p.getParent() );
      return *this;
    }

    /// @}


  public:

    /// Space dimension of the point
    size_t dim() { return 1; }

    
    /// @name Value accessors
    /// @{

    /// Get x value
    double x() const { return _x; }

    /// Set x value
    void setX(double x) { _x = x; }

    /// @todo Uniform "coords" accessor across all Scatters: returning fixed-size tuple?

    /// @}


    /// @name x error accessors
    /// @{

    /// Get x-error values
    const std::pair<double,double>& xErrs(  std::string source="") const {
      if (source!="") getVariationsFromParent();
      if (!_ex.count(source)) throw RangeError("xErrs has no such key: "+source);
      return _ex.at(source);
    }

    /// Get negative x-error value
    double xErrMinus( std::string source="") const {
      if (source!="") getVariationsFromParent();
      if (!_ex.count(source)) throw RangeError("xErrs has no such key: "+source);
      return _ex.at(source).first;
    }

    /// Get positive x-error value
    double xErrPlus( std::string source="") const {
      if (source!="") getVariationsFromParent();
      if (!_ex.count(source)) throw RangeError("xErrs has no such key: "+source);
      return _ex.at(source).second;
    }

    /// Get average x-error value
    double xErrAvg( std::string source="") const {
      if (source!="") getVariationsFromParent();
      if (!_ex.count(source)) throw RangeError("xErrs has no such key: "+source);
      return (_ex.at(source).first + _ex.at(source).second)/2.0;
    }

    /// Set negative x error
    void setXErrMinus(double exminus,  std::string source="") {
      if (!_ex.count(source)) _ex[source] = std::make_pair(0.,0.);
      _ex.at(source).first = exminus;
    }

    /// Set positive x error
    void setXErrPlus(double explus,  std::string source="") {
      if (!_ex.count(source)) _ex[source] = std::make_pair(0.,0.);
      _ex.at(source).second = explus;
    }

    /// Set symmetric x error
    void setXErr(double ex,  std::string source="") {
      setXErrMinus(ex, source);
      setXErrPlus(ex, source);
    }

    /// Set symmetric x error (alias)
    void setXErrs(double ex,  std::string source="") {
      setXErr(ex, source);
    }

    /// Set asymmetric x error
    void setXErrs(double exminus, double explus,  std::string source="") {
      setXErrMinus(exminus, source);
      setXErrPlus(explus, source);
    }

    /// Set asymmetric x error
    void setXErrs(const std::pair<double,double>& ex,  std::string source="") {
      _ex[source] = ex;
    }

    /// Get value minus negative x-error
    double xMin(std::string source="") const {
      if (source!="") getVariationsFromParent();
      if (!_ex.count(source)) throw RangeError("xErrs has no such key: "+source);
      return _x - _ex.at(source).first;
    }

    /// Get value plus positive x-error
    double xMax(std::string source="") const {
      if (source!="") getVariationsFromParent();
      if (!_ex.count(source)) throw RangeError("xErrs has no such key: "+source);
      return _x + _ex.at(source).second;
    }

    /// @}


    /// @name Combined x value and error setters
    /// @{

    /// Set x value and symmetric error
    void setX(double x, double ex, std::string source="") {
      setX(x);
      setXErr(ex, source);
    }

    /// Set x value and asymmetric error
    void setX(double x, double exminus, double explus, std::string source="") {
      setX(x);
      setXErrs(exminus, explus, source);
    }

    /// Set x value and asymmetric error
    void setX(double x, std::pair<double,double>& ex, std::string source="") {
      setX(x);
      setXErrs(ex, source);
    }

    /// @}


    // @name Manipulations
    /// @{

    /// Scaling of x axis
    void scaleX(double scalex) {
      setX(x()*scalex);
      for (const auto   &source : _ex){
        setXErrs(xErrMinus()*scalex, xErrPlus()*scalex, source.first);
      }
    }

    /// Scaling along direction @a i
    void scale(size_t i, double scale) {
      switch (i) {
      case 1: scaleX(scale); break;
      default: throw RangeError("Invalid axis int, must be in range 1..dim");
      }
    }

    /// @}


    /// @name Integer axis accessor equivalents
    /// @{

    /// Get the point value for direction @a i
    double val(size_t i) const {
      if (i == 0 || i > 1) throw RangeError("Invalid axis int, must be in range 1..dim");
      return x();
    }
    /// Set the point value for direction @a i
    void setVal(size_t i, double val) {
      if (i != 1) throw RangeError("Invalid axis int, must be in range 1..dim");
      setX(val);
    }

    /// Get error map for direction @a i
    const std::map< std::string, std::pair<double,double>> & errMap() const {
      getVariationsFromParent();
      return _ex;
    }
    
    /// Remove the parsed variations, but keep the total
    void rmVariations() { 
      std::map< std::string, std::pair<double,double> > tmp;
      const auto& it = _ex.find("");
      if (it != _ex.end())  tmp[""] = it->second;
      _ex.clear();  _ex = tmp;
    }

    // set the "" error source to the sum in quad of the existing variations 
    void updateTotalUncertainty() {
        float totalUp = 0.;
        float totalDn = 0.;
        for (const auto& variation : getParent()->variations()) {
          if (variation=="") continue;
          float thisUp =  xErrPlus(variation);
          float thisDn =  xErrMinus(variation);
          totalUp += thisUp*thisUp;
          totalDn += thisDn*thisDn;
        }
    
        totalUp = sqrt(totalUp);
        totalDn = sqrt(totalDn);
        setErrPlus(1, totalUp);
        setErrMinus(1, totalDn);
    }

    // Parse the variations from the parent AO if it exists
    void getVariationsFromParent() const;


    /// Get error values for direction @a i
    const std::pair<double,double>& errs(size_t i, std::string source="") const {
      if (i != 1) throw RangeError("Invalid axis int, must be in range 1..dim");
      return xErrs(source);
    }
    /// Get negative error value for direction @a i
    double errMinus(size_t i, std::string source="") const {
      if (i != 1) throw RangeError("Invalid axis int, must be in range 1..dim");
      return xErrMinus(source);
    }
    /// Get positive error value for direction @a i
    double errPlus(size_t i,  std::string source="") const {
      if (i != 1) throw RangeError("Invalid axis int, must be in range 1..dim");
      return xErrPlus(source);
    }
    /// Get average error value for direction @a i
    double errAvg(size_t i, std::string source="") const {
      if (i != 1) throw RangeError("Invalid axis int, must be in range 1..dim");
      return xErrAvg(source);
    }

    /// Set negative error for direction @a i
    void setErrMinus(size_t i, double eminus,  std::string source="") {
      if (i != 1) throw RangeError("Invalid axis int, must be in range 1..dim");
      setXErrMinus(eminus, source);
    }
    /// Set positive error for direction @a i
    void setErrPlus(size_t i, double eplus,  std::string source="") {
      if (i != 1) throw RangeError("Invalid axis int, must be in range 1..dim");
      setXErrPlus(eplus, source);
    }

    /// Set symmetric error for direction @a i
    void setErr(size_t i, double e, std::string source="") {
      if (i != 1) throw RangeError("Invalid axis int, must be in range 1..dim");
      setXErr(e, source);
    }
    /// Set asymmetric error for direction @a i
    void setErrs(size_t i, double eminus, double eplus,  std::string source="") {
      if (i != 1) throw RangeError("Invalid axis int, must be in range 1..dim");
      setXErrs(eminus, eplus, source);
    }
    /// Set asymmetric error for direction @a i
    void setErrs(size_t i, std::pair<double,double>& e, std::string source="") {
      if (i != 1) throw RangeError("Invalid axis int, must be in range 1..dim");
      setXErrs(e, source);
    }

    /// Set value and symmetric error for direction @a i
    void set(size_t i, double val, double e,  std::string source="") {
      if (i != 1) throw RangeError("Invalid axis int, must be in range 1..dim");
      setX(val, e, source);
    }
    /// Set value and asymmetric error for direction @a i
    void set(size_t i, double val, double eminus, double eplus,  std::string source="") {
      if (i != 1) throw RangeError("Invalid axis int, must be in range 1..dim");
      setX(val, eminus, eplus, source);
    }
    /// Set value and asymmetric error for direction @a i
    void set(size_t i, double val, std::pair<double,double>& e,  std::string source="") {
      if (i != 1) throw RangeError("Invalid axis int, must be in range 1..dim");
      setX(val, e, source);
    }

    /// @}


  protected:

    /// @name Value and error variables
    /// @{

    double _x;
    // a map of the errors for each source. Nominal stored under ""
    // to ensure backward compatibility
    std::map< std::string, std::pair<double,double> > _ex;

    /// @}

  };



  /// @name Comparison operators
  /// @{

  /// Equality test of x characteristics only
  /// @todo Base on a named fuzzyEquals(a,b,tol=1e-3) unbound function
  inline bool operator==(const YODA::Point1D& a, const YODA::Point1D& b) {
    if (!YODA::fuzzyEquals(a.x(), b.x()) ||
        !YODA::fuzzyEquals(a.xErrMinus(), b.xErrMinus()) ||
        !YODA::fuzzyEquals(a.xErrPlus(),  b.xErrPlus()) ) return false;
    return true;
  }

  /// Equality test of x characteristics only
  inline bool operator != (const YODA::Point1D& a, const YODA::Point1D& b) {
    return !(a == b);
  }

  /// Less-than operator used to sort bins by x-ordering
  inline bool operator < (const YODA::Point1D& a, const YODA::Point1D& b) {
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
  inline bool operator <= (const YODA::Point1D& a, const YODA::Point1D& b) {
    if (a == b) return true;
    return a < b;
  }

  /// Greater-than operator used to sort bins by x-ordering
  inline bool operator > (const YODA::Point1D& a, const YODA::Point1D& b) {
    return !(a <= b);
  }

  /// Greater-than-or-equals operator used to sort bins by x-ordering
  inline bool operator >= (const YODA::Point1D& a, const YODA::Point1D& b) {
    return !(a < b);
  }

  /// @}


}

#endif
