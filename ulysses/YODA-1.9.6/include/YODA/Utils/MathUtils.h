// -*- C++ -*-
//
// This file is part of YODA -- Yet more Objects for Data Analysis
// Copyright (C) 2008-2017 The YODA collaboration (see AUTHORS for details)
//
#ifndef YODA_MathUtils_H
#define YODA_MathUtils_H

/// @todo Add SFINAE math type stuff (see Rivet) and add inrange() and inrange_closed_closed() etc. aliases cf. MCUtils

#include "YODA/Config/BuildConfig.h"

#include <stdexcept>
#include <string>
#include <ostream>
#include <sstream>
#include <iostream>
#include <cmath>
#include <map>
#include <vector>
#include <utility>
#include <algorithm>
#include <functional>
#include <numeric>
#include <cassert>
#include <limits>
#include <climits>
#include <cfloat>

namespace YODA {


  /// Pre-defined numeric type limits
  /// @deprecated Prefer the standard DBL/INT_MAX
  const static double MAXDOUBLE = DBL_MAX; // was std::numeric_limits<double>::max(); -- warns in GCC5
  const static double MAXINT = INT_MAX; // was std::numeric_limits<int>::max(); -- warns in GCC5

  /// A pre-defined value of \f$ \pi \f$.
  static const double PI = M_PI;

  /// A pre-defined value of \f$ 2\pi \f$.
  static const double TWOPI = 2*M_PI;

  /// A pre-defined value of \f$ \pi/2 \f$.
  static const double HALFPI = M_PI_2;

  /// Enum for signs of numbers.
  enum Sign { MINUS = -1, ZERO = 0, PLUS = 1 };


  /// @name Comparison functions for safe floating point equality tests
  //@{

  /// Compare a floating point number to zero with a degree
  /// of fuzziness expressed by the absolute @a tolerance parameter.
  inline bool isZero(double val, double tolerance=1E-8) {
    return (fabs(val) < tolerance);
  }

  /// Compare an integral-type number to zero.
  ///
  /// Since there is no risk of floating point error, this function just exists
  /// in case @c isZero is accidentally used on an integer type, to avoid
  /// implicit type conversion. The @a tolerance parameter is ignored.
  inline bool isZero(long val, double=1E-8) {
    return val == 0;
  }


  /// @brief Compare two floating point numbers for equality with a degree of fuzziness
  ///
  /// The @a tolerance parameter is fractional, based on absolute values of the args.
  inline bool fuzzyEquals(double a, double b, double tolerance=1E-5) {
   const double absavg = (fabs(a) + fabs(b))/2.0;
   const double absdiff = fabs(a - b);
   const bool rtn = (isZero(a) && isZero(b)) || absdiff < tolerance*absavg;
   // cout << a << " == " << b << "? " << rtn << endl;
   return rtn;
  }


  /// @brief Compare two integral-type numbers for equality with a degree of fuzziness.
  ///
  /// Since there is no risk of floating point error with integral types,
  /// this function just exists in case @c fuzzyEquals is accidentally
  /// used on an integer type, to avoid implicit type conversion. The @a
  /// tolerance parameter is ignored, even if it would have an
  /// absolute magnitude greater than 1.
  inline bool fuzzyEquals(long a, long b, double=1E-5) {
    return a == b;
  }


  /// @brief Compare two floating point numbers for >= with a degree of fuzziness
  ///
  /// The @a tolerance parameter on the equality test is as for @c fuzzyEquals.
  inline bool fuzzyGtrEquals(double a, double b, double tolerance=1E-5) {
    return a > b || fuzzyEquals(a, b, tolerance);
  }

  /// @brief Compare two integral-type numbers for >= with a degree of fuzziness.
  ///
  /// Since there is no risk of floating point error with integral types,
  /// this function just exists in case @c fuzzyGtrEquals is accidentally
  /// used on an integer type, to avoid implicit type conversion. The @a
  /// tolerance parameter is ignored, even if it would have an
  /// absolute magnitude greater than 1.
  inline bool fuzzyGtrEquals(long a, long b, double=1E-5) {
    return a >= b;
  }


  /// @brief Compare two floating point numbers for <= with a degree of fuzziness
  ///
  /// The @a tolerance parameter on the equality test is as for @c fuzzyEquals.
  inline bool fuzzyLessEquals(double a, double b, double tolerance=1E-5) {
    return a < b || fuzzyEquals(a, b, tolerance);
  }

  /// @brief Compare two integral-type numbers for <= with a degree of fuzziness.
  ///
  /// Since there is no risk of floating point error with integral types,
  /// this function just exists in case @c fuzzyLessEquals is accidentally
  /// used on an integer type, to avoid implicit type conversion. The @a
  /// tolerance parameter is ignored, even if it would have an
  /// absolute magnitude greater than 1.
  inline bool fuzzyLessEquals(long a, long b, double=1E-5) {
    return a <= b;
  }

  /// Returns a number floored at the nth decimal place.
  inline double approx(double a, int n = 5) {
    double roundTo = pow(10.0,n);
    a *= roundTo;
    a = floor(a);
    return a/roundTo;
  }
  //@}


  /// @name Ranges and intervals
  //@{

  /// Represents whether an interval is open (non-inclusive) or closed (inclusive).
  ///
  /// For example, the interval \f$ [0, \pi) \f$ is closed (an inclusive
  /// boundary) at 0, and open (a non-inclusive boundary) at \f$ \pi \f$.
  enum RangeBoundary { OPEN=0, SOFT=0, CLOSED=1, HARD=1 };


  /// @brief Determine if @a value is in the range @a low to @a high, for floating point numbers
  ///
  /// Interval boundary types are defined by @a lowbound and @a highbound.
  /// @todo Optimise to one-line at compile time?
  template<typename NUM>
  inline bool inRange(NUM value, NUM low, NUM high,
                      RangeBoundary lowbound=CLOSED, RangeBoundary highbound=OPEN) {
    if (lowbound == OPEN && highbound == OPEN) {
      return (value > low && value < high);
    } else if (lowbound == OPEN && highbound == CLOSED) {
      return (value > low && value <= high);
    } else if (lowbound == CLOSED && highbound == OPEN) {
      return (value >= low && value < high);
    } else { // if (lowbound == CLOSED && highbound == CLOSED) {
      return (value >= low && value <= high);
    }
  }

  /// Alternative version of inRange for doubles, which accepts a pair for the range arguments.
  template<typename NUM>
  inline bool inRange(NUM value, std::pair<NUM, NUM> lowhigh,
                      RangeBoundary lowbound=CLOSED, RangeBoundary highbound=OPEN) {
    return inRange(value, lowhigh.first, lowhigh.second, lowbound, highbound);
  }


  /// @brief Determine if @a value is in the range @a low to @a high, for integer types
  ///
  /// Interval boundary types are defined by @a lowbound and @a highbound.
  /// @todo Optimise to one-line at compile time?
  inline bool inRange(int value, int low, int high,
                      RangeBoundary lowbound=CLOSED, RangeBoundary highbound=CLOSED) {
    if (lowbound == OPEN && highbound == OPEN) {
      return (value > low && value < high);
    } else if (lowbound == OPEN && highbound == CLOSED) {
      return (value > low && value <= high);
    } else if (lowbound == CLOSED && highbound == OPEN) {
      return (value >= low && value < high);
    } else { // if (lowbound == CLOSED && highbound == CLOSED) {
      return (value >= low && value <= high);
    }
  }

  /// Alternative version of @c inRange for ints, which accepts a pair for the range arguments.
  inline bool inRange(int value, std::pair<int, int> lowhigh,
                      RangeBoundary lowbound=CLOSED, RangeBoundary highbound=OPEN) {
    return inRange(value, lowhigh.first, lowhigh.second, lowbound, highbound);
  }

  //@}


  /// @name Miscellaneous numerical helpers
  //@{

  /// Named number-type squaring operation.
  template <typename NUM>
  inline NUM sqr(NUM a) {
    return a*a;
  }

  /// Named number-type addition in quadrature operation.
  template <typename Num>
  inline Num add_quad(Num a, Num b) {
    return sqrt(a*a + b*b);
  }

  /// Named number-type addition in quadrature operation.
  template <typename Num>
  inline Num add_quad(Num a, Num b, Num c) {
    return sqrt(a*a + b*b + c*c);
  }

  /// Find the sign of a number
  inline int sign(double val) {
    if (isZero(val)) return ZERO;
    const int valsign = (val > 0) ? PLUS : MINUS;
    return valsign;
  }

  /// Find the sign of a number
  inline int sign(int val) {
    if (val == 0) return ZERO;
    return (val > 0) ? PLUS : MINUS;
  }

  /// Find the sign of a number
  inline int sign(long val) {
    if (val == 0) return ZERO;
    return (val > 0) ? PLUS : MINUS;
  }

  //@}


  /// @name Binning helper functions
  //@{

  /// @brief Make a list of @a nbins + 1 values uniformly spaced between @a xmin and @a xmax inclusive.
  ///
  /// @note The arg ordering and the meaning of the nbins variable is "histogram-like",
  /// as opposed to the Numpy/Matlab version.
  inline std::vector<double> linspace(size_t nbins, double xmin, double xmax, bool include_end=true) {
    assert(xmax >= xmin);
    assert(nbins > 0);
    std::vector<double> rtn;
    const double interval = (xmax-xmin)/static_cast<double>(nbins);
    for (size_t i = 0; i < nbins; ++i) {
      rtn.push_back(xmin + i*interval);
    }
    assert(rtn.size() == nbins);
    if (include_end) rtn.push_back(xmax); // exact xmax, not result of n * interval
    return rtn;
  }


  /// @brief Make a list of @a nbins + 1 values uniformly spaced in log(x) between @a xmin and @a xmax inclusive.
  ///
  /// @note The arg ordering and the meaning of the nbins variable is "histogram-like",
  /// as opposed to the Numpy/Matlab version, and the xmin and xmax arguments are expressed
  /// in "normal" space, rather than as the logarithms of the xmin/xmax values as in Numpy/Matlab.
  inline std::vector<double> logspace(size_t nbins, double xmin, double xmax, bool include_end=true) {
    assert(xmax >= xmin);
    assert(xmin > 0);
    assert(nbins > 0);
    const double logxmin = std::log(xmin);
    const double logxmax = std::log(xmax);
    const std::vector<double> logvals = linspace(nbins, logxmin, logxmax);
    assert(logvals.size() == nbins+1);
    std::vector<double> rtn; rtn.reserve(logvals.size());
    rtn.push_back(xmin);
    for (size_t i = 1; i < logvals.size()-1; ++i) {
      rtn.push_back(std::exp(logvals[i]));
    }
    assert(rtn.size() == nbins);
    if (include_end) rtn.push_back(xmax);
    return rtn;
  }


  /// @todo fspace() for uniform sampling from f(x); requires ability to invert fn... how, in general?
  //inline std::vector<double> fspace(size_t nbins, double xmin, double xmax, std::function<double(double)>& fn) {


  /// @brief Make a list of @a nbins + 1 values spaced with *density* ~ f(x) between @a xmin and @a end inclusive.
  ///
  /// The density function @a fn will be evaluated at @a nsample uniformly
  /// distributed points between @a xmin and @a xmax, its integral approximated
  /// via the Trapezium Rule and used to normalize the distribution, and @a
  /// nbins + 1 edges then selected to (approximately) divide into bins each
  /// containing fraction 1/@a nbins of the integral.
  ///
  /// @note The function @a fn does not need to be a normalized pdf, but it should be non-negative.
  /// Any negative values returned by @a fn(x) will be truncated to zero and contribute nothing to
  /// the estimated integral and binning density.
  ///
  /// @note The arg ordering and the meaning of the nbins variable is "histogram-like",
  /// as opposed to the Numpy/Matlab version, and the xmin and xmax arguments are expressed
  /// in "normal" space, rather than in the function-mapped space as with Numpy/Matlab.
  ///
  /// @note The naming of this function differs from the other, Matlab-inspired ones: the bin spacing is
  /// uniform in the CDF of the density function given, rather than in the function itself. For
  /// most use-cases this is more intuitive.
  inline std::vector<double> pdfspace(size_t nbins, double xmin, double xmax, std::function<double(double)>& fn, size_t nsample=10000) {
    const double dx = (xmax-xmin)/(double)nsample;
    const std::vector<double> xs = linspace(nsample, xmin, xmax);
    std::vector<double> ys(0, nsample);
    auto posfn = [&](double x){return std::max(fn(x), 0.0);};
    std::transform(xs.begin(), xs.end(), ys.begin(), posfn);
    std::vector<double> areas; areas.reserve(nsample);
    double areasum = 0;
    for (size_t i = 0; i < ys.size()-1; ++i) {
      const double area = (ys[i] + ys[i+1])*dx/2.0;
      areas[i] = area;
      areasum += area;
    }
    const double df = areasum/(double)nbins;
    std::vector<double> xedges{xmin}; xedges.reserve(nbins+1);
    double fsum = 0;
    for (size_t i = 0; i < nsample-1; ++i) {
      fsum += areas[i];
      if (fsum > df) {
        fsum = 0;
        xedges.push_back(xs[i+1]);
      }
    }
    xedges.push_back(xmax);
    assert(xedges.size() == nbins+1);
    return xedges;
  }


  /// @brief Return the bin index of the given value, @a val, given a vector of bin edges
  ///
  /// NB. The @a binedges vector must be sorted
  template <typename NUM>
  inline int index_between(const NUM& val, const std::vector<NUM>& binedges) {
    if (!inRange(val, binedges.front(), binedges.back())) return -1; //< Out of histo range
    int index = -1;
    for (size_t i = 1; i < binedges.size(); ++i) {
      if (val < binedges[i]) {
        index = i-1;
        break;
      }
    }
    assert(inRange(index, -1, binedges.size()-1));
    return index;
  }

  //@}


  /// @name Statistics functions
  //@{

  /// Calculate the mean of a sample
  inline double mean(const std::vector<int>& sample) {
    double mean = 0.0;
    for (size_t i=0; i<sample.size(); ++i) {
      mean += sample[i];
    }
    return mean/sample.size();
  }


  /// Calculate the covariance (variance) between two samples
  inline double covariance(const std::vector<int>& sample1, const std::vector<int>& sample2) {
    const double mean1 = mean(sample1);
    const double mean2 = mean(sample2);
    const size_t N = sample1.size();
    double cov = 0.0;
    for (size_t i = 0; i < N; i++) {
      const double cov_i = (sample1[i] - mean1)*(sample2[i] - mean2);
      cov += cov_i;
    }
    if (N > 1) return cov/(N-1);
    else return 0.0;
  }


  /// Calculate the correlation strength between two samples
  inline double correlation(const std::vector<int>& sample1, const std::vector<int>& sample2) {
    const double cov = covariance(sample1, sample2);
    const double var1 = covariance(sample1, sample1);
    const double var2 = covariance(sample2, sample2);
    const double correlation = cov/sqrt(var1*var2);
    const double corr_strength = correlation*sqrt(var2/var1);
    return corr_strength;
  }

  //@}


}

#endif
