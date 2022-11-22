// -*- C++ -*-
//
// This file is part of YODA -- Yet more Objects for Data Analysis
// Copyright (C) 2008-2021 The YODA collaboration (see AUTHORS for details)
//
#include "YODA/Profile1D.h"
#include "YODA/Histo1D.h"
#include "YODA/Scatter2D.h"

#include <cmath>
using namespace std;

namespace YODA {


  void Profile1D::fill(double x, double y, double weight, double fraction) {
    if ( std::isnan(x) ) throw RangeError("X is NaN");
    if ( std::isnan(y) ) throw RangeError("Y is NaN");

    // Fill the overall distribution
    _axis.totalDbn().fill(x, y, weight, fraction);

    // Fill the bins and overflows
    /// Unify this with Histo1D's version, when binning and inheritance are reworked
    if (inRange(x, _axis.xMin(), _axis.xMax())) {
      try {
        /// @todo Replace try block with a check that there is a bin at x
        _binAt(x).fill(x, y, weight, fraction);
      } catch (const RangeError& re) {    }
    } else if (x < _axis.xMin()) {
      _axis.underflow().fill(x, y, weight, fraction);
    } else if (x >= _axis.xMax()) {
      _axis.overflow().fill(x, y, weight, fraction);
    }

    // Lock the axis now that a fill has happened
    _axis._setLock(true);
  }


  void Profile1D::fillBin(size_t i, double y, double weight, double fraction) {
    fill(bin(i).xMid(), y, weight, fraction);
  }



  /////////////// COMMON TO ALL BINNED

  double Profile1D::numEntries(bool includeoverflows) const {
    if (includeoverflows) return totalDbn().numEntries();
    unsigned long n = 0;
    for (const Bin& b : bins()) n += b.numEntries();
    return n;
  }


  double Profile1D::effNumEntries(bool includeoverflows) const {
    if (includeoverflows) return totalDbn().effNumEntries();
    double n = 0;
    for (const Bin& b : bins()) n += b.effNumEntries();
    return n;
  }


  double Profile1D::sumW(bool includeoverflows) const {
    if (includeoverflows) return _axis.totalDbn().sumW();
    double sumw = 0;
    for (const Bin& b : bins()) sumw += b.sumW();
    return sumw;
  }


  double Profile1D::sumW2(bool includeoverflows) const {
    if (includeoverflows) return _axis.totalDbn().sumW2();
    double sumw2 = 0;
    for (const Bin& b : bins()) sumw2 += b.sumW2();
    return sumw2;
  }

  // ^^^^^^^^^^^^^^^^^^


  double Profile1D::xMean(bool includeoverflows) const {
    if (includeoverflows) return _axis.totalDbn().xMean();
    Dbn2D dbn;
    for (const ProfileBin1D& b : bins()) dbn += b.dbn();
    return dbn.xMean();
  }


  double Profile1D::xVariance(bool includeoverflows) const {
    if (includeoverflows) return _axis.totalDbn().xVariance();
    Dbn2D dbn;
    for (const ProfileBin1D& b : bins()) dbn += b.dbn();
    return dbn.xVariance();
  }


  double Profile1D::xStdErr(bool includeoverflows) const {
    if (includeoverflows) return _axis.totalDbn().xStdErr();
    Dbn2D dbn;
    for (const ProfileBin1D& b : bins()) dbn += b.dbn();
    return dbn.xStdErr();
  }


  double Profile1D::xRMS(bool includeoverflows) const {
    if (includeoverflows) return _axis.totalDbn().xRMS();
    Dbn2D dbn;
    for (const ProfileBin1D& b : bins()) dbn += b.dbn();
    return dbn.xRMS();
  }


  ////////////////////////////////////////


  /// Copy constructor with optional new path
  Profile1D::Profile1D(const Profile1D& p, const std::string& path)
    : AnalysisObject("Profile1D", (path.size() == 0) ? p.path() : path, p, p.title())
  {
    _axis = p._axis;
  }


  /// Constructor from a Scatter2D's binning, with optional new path
  Profile1D::Profile1D(const Scatter2D& s, const std::string& path)
    : AnalysisObject("Profile1D", (path.size() == 0) ? s.path() : path, s, s.title())
  {
    std::vector<ProfileBin1D> bins;
    for (const Scatter2D::Point& p : s.points()) {
      bins.push_back(ProfileBin1D(p.xMin(), p.xMax()));
    }
    _axis = Profile1DAxis(bins);
  }


  /// Constructor from a Histo1D's binning, with optional new path
  Profile1D::Profile1D(const Histo1D& h, const std::string& path)
    : AnalysisObject("Profile1D", (path.size() == 0) ? h.path() : path, h, h.title())
  {
    Bins bins;
    for (const Histo1D::Bin& b : h.bins()) {
      bins.push_back(ProfileBin1D(b.xMin(), b.xMax()));
    }
    _axis = Profile1DAxis(bins);

  }


  ////////////////////////////////////////


  /// Divide two profile histograms
  Scatter2D divide(const Profile1D& numer, const Profile1D& denom) {
    Scatter2D rtn;

    for (size_t i = 0; i < numer.numBins(); ++i) {
      const ProfileBin1D& b1 = numer.bin(i);
      const ProfileBin1D& b2 = denom.bin(i);

      /// @todo Create a compatibleBinning function? Or just compare vectors of edges().
      if (!fuzzyEquals(b1.xMin(), b2.xMin()) || !fuzzyEquals(b1.xMax(), b2.xMax()))
        throw BinningError("x binnings are not equivalent in " + numer.path() + " / " + denom.path());

      // Assemble the x value and error
      // Use the midpoint of the "bin" for the new central x value, in the absence of better information
      const double x = b1.xMid();
      const double exminus = x - b1.xMin();
      const double explus  = b1.xMax() - x;

      // Assemble the y value and error
      double y = std::numeric_limits<double>::quiet_NaN();
      double ey = std::numeric_limits<double>::quiet_NaN();
      try {
        if (b2.mean() == 0 || (b1.mean() == 0 && b1.stdErr() != 0)) { ///< @todo Ok?
          // y = std::numeric_limits<double>::quiet_NaN();
          // ey = std::numeric_limits<double>::quiet_NaN();
          // throw LowStatsError("Requested division of empty bin");
        } else {
          y = b1.mean() / b2.mean();
          /// @todo Is this the exact error treatment for all (uncorrelated) cases? Behaviour around 0? +1 and -1 fills?
          const double relerr_1 = b1.stdErr() != 0 ? b1.stdErr()/b1.mean() : 0;
          const double relerr_2 = b2.stdErr() != 0 ? b2.stdErr()/b2.mean() : 0;
          ey = y * sqrt(sqr(relerr_1) + sqr(relerr_2));
        }
      } catch (const LowStatsError& e) {
        // y = std::numeric_limits<double>::quiet_NaN();
        // ey = std::numeric_limits<double>::quiet_NaN();
      }

      /// Deal with +/- errors separately, inverted for the denominator contributions:
      /// @todo check correctness with different signed numerator and denominator.
      //const double eyplus = y * sqrt( sqr(p1.yErrPlus()/p1.y()) + sqr(p2.yErrMinus()/p2.y()) );
      //const double eyminus = y * sqrt( sqr(p1.yErrMinus()/p1.y()) + sqr(p2.yErrPlus()/p2.y()) );
      rtn.addPoint(x, y, exminus, explus, ey, ey);
    }

    assert(rtn.numPoints() == numer.numBins());
    return rtn;
  }


  /// @todo Add asymm for profile histos?


}
