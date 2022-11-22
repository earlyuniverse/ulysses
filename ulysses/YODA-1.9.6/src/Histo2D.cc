// -*- C++ -*-
//
// This file is part of YODA -- Yet more Objects for Data Analysis
// Copyright (C) 2008-2021 The YODA collaboration (see AUTHORS for details)
//
#include "YODA/Histo2D.h"
#include "YODA/Profile2D.h"
#include "YODA/Scatter3D.h"
#include "YODA/Utils/StringUtils.h"

using namespace std;

namespace YODA {


  /// Copy constructor with optional new path
  Histo2D::Histo2D(const Histo2D& h, const std::string& path)
    : AnalysisObject("Histo2D", (path.size() == 0) ? h.path() : path, h, h.title()),
      _axis(h._axis)
  { }


  /// Constructor from a Scatter3D's binning, with optional new path
  Histo2D::Histo2D(const Scatter3D& s, const std::string& path)
    : AnalysisObject("Histo2D", (path.size() == 0) ? s.path() : path, s, s.title())
  {
    std::vector<HistoBin2D> bins;
    for (const Scatter3D::Point& p : s.points()) {
      bins.push_back(HistoBin2D(p.xMin(), p.xMax(), p.yMin(), p.yMax()));
    }
    _axis = Histo2DAxis(bins);
  }


  /// Constructor from a Profile2D's binning, with optional new path
  Histo2D::Histo2D(const Profile2D& p, const std::string& path)
    : AnalysisObject("Histo2D", (path.size() == 0) ? p.path() : path, p, p.title())
  {
    std::vector<HistoBin2D> bins;
    for (const ProfileBin2D& b : p.bins()) {
      bins.push_back(HistoBin2D(b.xMin(), b.xMax(), b.yMin(), b.yMax()));
    }
    _axis = Histo2DAxis(bins);
  }


  ////////////////////////////////////


  void Histo2D::fill(double x, double y, double weight, double fraction) {
    if ( std::isnan(x) ) throw RangeError("X is NaN");
    if ( std::isnan(y) ) throw RangeError("Y is NaN");

    // Fill the overall distribution
    _axis.totalDbn().fill(x, y, weight, fraction);

    // Fill the bins and overflows
    /// Unify this with Profile2D's version, when binning and inheritance are reworked
    if (inRange(x, _axis.xMin(), _axis.xMax()) && inRange(y, _axis.yMin(), _axis.yMax())) {
      try {
        /// @todo Replace try block with a check that there is a bin at x, y
        _binAt(x, y).fill(x, y, weight, fraction);
      } catch (const RangeError& re) {    }
    }
    /// @todo Reinstate! With outflow axis bin lookup
    // else {
    //   size_t ix(0), iy(0);
    //   if (x <  _axis.xMin()) ix = -1; else if (x >= _axis.xMax()) ix = 1;
    //   if (y <  _axis.yMin()) iy = -1; else if (y >= _axis.yMax()) iy = 1;
    //   _axis.outflow(ix, iy).fill(x, y, weight, fraction);
    // }

    // Lock the axis now that a fill has happened
    _axis._setLock(true);
  }


  void Histo2D::fillBin(size_t i, double weight, double fraction) {
    pair<double, double> mid = bin(i).xyMid();
    fill(mid.first, mid.second, weight, fraction);
  }



  /////////////// COMMON TO ALL BINNED

  double Histo2D::numEntries(bool includeoverflows) const {
    if (includeoverflows) return totalDbn().numEntries();
    unsigned long n = 0;
    for (const Bin& b : bins()) n += b.numEntries();
    return n;
  }


  double Histo2D::effNumEntries(bool includeoverflows) const {
    if (includeoverflows) return totalDbn().effNumEntries();
    double n = 0;
    for (const Bin& b : bins()) n += b.effNumEntries();
    return n;
  }


  double Histo2D::sumW(bool includeoverflows) const {
    if (includeoverflows) return _axis.totalDbn().sumW();
    double sumw = 0;
    for (const Bin& b : bins()) sumw += b.sumW();
    return sumw;
  }


  double Histo2D::sumW2(bool includeoverflows) const {
    if (includeoverflows) return _axis.totalDbn().sumW2();
    double sumw2 = 0;
    for (const Bin& b : bins()) sumw2 += b.sumW2();
    return sumw2;
  }


  ////////////////


  double Histo2D::xMean(bool includeoverflows) const {
    if (includeoverflows) return _axis.totalDbn().xMean();
    Dbn2D dbn;
    for (const HistoBin2D& b : bins()) dbn += b.dbn();
    return dbn.xMean();
  }


  double Histo2D::yMean(bool includeoverflows) const {
    if (includeoverflows) return _axis.totalDbn().yMean();
    Dbn2D dbn;
    for (const HistoBin2D& b : bins()) dbn += b.dbn();
    return dbn.yMean();
  }


  double Histo2D::xVariance(bool includeoverflows) const {
    if (includeoverflows) return _axis.totalDbn().xVariance();
    Dbn2D dbn;
    for (const HistoBin2D& b : bins()) dbn += b.dbn();
    return dbn.xVariance();
  }


  double Histo2D::yVariance(bool includeoverflows) const {
    if (includeoverflows) return _axis.totalDbn().yVariance();
    Dbn2D dbn;
    for (const HistoBin2D& b : bins()) dbn += b.dbn();
    return dbn.yVariance();
  }


  double Histo2D::xStdErr(bool includeoverflows) const {
    if (includeoverflows) return _axis.totalDbn().xStdErr();
    Dbn2D dbn;
    for (const HistoBin2D& b : bins()) dbn += b.dbn();
    return dbn.xStdErr();
  }


  double Histo2D::yStdErr(bool includeoverflows) const {
    if (includeoverflows) return _axis.totalDbn().yStdErr();
    Dbn2D dbn;
    for (const HistoBin2D& b : bins()) dbn += b.dbn();
    return dbn.yStdErr();
  }


  double Histo2D::xRMS(bool includeoverflows) const {
    if (includeoverflows) return _axis.totalDbn().xRMS();
    Dbn2D dbn;
    for (const HistoBin2D& b : bins()) dbn += b.dbn();
    return dbn.xRMS();
  }


  double Histo2D::yRMS(bool includeoverflows) const {
    if (includeoverflows) return _axis.totalDbn().yRMS();
    Dbn2D dbn;
    for (const HistoBin2D& b : bins()) dbn += b.dbn();
    return dbn.yRMS();
  }


  /////////////////////////////////////


  // Histo1D Histo2D::cutterX(double atY, const std::string& path, const std::string& title) {
  //   if (!_axis.isGrid()) throw GridError("Attempt to cut a Histo2D that is not a grid!");

  //   if (atY < yMin() || atY > highEdgeY()) throw RangeError("Y is outside the grid");
  //   vector<HistoBin1D> tempBins;

  //   for (double i = binByCoord(xMin(), atY).xMin(); i < highEdgeX(); i += binByCoord(i, atY).widthX()) {
  //     const HistoBin2D& b2 = binByCoord(i, atY);
  //     const Dbn1D dbn2(b2.numEntries(), b2.sumW(), b2.sumW2(), b2.sumWX(), b2.sumWX2());
  //     tempBins.push_back(HistoBin1D(b2.xMin(), b2.highEdgeX(), dbn2));
  //   }

  //   /// Setting under/over flows
  //   Dbn2D underflow;
  //   underflow += _axis.outflows()[7][_axis.getBinRow(_axis.getBinIndex(xMin(), atY))];

  //   Dbn2D overflow;
  //   overflow += _axis.outflows()[3][_axis.getBinRow(_axis.getBinIndex(xMin(), atY))];

  //   return Histo1D(tempBins, _axis.totalDbn().transformX(), underflow.transformX(), overflow.transformX(), path, title);

  // }


  // Histo1D Histo2D::cutterY(double atX, const std::string& path, const std::string& title) {
  //   if (!_axis.isGrid()) throw GridError("Attempt to cut a Histo2D that is not a grid!");

  //   if (atX < xMin() || atX > highEdgeX()) throw RangeError("X is outside the grid");
  //   vector<HistoBin1D> tempBins;

  //   for (double i = binByCoord(atX, yMin()).yMin(); i < highEdgeY(); i += binByCoord(atX, i).widthY()) {
  //     const HistoBin2D& b2 = binByCoord(atX, i);
  //     const Dbn1D dbn2(b2.numEntries(), b2.sumW(), b2.sumW2(), b2.sumWX(), b2.sumWX2());
  //     tempBins.push_back(HistoBin1D(b2.yMin(), b2.highEdgeY(), dbn2));
  //   }

  //   // Setting under/over flows
  //   Dbn2D underflow;
  //   underflow += _axis.outflows()[1][_axis.getBinColumn(_axis.getBinIndex(atX, yMin()))];

  //   Dbn2D overflow;
  //   overflow += _axis.outflows()[5][_axis.getBinColumn(_axis.getBinIndex(atX, yMin()))];
  //   Dbn2D total = _axis.totalDbn();

  //   // Making sure that we rotate our distributions, as we are cutting parallel to Y axis now
  //   total.flipXY();
  //   underflow.flipXY();
  //   overflow.flipXY();

  //   return Histo1D(tempBins, total.transformX(), underflow.transformX(), overflow.transformX(), path, title);
  // }


  // Profile1D Histo2D::mkProfileX() {
  //   if (!_axis.isGrid()) throw GridError("Profile1D cannot be made from a histogram that is not a grid!");

  //   vector<ProfileBin1D> prof;
  //   for(int i = xMin() + _axis.bin(0).xMid(); i < highEdgeX(); i+= _axis.bin(0).widthX()) {
  //     HistoBin2D& bin(_axis.binByCoord(i, yMin()));
  //     HistoBin2D composite(bin.xMin(), bin.xMax(), bin.yMin(), bin.yMax()) ;
  //     for(int j = yMin() + _axis.bin(0).yMid(); j < highEdgeY(); j += _axis.bin(0).widthY()) {
  //       composite += _axis.binByCoord(i, j);
  //     }
  //     prof.push_back(composite.transformX());
  //   }

  //   vector<vector<Dbn2D> >& outflows = _axis.outflows();

  //   /// Properly setting an underflow
  //   Dbn2D underflow;
  //   underflow += outflows[0][0]; underflow += outflows[6][0];
  //   for(size_t i = 0; i < outflows[7].size(); ++i) {
  //     underflow += outflows[7][i];
  //   }

  //   /// Setting an overflow
  //   Dbn2D overflow;
  //   overflow += outflows[2][0]; overflow += outflows[4][0];
  //   for(size_t i = 0; i < outflows[3].size(); ++i) {
  //     overflow += outflows[3][i];
  //   }

  //   /// And constructing a profile 1D from all this data
  //   Profile1D ret(prof, _axis.totalDbn(), underflow, overflow);
  //   return ret;

  // }

  // Profile1D Histo2D::mkProfileY() {
  //   if (!_axis.isGrid()) throw GridError("Profile1D cannot be made from a histogram that is not a grid!");

  //   vector<ProfileBin1D> prof;
  //   for(int i = yMin() + _axis.bin(0).yMid(); i < highEdgeY(); i+= _axis.bin(0).widthY()) {
  //     HistoBin2D& bin(_axis.binByCoord(i, yMin()));
  //     HistoBin2D composite(bin.xMin(), bin.xMax(), bin.yMin(), bin.yMax()) ;
  //     for(int j = xMin() + _axis.bin(0).xMid(); j < highEdgeX(); j += _axis.bin(0).widthX()) {
  //       composite += _axis.binByCoord(i, j);
  //     }
  //     prof.push_back(composite.transformY());
  //   }

  //   vector<vector<Dbn2D> >& outflows = _axis.outflows();

  //   /// Properly setting an underflow
  //   Dbn2D underflow;
  //   underflow += outflows[0][0]; underflow += outflows[2][0];
  //   for(size_t i = 0; i < outflows[1].size(); ++i) {
  //     underflow += outflows[1][i];
  //   }

  //   /// Setting an overflow
  //   Dbn2D overflow;
  //   overflow += outflows[6][0]; overflow += outflows[4][0];
  //   for(size_t i = 0; i < outflows[5].size(); ++i) {
  //     overflow += outflows[5][i];
  //   }

  //   /// Setting a flipped total distribution
  //   Dbn2D td = _axis.totalDbn();
  //   td.flipXY();

  //   /// And constructing a profile 1D from all this data
  //   Profile1D ret(prof, td, underflow, overflow);
  //   return ret;
  // }


  Scatter3D divide(const Histo2D& numer, const Histo2D& denom) {
    Scatter3D rtn;

    for (size_t i = 0; i < numer.numBins(); ++i) {
      const HistoBin2D& b1 = numer.bin(i);
      const HistoBin2D& b2 = denom.bin(i);

      /// @todo Create a compatibleBinning function? Or just compare vectors of edges().
      if (!fuzzyEquals(b1.xMin(), b2.xMin()) || !fuzzyEquals(b1.xMax(), b2.xMax()))
        throw BinningError("x binnings are not equivalent in " + numer.path() + " / " + denom.path());
      if (!fuzzyEquals(b1.yMin(), b2.yMin()) || !fuzzyEquals(b1.yMax(), b2.yMax()))
        throw BinningError("y binnings are not equivalent in " + numer.path() + " / " + denom.path());

      // Assemble the x value and error
      // Use the midpoint of the "bin" for the new central x value, in the absence of better information
      const double x = b1.xMid();
      const double exminus = x - b1.xMin();
      const double explus  = b1.xMax() - x;

      // Assemble the y value and error
      // Use the midpoint of the "bin" for the new central y value, in the absence of better information
      const double y = b1.yMid();
      const double eyminus = y - b1.yMin();
      const double eyplus  = b1.yMax() - y;

      // Assemble the z value and error
      double z = std::numeric_limits<double>::quiet_NaN();
      double ez = std::numeric_limits<double>::quiet_NaN();
      if (b2.height() == 0 || (b1.height() == 0 && b1.heightErr() != 0)) { ///< @todo Ok?
        // z = std::numeric_limits<double>::quiet_NaN();
        // ez = std::numeric_limits<double>::quiet_NaN();
        // throw LowStatsError("Requested division of empty bin");
      } else {
        z = b1.height() / b2.height();
        /// @todo Is this the exact error treatment for all (uncorrelated) cases? Behaviour around 0? +1 and -1 fills?
        const double relerr_1 = b1.heightErr() != 0 ? b1.relErr() : 0;
        const double relerr_2 = b2.heightErr() != 0 ? b2.relErr() : 0;
        ez = z * sqrt(sqr(relerr_1) + sqr(relerr_2));
      }
      /// Deal with +/- errors separately, inverted for the denominator contributions:
      /// @todo check correctness with different signed numerator and denominator.
      //const double eyplus = y * sqrt( sqr(p1.yErrPlus()/p1.y()) + sqr(p2.yErrMinus()/p2.y()) );
      //const double eyminus = y * sqrt( sqr(p1.yErrMinus()/p1.y()) + sqr(p2.yErrPlus()/p2.y()) );
      rtn.addPoint(x, y, z, exminus, explus, eyminus, eyplus, ez, ez);
    }

    assert(rtn.numPoints() == numer.numBins());
    return rtn;
  }


  Scatter3D efficiency(const Histo2D& accepted, const Histo2D& total) {
    Scatter3D tmp = divide(accepted, total);
    for (size_t i = 0; i < accepted.numBins(); ++i) {
      const HistoBin2D& b_acc = accepted.bin(i);
      const HistoBin2D& b_tot = total.bin(i);
      Point3D& point = tmp.point(i);

      /// BEGIN DIMENSIONALITY-INDEPENDENT BIT TO SHARE WITH H1

      // Check that the numerator is consistent with being a subset of the denominator
      /// @note Neither effNumEntries nor sumW are guaranteed to satisfy num <= den for general weights!
      if (b_acc.numEntries() > b_tot.numEntries())
        throw UserError("Attempt to calculate an efficiency when the numerator is not a subset of the denominator: "
                        + Utils::toStr(b_acc.numEntries()) + " entries / " + Utils::toStr(b_tot.numEntries()) + " entries");

      // If no entries on the denominator, set eff = err = 0 and move to the next bin
      double eff = std::numeric_limits<double>::quiet_NaN();
      double err = std::numeric_limits<double>::quiet_NaN();
      try {
        if (b_tot.sumW() != 0) {
          eff = b_acc.sumW() / b_tot.sumW(); //< Actually this is already calculated by the division...
          err = sqrt(abs( ((1-2*eff)*b_acc.sumW2() + sqr(eff)*b_tot.sumW2()) / sqr(b_tot.sumW()) ));
        }
      } catch (const LowStatsError& e) {
        //
      }

      /// END DIMENSIONALITY-INDEPENDENT BIT TO SHARE WITH H1

      point.setZ(eff, err);
    }
    return tmp;

  }


  /// @todo Add asymm for 2D histos?


}
