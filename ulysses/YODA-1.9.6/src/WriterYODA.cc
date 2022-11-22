// -*- C++ -*-
//
// This file is part of YODA -- Yet more Objects for Data Analysis
// Copyright (C) 2008-2021 The YODA collaboration (see AUTHORS for details)
//
#include "YODA/WriterYODA.h"

#include "yaml-cpp/yaml.h"
#ifdef YAML_NAMESPACE
#define YAML YAML_NAMESPACE
#endif

#include <iostream>
#include <iomanip>
using namespace std;

namespace YODA {


  /// Singleton creation function
  Writer& WriterYODA::create() {
    static WriterYODA _instance;
    _instance.setPrecision(6);
    return _instance;
  }


  /// YODA text-format version
  ///
  /// - V1/empty = make-plots annotations style
  /// - V2 = YAML annotations
  static const int YODA_FORMAT_VERSION = 2;


  // Version-formatting helper function
  inline string _iotypestr(const string& baseiotype) {
    ostringstream os;
    os << "YODA_" << Utils::toUpper(baseiotype) << "_V" << YODA_FORMAT_VERSION;
    return os.str();
  }


  void WriterYODA::_writeAnnotations(std::ostream& os, const AnalysisObject& ao) {
    os << scientific << setprecision(_aoprecision);
    for (const string& a : ao.annotations()) {
      if (a.empty()) continue;
      /// @todo Write out floating point annotations as scientific notation
      string ann = ao.annotation(a);
      // remove stpurious line returns at the end of a string so that we don't
      // end up with two line returns.
      ann.erase(std::remove(ann.begin(), ann.end(), '\n'), ann.end());
      os << a << ": " << ann << "\n";
    }
    os << "---\n";
  }


  void WriterYODA::writeCounter(std::ostream& os, const Counter& c) {
    ios_base::fmtflags oldflags = os.flags();
    os << scientific << showpoint << setprecision(_aoprecision);

    os << "BEGIN " << _iotypestr("COUNTER") << " " << c.path() << "\n";
    _writeAnnotations(os, c);
    os << "# sumW\t sumW2\t numEntries\n";
    os << c.sumW()  << "\t" << c.sumW2() << "\t" << c.numEntries() << "\n";
    os << "END " << _iotypestr("COUNTER") << "\n\n";

    os.flags(oldflags);
  }


  void WriterYODA::writeHisto1D(std::ostream& os, const Histo1D& h) {
    ios_base::fmtflags oldflags = os.flags();
    os << scientific << showpoint << setprecision(_aoprecision);

    os << "BEGIN " << _iotypestr("HISTO1D") << " " << h.path() << "\n";
    _writeAnnotations(os, h);
    try {
      //if ( h.totalDbn().effNumEntries() > 0 ) {
      os << "# Mean: " << h.xMean() << "\n";
      os << "# Area: " << h.integral() << "\n";
    } catch (LowStatsError& e) {
      //
    }
    os << "# ID\t ID\t sumw\t sumw2\t sumwx\t sumwx2\t numEntries\n";
    os << "Total   \tTotal   \t";
    os << h.totalDbn().sumW()  << "\t" << h.totalDbn().sumW2()  << "\t";
    os << h.totalDbn().sumWX() << "\t" << h.totalDbn().sumWX2() << "\t";
    os << h.totalDbn().numEntries() << "\n";
    os << "Underflow\tUnderflow\t";
    os << h.underflow().sumW()  << "\t" << h.underflow().sumW2()  << "\t";
    os << h.underflow().sumWX() << "\t" << h.underflow().sumWX2() << "\t";
    os << h.underflow().numEntries() << "\n";
    os << "Overflow\tOverflow\t";
    os << h.overflow().sumW()  << "\t" << h.overflow().sumW2()  << "\t";
    os << h.overflow().sumWX() << "\t" << h.overflow().sumWX2() << "\t";
    os << h.overflow().numEntries() << "\n";
    os << "# xlow\t xhigh\t sumw\t sumw2\t sumwx\t sumwx2\t numEntries\n";
    for (const HistoBin1D& b : h.bins()) {
      os << b.xMin() << "\t" << b.xMax() << "\t";
      os << b.sumW()    << "\t" << b.sumW2()    << "\t";
      os << b.sumWX()   << "\t" << b.sumWX2()   << "\t";
      os << b.numEntries() << "\n";
    }
    os << "END " << _iotypestr("HISTO1D") << "\n\n";

    os.flags(oldflags);
  }


  void WriterYODA::writeHisto2D(std::ostream& os, const Histo2D& h) {
    ios_base::fmtflags oldflags = os.flags();
    os << scientific << showpoint << setprecision(_aoprecision);
    os << "BEGIN " << _iotypestr("HISTO2D") << " " << h.path() << "\n";
    _writeAnnotations(os, h);
    try {
      //if ( h.totalDbn().numEntries() > 0 )
      os << "# Mean: (" << h.xMean() << ", " << h.yMean() << ")\n";
      os << "# Volume: " << h.integral() << "\n";
    } catch (LowStatsError& e) {
      //
    }
    os << "# ID\t ID\t sumw\t sumw2\t sumwx\t sumwx2\t sumwy\t sumwy2\t sumwxy\t numEntries\n";
    // Total distribution
    const Dbn2D& td = h.totalDbn();
    os << "Total   \tTotal   \t";
    os << td.sumW()   << "\t" << td.sumW2()  << "\t";
    os << td.sumWX()  << "\t" << td.sumWX2() << "\t";
    os << td.sumWY()  << "\t" << td.sumWY2() << "\t";
    os << td.sumWXY() << "\t";
    os << td.numEntries() << "\n";
    // Outflows
    /// @todo Disabled for now, reinstate with a *full* set of outflow info to allow marginalisation
    os << "# 2D outflow persistency not currently supported until API is stable\n";
    // for (int ix = -1; ix <= 1; ++ix) {
    //   for (int iy = -1; iy <= 1; ++iy) {
    //     if (ix == 0 && iy == 0) continue;
    //     os << "Outflow\t" << ix << ":" << iy << "\t";
    //     const Dbn2D& d = h.outflow(ix, iy);
    //     os << d.sumW()   << "\t" << d.sumW2()  << "\t";
    //     os << d.sumWX()  << "\t" << d.sumWX2() << "\t";
    //     os << d.sumWY()  << "\t" << d.sumWY2() << "\t";
    //     os << d.sumWXY() << "\t";
    //     os << d.numEntries() << "\n";
    //   }
    // }
    // Bins
    os << "# xlow\t xhigh\t ylow\t yhigh\t sumw\t sumw2\t sumwx\t sumwx2\t sumwy\t sumwy2\t sumwxy\t numEntries\n";
    for (const HistoBin2D& b : h.bins()) {
      os << b.xMin() << "\t" << b.xMax() << "\t";
      os << b.yMin() << "\t" << b.yMax() << "\t";
      os << b.sumW()     << "\t" << b.sumW2()     << "\t";
      os << b.sumWX()    << "\t" << b.sumWX2()    << "\t";
      os << b.sumWY()    << "\t" << b.sumWY2()    << "\t";
      os << b.sumWXY()   << "\t";
      os << b.numEntries() << "\n";
    }
    os << "END " << _iotypestr("HISTO2D") << "\n\n";

    os.flags(oldflags);
  }


  void WriterYODA::writeProfile1D(std::ostream& os, const Profile1D& p) {
    ios_base::fmtflags oldflags = os.flags();
    os << scientific << showpoint << setprecision(_aoprecision);

    os << "BEGIN " << _iotypestr("PROFILE1D") << " " << p.path() << "\n";
    _writeAnnotations(os, p);
    os << "# ID\t ID\t sumw\t sumw2\t sumwx\t sumwx2\t sumwy\t sumwy2\t numEntries\n";
    os << "Total   \tTotal   \t";
    os << p.totalDbn().sumW()  << "\t" << p.totalDbn().sumW2()  << "\t";
    os << p.totalDbn().sumWX() << "\t" << p.totalDbn().sumWX2() << "\t";
    os << p.totalDbn().sumWY() << "\t" << p.totalDbn().sumWY2() << "\t";
    os << p.totalDbn().numEntries() << "\n";
    os << "Underflow\tUnderflow\t";
    os << p.underflow().sumW()  << "\t" << p.underflow().sumW2()  << "\t";
    os << p.underflow().sumWX() << "\t" << p.underflow().sumWX2() << "\t";
    os << p.underflow().sumWY() << "\t" << p.underflow().sumWY2() << "\t";
    os << p.underflow().numEntries() << "\n";
    os << "Overflow\tOverflow\t";
    os << p.overflow().sumW()  << "\t" << p.overflow().sumW2()  << "\t";
    os << p.overflow().sumWX() << "\t" << p.overflow().sumWX2() << "\t";
    os << p.overflow().sumWY() << "\t" << p.overflow().sumWY2() << "\t";
    os << p.overflow().numEntries() << "\n";
    os << "# xlow\t xhigh\t sumw\t sumw2\t sumwx\t sumwx2\t sumwy\t sumwy2\t numEntries\n";
    for (const ProfileBin1D& b : p.bins()) {
      os << b.xMin() << "\t" << b.xMax() << "\t";
      os << b.sumW()    << "\t" << b.sumW2()    << "\t";
      os << b.sumWX()   << "\t" << b.sumWX2()   << "\t";
      os << b.sumWY()   << "\t" << b.sumWY2()   << "\t";
      os << b.numEntries() << "\n";
    }
    os << "END " << _iotypestr("PROFILE1D") << "\n\n";

    os.flags(oldflags);
  }


  void WriterYODA::writeProfile2D(std::ostream& os, const Profile2D& p) {
    ios_base::fmtflags oldflags = os.flags();
    os << scientific << showpoint << setprecision(_aoprecision);

    os << "BEGIN " << _iotypestr("PROFILE2D") << " " << p.path() << "\n";
    _writeAnnotations(os, p);
    os << "# sumw\t sumw2\t sumwx\t sumwx2\t sumwy\t sumwy2\t sumwz\t sumwz2\t sumwxy\t numEntries\n";
    // Total distribution
    const Dbn3D& td = p.totalDbn();
    os << "Total   \tTotal   \t";
    os << td.sumW()   << "\t" << td.sumW2()  << "\t";
    os << td.sumWX()  << "\t" << td.sumWX2() << "\t";
    os << td.sumWY()  << "\t" << td.sumWY2() << "\t";
    os << td.sumWZ()  << "\t" << td.sumWZ2() << "\t";
    os << td.sumWXY() << "\t"; // << td.sumWXZ() << "\t" << td.sumWYZ() << "\t";
    os << td.numEntries() << "\n";
    // Outflows
    /// @todo Disabled for now, reinstate with a *full* set of outflow info to allow marginalisation
    os << "# 2D outflow persistency not currently supported until API is stable\n";
    // for (int ix = -1; ix <= 1; ++ix) {
    //   for (int iy = -1; iy <= 1; ++iy) {
    //     if (ix == 0 && iy == 0) continue;
    //     os << "Outflow\t" << ix << ":" << iy << "\t";
    //     const Dbn3D& d = p.outflow(ix, iy);
    //     os << d.sumW()   << "\t" << d.sumW2()  << "\t";
    //     os << d.sumWX()  << "\t" << d.sumWX2() << "\t";
    //     os << d.sumWY()  << "\t" << d.sumWY2() << "\t";
    //     os << d.sumWZ()  << "\t" << d.sumWZ2() << "\t";
    //     os << d.sumWXY() << "\t"; // << d.sumWXZ() << "\t" << d.sumWYZ() << "\t";
    //     os << d.numEntries() << "\n";
    //   }
    // }
    // Bins
    os << "# xlow\t xhigh\t ylow\t yhigh\t sumw\t sumw2\t sumwx\t sumwx2\t sumwy\t sumwy2\t sumwz\t sumwz2\t sumwxy\t numEntries\n";
    for (const ProfileBin2D& b : p.bins()) {
      os << b.xMin() << "\t" << b.xMax() << "\t";
      os << b.yMin() << "\t" << b.yMax() << "\t";
      os << b.sumW()     << "\t" << b.sumW2()     << "\t";
      os << b.sumWX()    << "\t" << b.sumWX2()    << "\t";
      os << b.sumWY()    << "\t" << b.sumWY2()    << "\t";
      os << b.sumWZ()    << "\t" << b.sumWZ2()    << "\t";
      os << b.sumWXY()   << "\t"; // << b.sumWXZ()    << "\t" << b.sumWYZ() << "\t";
      os << b.numEntries() << "\n";
    }
    os << "END " << _iotypestr("PROFILE2D") << "\n\n";

    os.flags(oldflags);
  }


  void WriterYODA::writeScatter1D(std::ostream& os, const Scatter1D& s) {
    ios_base::fmtflags oldflags = os.flags();
    os << scientific << showpoint << setprecision(_aoprecision);

    // we promised not to modify const s, but we want to add an annotation
    // we did not promise to not modify the *clone* of s...
    auto sclone =  s.clone();
    sclone.writeVariationsToAnnotations();

    os << "BEGIN " << _iotypestr("SCATTER1D") << " " << s.path() << "\n";
    _writeAnnotations(os, sclone);

    //write headers
    std::string headers="# xval\t xerr-\t xerr+\t";
    os << headers << "\n";

    //write points
    for (const Point1D& pt : s.points()) {
      // fill central value
      os << pt.x() << "\t" << pt.xErrMinus() << "\t" << pt.xErrPlus() ;
      os <<  "\n";
    }
    os << "END " << _iotypestr("SCATTER1D") << "\n\n";

    os << flush;
    os.flags(oldflags);
  }


  void WriterYODA::writeScatter2D(std::ostream& os, const Scatter2D& s) {
    ios_base::fmtflags oldflags = os.flags();
    os << scientific << showpoint << setprecision(_aoprecision);
    os << "BEGIN " << _iotypestr("SCATTER2D") << " " << s.path() << "\n";

    // Write annotations.
    // We promised not to modify const s, but we want to add an annotation;
    // We did not promise to not modify the *clone* of s;
    // Judge not, lest ye be judged
    auto sclone = s.clone();
    sclone.writeVariationsToAnnotations();
    _writeAnnotations(os, sclone);

    //write headers
    /// @todo Change ordering to {vals} {errs} {errs} ...
    std::string headers="# xval\t xerr-\t xerr+\t yval\t yerr-\t yerr+\t";
    os << headers << "\n";

    //write points
    for (const Point2D& pt : s.points()) {
      /// @todo Change ordering to {vals} {errs} {errs} ...
      // fill central value
      os << pt.x() << "\t" << pt.xErrMinus() << "\t" << pt.xErrPlus() << "\t";
      os << pt.y() << "\t" << pt.yErrMinus() << "\t" << pt.yErrPlus() ;
      os <<  "\n";
    }
    os << "END " << _iotypestr("SCATTER2D") << "\n\n";

    os << flush;
    os.flags(oldflags);
  }


  void WriterYODA::writeScatter3D(std::ostream& os, const Scatter3D& s) {
    ios_base::fmtflags oldflags = os.flags();
    os << scientific << showpoint << setprecision(_aoprecision);
    os << "BEGIN " << _iotypestr("SCATTER3D") << " " << s.path() << "\n";

    // write annotations
    // we promised not to modify const s, but we want to add an annotation
    // we did not promise to not modify the *clone* of s...
    auto sclone =  s.clone();
    sclone.writeVariationsToAnnotations();
    _writeAnnotations(os, sclone);

   // std::vector<std::string> variations= s.variations();
    //write headers
    /// @todo Change ordering to {vals} {errs} {errs} ...
    std::string headers="# xval\t xerr-\t xerr+\t yval\t yerr-\t yerr+\t zval\t zerr-\t zerr+\t";
    os << headers << "\n";

    //write points
    for (const Point3D& pt : s.points()) {
      /// @todo Change ordering to {vals} {errs} {errs} ...
      // fill central value
      os << pt.x() << "\t" << pt.xErrMinus() << "\t" << pt.xErrPlus() << "\t";
      os << pt.y() << "\t" << pt.yErrMinus() << "\t" << pt.yErrPlus() << "\t";
      os << pt.z() << "\t" << pt.zErrMinus() << "\t" << pt.zErrPlus() ;
      os <<  "\n";
    }
    os << "END " << _iotypestr("SCATTER3D") << "\n\n";

    os << flush;
    os.flags(oldflags);
  }


}
