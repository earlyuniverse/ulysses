// -*- C++ -*-
//
// This file is part of YODA -- Yet more Objects for Data Analysis
// Copyright (C) 2008-2021 The YODA collaboration (see AUTHORS for details)
//
#include "YODA/WriterFLAT.h"

#include "YODA/Counter.h"
#include "YODA/Histo1D.h"
#include "YODA/Histo2D.h"
#include "YODA/Profile1D.h"
#include "YODA/Profile2D.h"
#include "YODA/Scatter1D.h"
#include "YODA/Scatter2D.h"
#include "YODA/Scatter3D.h"

#include <iostream>
#include <iomanip>

using namespace std;

namespace YODA {

  /// Singleton creation function
  Writer& WriterFLAT::create() {
    static WriterFLAT _instance;
    _instance.setPrecision(6);
    return _instance;
  }

  void WriterFLAT::_writeAnnotations(std::ostream& os, const AnalysisObject& ao) {
    os << scientific << setprecision(_aoprecision);
    for (const string& a : ao.annotations()) {
      if (a.empty()) continue;
      if (a == "Type") continue;
      /// @todo Should write out floating point annotations as scientific notation...
      os << a << "=" << ao.annotation(a) << "\n";
    }
  }


  void WriterFLAT::writeCounter(std::ostream& os, const Counter& c) {
    ios_base::fmtflags oldflags = os.flags();
    os << scientific << showpoint << setprecision(_aoprecision);

    os << "# BEGIN COUNTER " << c.path() << "\n";
    _writeAnnotations(os, c);
    os << "# value\t error\n";
    os << c.val() << "\t" << c.err() << "\n";
    os << "# END COUNTER\n\n";

    os << flush;
    os.flags(oldflags);
  }


  void WriterFLAT::writeHisto1D(std::ostream& os, const Histo1D& h) {
    Scatter2D tmp = mkScatter(h);
    tmp.setAnnotation("Type", "Histo1D");
    writeScatter2D(os, tmp);
  }


  void WriterFLAT::writeHisto2D(std::ostream& os, const Histo2D& h) {
    Scatter3D tmp = mkScatter(h);
    tmp.setAnnotation("Type", "Histo2D");
    writeScatter3D(os, tmp);
  }


  void WriterFLAT::writeProfile1D(std::ostream& os, const Profile1D& p) {
    Scatter2D tmp = mkScatter(p);
    tmp.setAnnotation("Type", "Profile1D");
    writeScatter2D(os, tmp);
  }


  void WriterFLAT::writeProfile2D(std::ostream& os, const Profile2D& h) {
    Scatter3D tmp = mkScatter(h);
    tmp.setAnnotation("Type", "Profile2D");
    writeScatter3D(os, tmp);
  }




  void WriterFLAT::writeScatter1D(std::ostream& os, const Scatter1D& s) {
    ios_base::fmtflags oldflags = os.flags();
    os << scientific << showpoint << setprecision(_aoprecision);

    os << "# BEGIN VALUE " << s.path() << "\n";
    _writeAnnotations(os, s);
    os << "# value\t errminus\t errplus\n";
    for (const Point1D& pt : s.points()) {
      os << pt.x() << "\t" << pt.xErrMinus() << "\t" << pt.xErrPlus() << "\n";
    }
    os << "# END VALUE\n\n";

    os << flush;
    os.flags(oldflags);
  }


  void WriterFLAT::writeScatter2D(std::ostream& os, const Scatter2D& s) {
    ios_base::fmtflags oldflags = os.flags();
    os << scientific << showpoint << setprecision(_aoprecision);

    os << "# BEGIN HISTO1D " << s.path() << "\n";
    _writeAnnotations(os, s);
    os << "# xlow\t xhigh\t val\t errminus\t errplus\n";
    for (const Point2D& pt : s.points()) {
      os << pt.x()-pt.xErrMinus() << "\t" << pt.x()+pt.xErrPlus() << "\t";
      os << pt.y() << "\t" << pt.yErrMinus() << "\t" << pt.yErrPlus() << "\n";
    }
    os << "# END HISTO1D\n\n";

    os << flush;
    os.flags(oldflags);
  }


  void WriterFLAT::writeScatter3D(std::ostream& os, const Scatter3D& s) { // , bool asHist2D) {
    ios_base::fmtflags oldflags = os.flags();
    os << scientific << showpoint << setprecision(_aoprecision);

    os << "# BEGIN HISTO2D " << s.path() << "\n";
    _writeAnnotations(os, s);

    // if (asHist2D) { // Extension of what writeScatter2D does
    os << "# xlow\t xhigh\t ylow\t yhigh\t val\t errminus\t errplus\n";
    for (const Point3D& pt : s.points()) {
      os << pt.x()-pt.xErrMinus() << "\t" << pt.x()+pt.xErrPlus() << "\t";
      os << pt.y()-pt.yErrMinus() << "\t" << pt.y()+pt.yErrPlus() << "\t";
      os << pt.z() << "\t" << pt.zErrMinus() << "\t" << pt.zErrPlus() << "\n";
    }

    // } else { // What writerYODA should do... let's just put this in there (generalised for multiple errs).
    //   os << "# xval\t xerr-\t xerr+\t yval\t yerr-\t yerr+\t zval\t zerr-\t zerr+\n";
    //   for (const Point3D& pt : s.points()) {
    //     os << pt.x() << "\t" << pt.xErrMinus() << "\t" << pt.xErrMinus() << "\t";
    //     os << pt.y() << "\t" << pt.yErrMinus() << "\t" << pt.yErrMinus() << "\t";
    //     os << pt.z() << "\t" << pt.zErrMinus() << "\t" << pt.zErrMinus() << "\n";
    //   }
    // }

    os << "# END HISTO2D\n\n";

    os << flush;
    os.flags(oldflags);
  }


}
