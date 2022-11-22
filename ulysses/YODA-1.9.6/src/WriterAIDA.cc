// -*- C++ -*-
//
// This file is part of YODA -- Yet more Objects for Data Analysis
// Copyright (C) 2008-2021 The YODA collaboration (see AUTHORS for details)
//
#include "YODA/WriterAIDA.h"
#include "YODA/Utils/StringUtils.h"

#include <iostream>
#include <iomanip>

using namespace std;

namespace YODA {

  /// Singleton creation function
  Writer& WriterAIDA::create() {
    static WriterAIDA _instance;
    _instance.setPrecision(6);
    return _instance;
  }


  void WriterAIDA::writeHead(std::ostream& stream) {
    stream << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    stream << "<!DOCTYPE aida SYSTEM \"http://aida.freehep.org/schemas/3.0/aida.dtd\">\n";
    stream << "<aida>\n";
    stream << "  <implementation version=\"1.0\" package=\"YODA\"/>\n";
  }


  void WriterAIDA::writeFoot(std::ostream& stream) {
    stream << "</aida>\n" << flush;
  }


  void WriterAIDA::writeCounter(std::ostream& os, const Counter&) {
    os << endl << "<!-- COUNTER WRITING TO AIDA IS CURRENTLY UNSUPPORTED! -->" << endl << endl;
  }


  void WriterAIDA::writeHisto1D(std::ostream& os, const Histo1D& h) {
    Scatter2D tmp = mkScatter(h);
    tmp.setAnnotation("Type", "Histo1D");
    writeScatter2D(os, tmp);
  }


  void WriterAIDA::writeHisto2D(std::ostream& os, const Histo2D&) {
    os << endl << "<!-- HISTO2D WRITING TO AIDA IS CURRENTLY UNSUPPORTED! -->" << endl << endl;
    // Scatter3D tmp = mkScatter(h);
    // tmp.setAnnotation("Type", "Histo2D");
    // writeScatter3D(os, tmp);
  }


  void WriterAIDA::writeProfile1D(std::ostream& os, const Profile1D& p) {
    Scatter2D tmp = mkScatter(p);
    tmp.setAnnotation("Type", "Profile1D");
    writeScatter2D(os, tmp);
  }


  void WriterAIDA::writeProfile2D(std::ostream& os, const Profile2D& p) {
    os << endl << "<!-- PROFILE2D WRITING TO AIDA IS CURRENTLY UNSUPPORTED! -->" << endl << endl;
    // Scatter3D tmp = mkScatter(p);
    // tmp.setAnnotation("Type", "Profile2D");
    // writeScatter3D(os, tmp);
  }


  void WriterAIDA::writeScatter1D(std::ostream& os, const Scatter1D& s) {
    os << endl << "<!-- SCATTER1D WRITING TO AIDA IS CURRENTLY UNSUPPORTED! -->" << endl << endl;
  }


  void WriterAIDA::writeScatter2D(std::ostream& os, const Scatter2D& s) {
    ios_base::fmtflags oldflags = os.flags();
    // const int precision = 8;
    os << scientific << showpoint << setprecision(_precision);

    string name = "";
    string path = "/";
    const size_t slashpos = s.path().rfind("/");
    if (slashpos != string::npos) {
      name = s.path().substr(slashpos+1, s.path().length() - slashpos - 1);
      if (slashpos > 0) path = s.path().substr(0, slashpos);
    }
    os << "  <dataPointSet name=\"" << Utils::encodeForXML(name) << "\"\n"
       << "    title=\"" << Utils::encodeForXML(s.title()) << "\""
       << " path=\"" << Utils::encodeForXML(path) << "\" dimension=\"2\">\n";
    os << "    <dimension dim=\"0\" title=\"\" />\n";
    os << "    <dimension dim=\"1\" title=\"\" />\n";
    os << "    <annotation>\n";
    for (const string& a : s.annotations()) {
      if (a.empty()) continue;
      os << "      <item key=\"" << Utils::encodeForXML(a)
         << "\" value=\"" << Utils::encodeForXML(s.annotation(a)) << "\" />\n";
    }
    if (!s.hasAnnotation("Type")) {
      os << "      <item key=\"Type\" value=\"Scatter2D\" />\n";
    }
    os << "    </annotation>\n";
    for (const Point2D& pt : s.points()) {
      os << "    <dataPoint>\n";
      os << "      <measurement value=\"" << pt.x()
	 << "\" errorPlus=\"" << pt.xErrPlus()
         << "\" errorMinus=\"" << pt.xErrMinus()
	 << "\"/>\n";
      os << "      <measurement value=\"" << pt.y()
	 << "\" errorPlus=\"" << pt.yErrPlus()
         << "\" errorMinus=\"" << pt.yErrMinus()
	 << "\"/>\n";
      os << "    </dataPoint>\n";
    }
    os << "  </dataPointSet>\n";
    os << flush;

    os.flags(oldflags);
  }


  void WriterAIDA::writeScatter3D(std::ostream& os, const Scatter3D& s) {
    os << endl << "<!-- SCATTER3D WRITING TO AIDA IS CURRENTLY UNSUPPORTED! -->" << endl << endl;
  }


}
