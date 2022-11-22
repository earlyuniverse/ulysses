// -*- C++ -*-
//
// This file is part of YODA -- Yet more Objects for Data Analysis
// Copyright (C) 2008-2021 The YODA collaboration (see AUTHORS for details)
//
#include "YODA/Writer.h"
#include "YODA/WriterYODA.h"
#include "YODA/WriterAIDA.h"
#include "YODA/WriterFLAT.h"
#include "YODA/Config/BuildConfig.h"

#ifdef HAVE_LIBZ
#define _XOPEN_SOURCE 700
#include "zstr/zstr.hpp"
#endif

#include <iostream>
#include <locale>
#include <typeinfo>
#include <sstream>
using namespace std;

namespace YODA {


  Writer& mkWriter(const string& name) {
    // Determine the format from the string (a file or file extension)
    const size_t lastdot = name.find_last_of(".");
    string fmt = Utils::toLower(lastdot == string::npos ? name : name.substr(lastdot+1));
    const bool compress = (fmt == "gz");
    //cout << "***" << compress << endl;
    if (compress) {
      #ifndef HAVE_LIBZ
      throw UserError("YODA was compiled without zlib support: can't write " + name);
      #endif
      const size_t lastbutonedot = (lastdot == string::npos) ? string::npos : name.find_last_of(".", lastdot-1);
      fmt = Utils::toLower(lastbutonedot == string::npos ? name : name.substr(lastbutonedot+1));
    }
    // Create the appropriate Writer
    Writer* w = nullptr;
    if (Utils::startswith(fmt, "yoda")) w = &WriterYODA::create();
    if (Utils::startswith(fmt, "aida")) w = &WriterAIDA::create();
    if (Utils::startswith(fmt, "dat" )) w = &WriterFLAT::create(); ///< @todo Improve/remove... .ydat?
    if (Utils::startswith(fmt, "flat")) w = &WriterFLAT::create();
    if (!w) throw UserError("Format cannot be identified from string '" + name + "'");
    w->useCompression(compress);
    return *w;
  }


  void Writer::write(const std::string& filename, const AnalysisObject& ao) {
    std::vector<const AnalysisObject*> vec{&ao};
    write(filename, vec);
  }

  // Canonical writer function, including compression handling
  void Writer::write(ostream& stream, const vector<const AnalysisObject*>& aos) {
    std::unique_ptr<std::ostream> zos;
    std::ostream* os = &stream;

    // Write numbers in the "C" locale
    std::locale prev_locale = os->getloc();
    os->imbue(std::locale::classic());

    // Wrap the stream if needed
    if (_compress) {
      #ifdef HAVE_LIBZ
      // Doesn't work to always create zstr wrapper: have to only create if compressing :-/
      // zstr::ostream zstream(stream);
      // ostream& os = _compress ? zstream : stream;
      os = new zstr::ostream(stream);
      zos.reset(os);
      #else
      throw UserError("YODA was compiled without zlib support: can't write to a compressed stream");
      #endif
    }

    // Write the data components
    /// @todo Remove the head/body/foot distinction?
    writeHead(*os);
    bool first = true;
    for (const AnalysisObject* aoptr : aos) {
      setAOPrecision( aoptr->annotation("WriterDoublePrecision", 0) );
      try {
        if (!first) *os << "\n"; //< blank line between items
        writeBody(*os, aoptr);
        first = false;
      } catch (const LowStatsError& ex) {
        /// @todo Why specifically LowStatsError?
        std::cerr << "LowStatsError in writing AnalysisObject " << aoptr->title() << ":\n" << ex.what() << "\n";
      }
    }
    writeFoot(*os);
    *os << flush;

    os->imbue(prev_locale);
  }


  void Writer::writeBody(ostream& stream, const AnalysisObject* ao) {
    if (!ao) throw WriteError("Attempting to write a null AnalysisObject*");
    writeBody(stream, *ao);
  }

  void Writer::writeBody(ostream& stream, const AnalysisObject& ao) {
    const string aotype = ao.type();
    if (aotype == "Counter") {
      writeCounter(stream, dynamic_cast<const Counter&>(ao));
    } else if (aotype == "Histo1D") {
      writeHisto1D(stream, dynamic_cast<const Histo1D&>(ao));
    } else if (aotype == "Histo2D") {
      writeHisto2D(stream, dynamic_cast<const Histo2D&>(ao));
    } else if (aotype == "Profile1D") {
      writeProfile1D(stream, dynamic_cast<const Profile1D&>(ao));
    } else if (aotype == "Profile2D") {
      writeProfile2D(stream, dynamic_cast<const Profile2D&>(ao));
    } else if (aotype == "Scatter1D") {
			writeScatter1D(stream, dynamic_cast<const Scatter1D&>(ao));
    } else if (aotype == "Scatter2D") {
      writeScatter2D(stream, dynamic_cast<const Scatter2D&>(ao));
    } else if (aotype == "Scatter3D") {
      writeScatter3D(stream, dynamic_cast<const Scatter3D&>(ao));
    } else if (aotype[0] == '_') {
      // Skip writing AO types with underscore prefixes (needed e.g. for Rivet wrappers)
      // maybe write a comment line in the output
    } else {
      ostringstream oss;
      oss << "Unrecognised analysis object type " << aotype << " in Writer::write";
      throw Exception(oss.str());
    }
  }


}
