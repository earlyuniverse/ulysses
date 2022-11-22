// -*- C++ -*-
//
// This file is part of YODA -- Yet more Objects for Data Analysis
// Copyright (C) 2008-2021 The YODA collaboration (see AUTHORS for details)
//
#include "YODA/Reader.h"
#include "YODA/ReaderYODA.h"
#include "YODA/ReaderAIDA.h"
#include "YODA/ReaderFLAT.h"
#include "YODA/Config/DummyConfig.h"

using namespace std;

namespace YODA {


  Reader& mkReader(const string& name) {
    // Determine the format from the string (a file or file extension)
    const size_t lastdot = name.find_last_of(".");
    string fmt = Utils::toLower(lastdot == string::npos ? name : name.substr(lastdot+1));
    if (fmt == "gz") {
      #ifndef HAVE_LIBZ
      throw UserError("YODA was compiled without zlib support: can't read " + name);
      #endif
      const size_t lastbutonedot = (lastdot == string::npos) ? string::npos : name.find_last_of(".", lastdot-1);
      fmt = Utils::toLower(lastbutonedot == string::npos ? name : name.substr(lastbutonedot+1));
    }
    // Create the appropriate Reader
    if (Utils::startswith(fmt, "yoda")) return ReaderYODA::create();
    if (Utils::startswith(fmt, "aida")) return ReaderAIDA::create();
    if (Utils::startswith(fmt, "dat" )) return ReaderFLAT::create(); ///< @todo Improve/remove... .ydat?
    if (Utils::startswith(fmt, "flat")) return ReaderFLAT::create();
    throw UserError("Format cannot be identified from string '" + name + "'");
  }


}
