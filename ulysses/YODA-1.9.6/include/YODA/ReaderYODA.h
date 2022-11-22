// -*- C++ -*-
//
// This file is part of YODA -- Yet more Objects for Data Analysis
// Copyright (C) 2008-2021 The YODA collaboration (see AUTHORS for details)
//
#ifndef YODA_READERYODA_H
#define YODA_READERYODA_H

#include "YODA/AnalysisObject.h"
#include "YODA/Reader.h"
#include "YODA/Index.h"

namespace YODA {


  /// Persistency reader from YODA text data format.
  class ReaderYODA : public Reader {
  public:

    /// Singleton creation function
    static Reader& create();

    void read(std::istream& stream, std::vector<AnalysisObject*>& aos);
    
    Index mkIndex(std::istream& stream);

    // Include definitions of all read methods (all fulfilled by Reader::read(...))
    #include "YODA/ReaderMethods.icc"

  private:

    /// Private constructor, since it's a singleton.
    ReaderYODA() { }

  };


}

#endif
