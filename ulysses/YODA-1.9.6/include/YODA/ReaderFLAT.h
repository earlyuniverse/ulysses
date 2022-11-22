// -*- C++ -*-
//
// This file is part of YODA -- Yet more Objects for Data Analysis
// Copyright (C) 2008-2021 The YODA collaboration (see AUTHORS for details)
//
#ifndef YODA_READERFLAT_H
#define YODA_READERFLAT_H

#include "YODA/AnalysisObject.h"
#include "YODA/Reader.h"
#include "YODA/Index.h"
#include <YODA/Counter.h>
#include <YODA/Scatter1D.h>
#include <YODA/Scatter2D.h>
#include <YODA/Scatter3D.h>

namespace YODA {


  /// Persistency reader from YODA flat text data format.
  class ReaderFLAT : public Reader {
  public:

    /// Singleton creation function
    static Reader& create();
    
    void read(std::istream& stream, std::vector<AnalysisObject*>& aos);

    /// @brief Not implemented.
    Index mkIndex(std::istream& stream) {
      // Not supported
      return Index();
    }

  private:

    /// Private constructor, since it's a singleton.
    ReaderFLAT() { }

  };

}

#endif
