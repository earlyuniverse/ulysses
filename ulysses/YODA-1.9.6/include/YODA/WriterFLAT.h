// -*- C++ -*-
//
// This file is part of YODA -- Yet more Objects for Data Analysis
// Copyright (C) 2008-2021 The YODA collaboration (see AUTHORS for details)
//
#ifndef YODA_WRITERFLAT_H
#define YODA_WRITERFLAT_H

#include "YODA/AnalysisObject.h"
#include "YODA/Writer.h"

namespace YODA {


  /// Persistency writer for flat text format.
  class WriterFLAT : public Writer {
  public:

    /// Singleton creation function
    static Writer& create();

    // Include definitions of all write methods (all fulfilled by Writer::write(...))
    #include "YODA/WriterMethods.icc"


  protected:

    void writeCounter(std::ostream& stream, const Counter& c);
    void writeHisto1D(std::ostream& stream, const Histo1D& h);
    void writeHisto2D(std::ostream& stream, const Histo2D& h);
    void writeProfile1D(std::ostream& stream, const Profile1D& p);
    void writeProfile2D(std::ostream& stream, const Profile2D& p);
    void writeScatter1D(std::ostream& stream, const Scatter1D& s);
    void writeScatter2D(std::ostream& stream, const Scatter2D& s);
    void writeScatter3D(std::ostream& stream, const Scatter3D& s);


  private:

    void _writeAnnotations(std::ostream& os, const AnalysisObject& ao);

    /// Private since it's a singleton.
    WriterFLAT() { }

  };


}

#endif
