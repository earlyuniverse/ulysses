// -*- C++ -*-
//
// This file is part of YODA -- Yet more Objects for Data Analysis
// Copyright (C) 2008-2021 The YODA collaboration (see AUTHORS for details)
//
#include "YODA/ReaderFLAT.h"
#include "YODA/Utils/StringUtils.h"
#include "YODA/Utils/getline.h"
#include "YODA/Exceptions.h"

#include "YODA/Counter.h"
#include "YODA/Scatter1D.h"
#include "YODA/Scatter2D.h"
#include "YODA/Scatter3D.h"

#include <iostream>
#include <locale>
using namespace std;

namespace YODA {

  /// Singleton creation function
  Reader& ReaderFLAT::create() {
    static ReaderFLAT _instance;
    return _instance;
  }


  void ReaderFLAT::read(istream& stream, vector<AnalysisObject*>& aos) {

    // Data format parsing states, representing current data type
    enum Context { NONE, //< outside any data block
                   SCATTER1D, SCATTER2D, SCATTER3D };

    /// State of the parser: line number, line, parser context, and pointer(s) to the object currently being assembled
    //unsigned int nline = 0;
    string s;
    Context context = NONE;
    //
    AnalysisObject* aocurr = NULL; //< Generic current AO pointer (useful or not?)
    Scatter1D* s1curr = NULL;
    Scatter2D* s2curr = NULL;
    Scatter3D* s3curr = NULL;

    // Loop over all lines of the input file
    while (Utils::getline(stream, s)) {
      //nline += 1;

      // Trim the line
      Utils::itrim(s);

      // Ignore blank lines
      if (s.empty()) continue;

      // Ignore comments (whole-line only, without indent, and still allowed for compatibility on BEGIN/END lines)
      if (s.find("#") == 0 && s.find("BEGIN") == string::npos && s.find("END") == string::npos) continue;


      // Now the context-sensitive part
      if (context == NONE) {

        // We require a BEGIN line to start a context
        if (s.find("BEGIN ") == string::npos) throw ReadError("Unexpected line in FLAT format parsing when BEGIN expected");

        // Split into parts
        vector<string> parts;
        istringstream iss(s); string tmp;

        iss.imbue(std::locale::classic()); // Interpret numbers in the "C" locale

        while (iss >> tmp) {
          if (tmp != "#") parts.push_back(tmp);
        }

        // Extract context from BEGIN type
        assert(parts.size() >= 2 && parts[0] == "BEGIN");
        const string ctxstr = parts[1];

        // Get block path if possible
        const string path = (parts.size() >= 3) ? parts[2] : "";

        // Set the new context and create a new AO to populate
        if (ctxstr == "VALUE") {
          context = SCATTER1D;
          s1curr = new Scatter1D(path);
          aocurr = s1curr;
        } else if (ctxstr == "HISTO1D" || ctxstr == "HISTOGRAM") {
          context = SCATTER2D;
          s2curr = new Scatter2D(path);
          aocurr = s2curr;
        } else if (ctxstr == "HISTO2D" || ctxstr == "HISTOGRAM2D") {
          context = SCATTER3D;
          s3curr = new Scatter3D(path);
          aocurr = s3curr;
        }
        // cout << aocurr->path() << " " << nline << " " << context << endl;

      } else {
        /// @todo Flatten conditional blocks with more else-ifs?

        // Throw error if a BEGIN line is found
        if (s.find("BEGIN ") != string::npos) throw ReadError("Unexpected BEGIN line in FLAT format parsing before ending current BEGIN..END block");

        // Clear/reset context and register AO if END line is found
        /// @todo Throw error if mismatch between BEGIN (context) and END types
        if (s.find("END ") != string::npos) {
          aos.push_back(aocurr);
          context = NONE;
          aocurr = NULL; s1curr = NULL; s2curr = NULL; s3curr = NULL;
          continue; ///< @todo Improve... would be good to avoid these continues
        }

        // Extract annotations for all types
        const size_t ieq = s.find("=");
        if (ieq != string::npos) {
          const string akey = s.substr(0, ieq);
          const string aval = s.substr(ieq+1);
          aocurr->setAnnotation(akey, aval);
          continue; ///< @todo Improve... would be good to avoid these continues
        }

        // Populate the data lines for points
        istringstream iss(s);
        switch (context) {

        case SCATTER1D:
          {
            double x(0), exm(0), exp(0);
            iss >> x >> exm >> exp;
            s1curr->addPoint(Point1D(x, exm, exp));
          }
          break;

        case SCATTER2D:
          {
            double xlow(0), xhigh(0), y(0), eym(0), eyp(0);
            iss >> xlow >> xhigh >> y >> eym >> eyp;
            const double x = (xlow + xhigh)/2.0;
            const double ex = (xhigh - xlow)/2.0;
            s2curr->addPoint(Point2D(x, y, ex, ex, eym, eyp));
          }
          break;

        case SCATTER3D:
          {
            double xlow(0), xhigh(0), ylow(0), yhigh(0), z(0), ezm(0), ezp(0);
            iss >> xlow >> xhigh >> ylow >> yhigh >> z >> ezm >> ezp;
            const double x = (xlow + xhigh)/2.0;
            const double ex = (xhigh - xlow)/2.0;
            const double y = (ylow + yhigh)/2.0;
            const double ey = (yhigh - ylow)/2.0;
            s3curr->addPoint(Point3D(x, y, z, ex, ex, ey, ey, ezm, ezp));
          }
          break;

        default:
          throw ReadError("Unknown context in FLAT format parsing: how did this happen?");
        }

      }
    }

  }


}
