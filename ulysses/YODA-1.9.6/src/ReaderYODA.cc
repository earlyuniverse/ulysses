// -*- C++ -*-
//
// This file is part of YODA -- Yet more Objects for Data Analysis
// Copyright (C) 2008-2021 The YODA collaboration (see AUTHORS for details)
//
#include "YODA/ReaderYODA.h"
#include "YODA/Reader.h"
#include "YODA/Index.h"
#include "YODA/Utils/StringUtils.h"
#include "YODA/Utils/getline.h"
#include "YODA/Exceptions.h"
#include "YODA/Config/DummyConfig.h"

#include "YODA/Counter.h"
#include "YODA/Histo1D.h"
#include "YODA/Histo2D.h"
#include "YODA/Profile1D.h"
#include "YODA/Profile2D.h"
#include "YODA/Scatter1D.h"
#include "YODA/Scatter2D.h"
#include "YODA/Scatter3D.h"

#include "yaml-cpp/yaml.h"
#ifdef YAML_NAMESPACE
#define YAML YAML_NAMESPACE
#endif

#ifdef HAVE_LIBZ
#define _XOPEN_SOURCE 700
#include "zstr/zstr.hpp"
#endif

#include <iostream>
#include <locale>
#include <string>
#include <cstring>
using namespace std;

namespace YODA {

  /// Singleton creation function
  Reader& ReaderYODA::create() {
    static ReaderYODA _instance;
    return _instance;
  }

  namespace {

    /// Fast ASCII tokenizer, extended from FastIStringStream by Gavin Salam.
    class aistringstream {
    public:
      // Constructor from char*
      aistringstream(const char* line=0) {
        reset(line);
        _set_locale();
      }
      // Constructor from std::string
      aistringstream(const string& line) {
        reset(line);
        _set_locale();
      }
      ~aistringstream() {
        _reset_locale();
      }

      // Re-init to new line as char*
      void reset(const char* line=0) {
        _next = const_cast<char*>(line);
        _new_next = _next;
        _error = false;
      }
      // Re-init to new line as std::string
      void reset(const string& line) { reset(line.c_str()); }

      // Tokenizing stream operator (forwards to specialisations)
      template<class T>
      aistringstream& operator >> (T& value) {
        _get(value);
        if (_new_next == _next) _error = true; // handy error condition behaviour!
        _next = _new_next;
        return *this;
      }

      // Allow use of operator>> in a while loop
      operator bool() const { return !_error; }


    private:

      // Changes the thread-local locale to interpret numbers in the "C" locale
      void _set_locale() {
        _locale_set = newlocale(LC_NUMERIC_MASK, "C", NULL);
        _locale_prev = uselocale(_locale_set);
        if (!_locale_prev) {
          throw ReadError(std::string("Error setting locale: ") + strerror(errno));
        }
      }
      void _reset_locale() {
        if (!uselocale(_locale_prev)) {
          throw ReadError(std::string("Error setting locale: ") + strerror(errno));
        }
        freelocale(_locale_set);
      }

      void _get(double& x) { x = std::strtod(_next, &_new_next); }
      void _get(float& x) { x = std::strtof(_next, &_new_next); }
      void _get(int& i) { i = std::strtol(_next, &_new_next, 10); } // force base 10!
      void _get(long& i) { i = std::strtol(_next, &_new_next, 10); } // force base 10!
      void _get(unsigned int& i) { i = std::strtoul(_next, &_new_next, 10); } // force base 10!
      void _get(long unsigned int& i) { i = std::strtoul(_next, &_new_next, 10); } // force base 10!
      void _get(string& x) {
        /// @todo If _next and _new_next become null?
        while (std::isspace(*_next)) _next += 1;
        _new_next = _next;
        while (!std::isspace(*_new_next)) _new_next += 1;
        x = string(_next, _new_next-_next);
      }

      locale_t _locale_set, _locale_prev;
      char *_next, *_new_next;
      bool _error;
    };

  }


  void ReaderYODA::read(istream& inputStream, vector<AnalysisObject*>& aos) {

    #ifdef HAVE_LIBZ
    // NB. zstr auto-detects if file is deflated or plain-text
    zstr::istream stream(inputStream);
    #else
    istream& stream = inputStream;
    #endif
    // Data format parsing states, representing current data type
    /// @todo Extension to e.g. "bar" or multi-counter or binned-value types, and new formats for extended Scatter types
    enum Context { NONE, //< outside any data block
                   SCATTER1D, SCATTER2D, SCATTER3D,
                   COUNTER,
                   HISTO1D, HISTO2D,
                   PROFILE1D, PROFILE2D };

    /// State of the parser: line number, line, parser context, and pointer(s) to the object currently being assembled
    unsigned int nline = 0;
    string s;
    Context context = NONE;
    //
    AnalysisObject* aocurr = NULL; //< Generic current AO pointer
    vector<HistoBin1D> h1binscurr; //< Current H1 bins container
    vector<HistoBin2D> h2binscurr; //< Current H2 bins container
    vector<ProfileBin1D> p1binscurr; //< Current P1 bins container
    vector<ProfileBin2D> p2binscurr; //< Current P2 bins container
    vector<Point1D> pt1scurr; //< Current Point1Ds container
    vector<Point2D> pt2scurr; //< Current Point2Ds container
    vector<Point3D> pt3scurr; //< Current Point3Ds container
    Counter* cncurr = NULL;
    Histo1D* h1curr = NULL;
    Histo2D* h2curr = NULL;
    Profile1D* p1curr = NULL;
    Profile2D* p2curr = NULL;
    Scatter1D* s1curr = NULL;
    Scatter2D* s2curr = NULL;
    Scatter3D* s3curr = NULL;
    //std::vector<std::string> variationscurr;
    string annscurr;

    // Loop over all lines of the input file
    aistringstream aiss;
    bool in_anns = false;
    string fmt = "1";
    //int nfmt = 1;
    while (Utils::getline(stream, s)) {
      nline += 1;


      // CLEAN LINES IF NOT IN ANNOTATION MODE
      if (!in_anns) {
        // Trim the line
        Utils::itrim(s);

        // Ignore blank lines
        if (s.empty()) continue;

        // Ignore comments (whole-line only, without indent, and still allowed for compatibility on BEGIN/END lines)
        if (s.find("#") == 0 && s.find("BEGIN") == string::npos && s.find("END") == string::npos) continue;
      }


      // STARTING A NEW CONTEXT
      if (context == NONE) {

        // We require a BEGIN line to start a context
        if (s.find("BEGIN ") == string::npos) {
          stringstream ss;
          ss << "Unexpected line in YODA format parsing when BEGIN expected: '" << s << "' on line " << nline;
          throw ReadError(ss.str());
        }

        // Remove leading #s from the BEGIN line if necessary
        while (s.find("#") == 0) s = Utils::trim(s.substr(1));

        // Split into parts
        vector<string> parts;
        istringstream iss(s); string tmp;
        while (iss >> tmp) parts.push_back(tmp);

        // Extract context from BEGIN type
        if (parts.size() < 2 || parts[0] != "BEGIN") {
          stringstream ss;
          ss << "Unexpected BEGIN line structure when BEGIN expected: '" << s << "' on line " << nline;
          throw ReadError(ss.str());
        }

        // Second part is the context name
        const string ctxstr = parts[1];

        // Get block path if possible
        const string path = (parts.size() >= 3) ? parts[2] : "";

        // Set the new context and create a new AO to populate
        /// @todo Use the block format version for (occasional, careful) format evolution
        if (Utils::startswith(ctxstr, "YODA_COUNTER")) {
          context = COUNTER;
          cncurr = new Counter(path);
          aocurr = cncurr;
        } else if (Utils::startswith(ctxstr, "YODA_SCATTER1D")) {
          context = SCATTER1D;
          s1curr = new Scatter1D(path);
          aocurr = s1curr;
        } else if (Utils::startswith(ctxstr, "YODA_SCATTER2D")) {
          context = SCATTER2D;
          s2curr = new Scatter2D(path);
          aocurr = s2curr;
        } else if (Utils::startswith(ctxstr, "YODA_SCATTER3D")) {
          context = SCATTER3D;
          s3curr = new Scatter3D(path);
          aocurr = s3curr;
        } else if (Utils::startswith(ctxstr, "YODA_HISTO1D")) {
          context = HISTO1D;
          h1curr = new Histo1D(path);
          aocurr = h1curr;
        } else if (Utils::startswith(ctxstr, "YODA_HISTO2D")) {
          context = HISTO2D;
          h2curr = new Histo2D(path);
          aocurr = h2curr;
        } else if (Utils::startswith(ctxstr, "YODA_PROFILE1D")) {
          context = PROFILE1D;
          p1curr = new Profile1D(path);
          aocurr = p1curr;
        } else if (Utils::startswith(ctxstr, "YODA_PROFILE2D")) {
          context = PROFILE2D;
          p2curr = new Profile2D(path);
          aocurr = p2curr;
        }
        // cout << aocurr->path() << " " << nline << " " << context << endl;

        // Get block format version if possible (assume version=1 if none found)
        const size_t vpos = ctxstr.find_last_of("V");
        fmt = vpos != string::npos ? ctxstr.substr(vpos+1) : "1";
        // cout << fmt << endl;

        // From version 2 onwards, use the in_anns state from BEGIN until ---
        if (fmt != "1") in_anns = true;


      } else { //< not a BEGIN line


        // Throw error if a BEGIN line is found
        if (s.find("BEGIN ") != string::npos) ///< @todo require pos = 0 from fmt=V2
          throw ReadError("Unexpected BEGIN line in YODA format parsing before ending current BEGIN..END block");


        // FINISHING THE CURRENT CONTEXT
        // Clear/reset context and register AO
        /// @todo Throw error if mismatch between BEGIN (context) and END types
        if (s.find("END ") != string::npos) { ///< @todo require pos = 0 from fmt=V2
          switch (context) {
          case COUNTER:
            break;
          case HISTO1D:
            h1curr->addBins(h1binscurr);
            h1binscurr.clear();
            break;
          case HISTO2D:
            h2curr->addBins(h2binscurr);
            h2binscurr.clear();
            break;
          case PROFILE1D:
            p1curr->addBins(p1binscurr);
            p1binscurr.clear();
            break;
          case PROFILE2D:
            p2curr->addBins(p2binscurr);
            p2binscurr.clear();
            break;
          case SCATTER1D:
            for (auto& p : pt1scurr) p.setParent(s1curr);
            s1curr->addPoints(pt1scurr);
            pt1scurr.clear();
            break;
          case SCATTER2D:
            for (auto& p : pt2scurr) p.setParent(s2curr);
            s2curr->addPoints(pt2scurr);
            pt2scurr.clear();
            break;
          case SCATTER3D:
            for (auto& p : pt3scurr) p.setParent(s3curr);
            s3curr->addPoints(pt3scurr);
            pt3scurr.clear();
            break;
          case NONE:
            break;
          }

          // Set all annotations
          try {
            YAML::Node anns = YAML::Load(annscurr);
            // for (YAML::const_iterator it = anns.begin(); it != anns.end(); ++it) {
            for (const auto& it : anns) {
              const string key = it.first.as<string>();
              // const string val = it.second.as<string>();
              YAML::Emitter em;
              em << YAML::Flow << it.second; //< use single-line formatting, for lists & maps
              const string val = em.c_str();
             // if (!(key.find("ErrorBreakdown") != string::npos))
              aocurr->setAnnotation(key, val);
            }
          } catch (...) {
            /// @todo Is there a case for just giving up on these annotations, printing the error msg, and keep going? As an option?
            const string err = "Problem during annotation parsing of YAML block:\n'''\n" + annscurr + "\n'''";
            // cerr << err << endl;
            throw ReadError(err);
          }
          annscurr.clear();
          //variationscurr.clear();
          in_anns = false;

          // Put this AO in the completed stack
          aos.push_back(aocurr);

          // Clear all current-object pointers
          aocurr = nullptr;
          cncurr = nullptr;
          h1curr = nullptr; h2curr = nullptr;
          p1curr = nullptr; p2curr = nullptr;
          s1curr = nullptr; s2curr = nullptr; s3curr = nullptr;
          context = NONE;
          continue;
        }


        // ANNOTATIONS PARSING
        if (fmt == "1") {
          // First convert to one-key-per-line YAML syntax
          const size_t ieq = s.find("=");
          if (ieq != string::npos) s.replace(ieq, 1, ": ");
          // Special-case treatment for syntax clashes
          const size_t icost = s.find(": *");
          if (icost != string::npos) {
            s.replace(icost, 1, ": '*");
            s += "'";
          }
          // Store reformatted annotation
          const size_t ico = s.find(":");
          if (ico != string::npos) {
            annscurr += (annscurr.empty() ? "" : "\n") + s;
            continue;
          }
        } else if (in_anns) {
          if (s == "---") {
            in_anns = false;
          } else {
            annscurr += (annscurr.empty() ? "" : "\n") + s;
            // In order to handle multi-error points in scatters, we need to know which variations are stored, if any
            // can't wait until we process the annotations at the end, since need to know when filling points.
            // This is a little inelegant though...
            //if (s.find("ErrorBreakdown") != string::npos) {
             // errorBreakdown = YAML::Load(s)["ErrorBreakdown"];
              //for (const auto& it : errorBreakdown) {
              //  const string val0 = it.first.as<string>();
                //for (const auto& it2 : it.second) {
                //  const string val = it2.as<string>();
                //}
             // }
           // }
          }
          continue;
        }


        // DATA PARSING
        aiss.reset(s);
        // double sumw(0), sumw2(0), sumwx(0), sumwx2(0), sumwy(0), sumwy2(0), sumwz(0), sumwz2(0), sumwxy(0), sumwxz(0), sumwyz(0), n(0);
        switch (context) {

        case COUNTER:
          {
            double sumw(0), sumw2(0), n(0);
            aiss >> sumw >> sumw2 >> n;
            cncurr->setDbn(Dbn0D(n, sumw, sumw2));
          }
          break;

        case HISTO1D:
          {
            string xoflow1, xoflow2; double xmin(0), xmax(0);
            double sumw(0), sumw2(0), sumwx(0), sumwx2(0), n(0);
            /// @todo Improve/factor this "bin" string-or-float parsing... esp for mixed case of 2D overflows
            /// @todo When outflows are treated as "infinity bins" and don't require a distinct type, string replace under/over -> -+inf
            if (s.find("Total") != string::npos || s.find("Underflow") != string::npos || s.find("Overflow") != string::npos) {
              aiss >> xoflow1 >> xoflow2;
            } else {
              aiss >> xmin >> xmax;
            }
            // The rest is the same for overflows and in-range bins
            aiss >> sumw >> sumw2 >> sumwx >> sumwx2 >> n;
            const Dbn1D dbn(n, sumw, sumw2, sumwx, sumwx2);
            if (xoflow1 == "Total") h1curr->setTotalDbn(dbn);
            else if (xoflow1 == "Underflow") h1curr->setUnderflow(dbn);
            else if (xoflow1 == "Overflow")  h1curr->setOverflow(dbn);
            // else h1curr->addBin(HistoBin1D(std::make_pair(xmin,xmax), dbn));
            else h1binscurr.push_back(HistoBin1D(std::make_pair(xmin,xmax), dbn));
          }
          break;

        case HISTO2D:
          {
            string xoflow1, xoflow2, yoflow1, yoflow2; double xmin(0), xmax(0), ymin(0), ymax(0);
            double sumw(0), sumw2(0), sumwx(0), sumwx2(0), sumwy(0), sumwy2(0), sumwxy(0), n(0);
            /// @todo Improve/factor this "bin" string-or-float parsing... esp for mixed case of 2D overflows
            /// @todo When outflows are treated as "infinity bins" and don't require a distinct type, string replace under/over -> -+inf
            if (s.find("Total") != string::npos) {
              aiss >> xoflow1 >> xoflow2; // >> yoflow1 >> yoflow2;
            } else if (s.find("Underflow") != string::npos || s.find("Overflow") != string::npos) {
              throw ReadError("2D histogram overflow syntax is not yet defined / handled");
            } else {
              aiss >> xmin >> xmax >> ymin >> ymax;
            }
            // The rest is the same for overflows and in-range bins
            aiss >> sumw >> sumw2 >> sumwx >> sumwx2 >> sumwy >> sumwy2 >> sumwxy >> n;
            const Dbn2D dbn(n, sumw, sumw2, sumwx, sumwx2, sumwy, sumwy2, sumwxy);
            if (xoflow1 == "Total") h2curr->setTotalDbn(dbn);
            // else if (xoflow1 == "Underflow") p1curr->setUnderflow(dbn);
            // else if (xoflow1 == "Overflow")  p1curr->setOverflow(dbn);
            else {
              assert(xoflow1.empty());
              // h2curr->addBin(HistoBin2D(std::make_pair(xmin,xmax), std::make_pair(ymin,ymax), dbn));
              h2binscurr.push_back(HistoBin2D(std::make_pair(xmin,xmax), std::make_pair(ymin,ymax), dbn));
            }
          }
          break;

        case PROFILE1D:
          {
            string xoflow1, xoflow2; double xmin(0), xmax(0);
            double sumw(0), sumw2(0), sumwx(0), sumwx2(0), sumwy(0), sumwy2(0), n(0);
            /// @todo Improve/factor this "bin" string-or-float parsing... esp for mixed case of 2D overflows
            /// @todo When outflows are treated as "infinity bins" and don't require a distinct type, string replace under/over -> -+inf
            if (s.find("Total") != string::npos || s.find("Underflow") != string::npos || s.find("Overflow") != string::npos) {
              aiss >> xoflow1 >> xoflow2;
            } else {
              aiss >> xmin >> xmax;
            }
            // The rest is the same for overflows and in-range bins
            aiss >> sumw >> sumw2 >> sumwx >> sumwx2 >> sumwy >> sumwy2 >> n;
            const double DUMMYWXY = 0;
            const Dbn2D dbn(n, sumw, sumw2, sumwx, sumwx2, sumwy, sumwy2, DUMMYWXY);
            if (xoflow1 == "Total") p1curr->setTotalDbn(dbn);
            else if (xoflow1 == "Underflow") p1curr->setUnderflow(dbn);
            else if (xoflow1 == "Overflow")  p1curr->setOverflow(dbn);
            // else p1curr->addBin(ProfileBin1D(std::make_pair(xmin,xmax), dbn));
            else p1binscurr.push_back(ProfileBin1D(std::make_pair(xmin,xmax), dbn));
          }
          break;

        case PROFILE2D:
          {
            string xoflow1, xoflow2, yoflow1, yoflow2; double xmin(0), xmax(0), ymin(0), ymax(0);
            double sumw(0), sumw2(0), sumwx(0), sumwx2(0), sumwy(0), sumwy2(0), sumwz(0), sumwz2(0), sumwxy(0), sumwxz(0), sumwyz(0), n(0);
            /// @todo Improve/factor this "bin" string-or-float parsing... esp for mixed case of 2D overflows
            /// @todo When outflows are treated as "infinity bins" and don't require a distinct type, string replace under/over -> -+inf
            if (s.find("Total") != string::npos) {
              aiss >> xoflow1 >> xoflow2; // >> yoflow1 >> yoflow2;
            } else if (s.find("Underflow") != string::npos || s.find("Overflow") != string::npos) {
              throw ReadError("2D profile overflow syntax is not yet defined / handled");
            } else {
              aiss >> xmin >> xmax >> ymin >> ymax;
            }
            // The rest is the same for overflows and in-range bins
            aiss >> sumw >> sumw2 >> sumwx >> sumwx2 >> sumwy >> sumwy2 >> sumwz >> sumwz2 >> sumwxy >> sumwxz >> sumwyz >> n;
            const Dbn3D dbn(n, sumw, sumw2, sumwx, sumwx2, sumwy, sumwy2, sumwz, sumwz2, sumwxy, sumwxz, sumwyz);
            if (xoflow1 == "Total") p2curr->setTotalDbn(dbn);
            // else if (xoflow1 == "Underflow") p2curr->setUnderflow(dbn);
            // else if (xoflow1 == "Overflow")  p2curr->setOverflow(dbn);
            else {
              assert(xoflow1.empty());
              // p2curr->addBin(ProfileBin2D(std::make_pair(xmin,xmax), std::make_pair(ymin,ymax), dbn));
              p2binscurr.push_back(ProfileBin2D(std::make_pair(xmin,xmax), std::make_pair(ymin,ymax), dbn));
            }
          }
          break;

        case SCATTER1D:
          {
            double x(0), exm(0), exp(0);
            aiss >> x >> exm >> exp;
            // set nominal point
            Point1D thispoint=Point1D(x, exm, exp);
            // check if we stored variations of this point
            //if (variationscurr.size()>0){
            //  // for each variation, store the alt errors.
            //  // start at 1 since we have already filled nominal !
            //  for (unsigned int ivar=1; ivar<variationscurr.size(); ivar++){
            //   std::string thisvariation=variationscurr[ivar];
            //    aiss >> exm >> exp;
            //    thispoint.setXErrs(exm,exp,thisvariation);
            //  }
            //}
            pt1scurr.push_back(thispoint);
          }
          break;

        case SCATTER2D:
          {
            double x(0), y(0), exm(0), exp(0), eym(0), eyp(0);
            aiss >> x >> exm >> exp >> y >> eym >> eyp;
            // set nominal point
            Point2D thispoint=Point2D(x, y, exm, exp, eym, eyp);
            // check if we stored variations of this point
            // for each variation, store the alt errors.
            // start at 1 since we have already filled nominal !
            pt2scurr.push_back(thispoint);
          }
          break;

        case SCATTER3D:
          {
            double x(0), y(0), z(0), exm(0), exp(0), eym(0), eyp(0), ezm(0), ezp(0);
            aiss >> x >> exm >> exp >> y >> eym >> eyp >> z >> ezm >> ezp;
            // set nominal point
            Point3D thispoint=Point3D(x, y, z, exm, exp, eym, eyp, ezm, ezp);
            pt3scurr.push_back(thispoint);
          }
          break;

        default:
          throw ReadError("Unknown context in YODA format parsing: how did this happen?");

          }
        // cout << "AO CONTENT " << nline << endl;
        // cout << "  " << xmin << " " << xmax << " " << ymin << " " << ymax << " / '" << xoflow1 << "' '" << xoflow2 << "' '" << yoflow1 << "' '" << yoflow2 << "'" << endl;
        // cout << "  " << sumw << " " << sumw2 << " " << sumwx << " " << sumwx2 << " " << sumwy << " " << sumwy2 << " " << sumwz << " " << sumwz2 << " " << sumwxy << " " << sumwxz << " " << sumwyz << " " << n << endl;
        // cout << "  " << x << " " << y << " " << z << " " << exm << " " << exp << " " << eym << " " << eyp << " " << ezm << " " << ezp << endl;
        }
     }
  }


  Index ReaderYODA::mkIndex(std::istream& inputStream) {
    Index::AOIndex hmap;
    hmap.insert({"Histo1D", unordered_map<string, int>()});
    hmap.insert({"Histo2D", unordered_map<string, int>()});
    hmap.insert({"Profile1D", unordered_map<string, int>()});
    hmap.insert({"Profile2D", unordered_map<string, int>()});
    hmap.insert({"Scatter1D", unordered_map<string, int>()});
    hmap.insert({"Scatter2D", unordered_map<string, int>()});
    hmap.insert({"Scatter3D", unordered_map<string, int>()});
    hmap.insert({"Counter", unordered_map<string, int>()});

    //#ifdef HAVE_LIBZ
    //// NB. zstr auto-detects if file is deflated or plain-text
    // zstr::istream stream(inputStream);
    //#else
    istream& stream = inputStream;
    //#endif
    enum Context {
      NONE,
      SCATTER1D,
      SCATTER2D,
      SCATTER3D,
      COUNTER,
      HISTO1D,
      HISTO2D,
      PROFILE1D,
      PROFILE2D
    };

    /// State of the parser: line number, line, parser context, and pointer(s)
    /// to the object currently being assembled
    unsigned int nline = 0;
    string s;
    Context context = NONE;
    std::string curpath = "";
    int nbins = 0;

    // Loop over all lines of the input file
    bool in_anns = false;
    string fmt = "1";
    while (std::getline(stream, s)) {
      nline += 1;
      if (!in_anns) {
        Utils::itrim(s);
        if (s.empty())
          continue;
        if (s.find("#") == 0 && s.find("BEGIN") == string::npos &&
            s.find("END") == string::npos)
          continue;
      }
      if (context == NONE) {
        if (s.find("BEGIN ") == string::npos) {
          stringstream ss;
          ss << "Unexpected line in YODA format parsing when BEGIN expected: '"
             << s << "' on line " << nline;
          throw ReadError(ss.str());
        }
        while (s.find("#") == 0)
          s = Utils::trim(s.substr(1));
        vector<string> parts;
        istringstream iss(s);
        string tmp;
        while (iss >> tmp)
          parts.push_back(tmp);

        if (parts.size() < 2 || parts[0] != "BEGIN") {
          stringstream ss;
          ss << "Unexpected BEGIN line structure when BEGIN expected: '" << s
             << "' on line " << nline;
          throw ReadError(ss.str());
        }
        const string ctxstr = parts[1];
        curpath = (parts.size() >= 3) ? parts[2] : "";
        nbins = 0;
        if (Utils::startswith(ctxstr, "YODA_COUNTER")) {
          context = COUNTER;
        } else if (Utils::startswith(ctxstr, "YODA_SCATTER1D")) {
          context = SCATTER1D;
        } else if (Utils::startswith(ctxstr, "YODA_SCATTER2D")) {
          context = SCATTER2D;
        } else if (Utils::startswith(ctxstr, "YODA_SCATTER3D")) {
          context = SCATTER3D;
        } else if (Utils::startswith(ctxstr, "YODA_HISTO1D")) {
          context = HISTO1D;
        } else if (Utils::startswith(ctxstr, "YODA_HISTO2D")) {
          context = HISTO2D;
        } else if (Utils::startswith(ctxstr, "YODA_PROFILE1D")) {
          context = PROFILE1D;
        } else if (Utils::startswith(ctxstr, "YODA_PROFILE2D")) {
          context = PROFILE2D;
        }
        const size_t vpos = ctxstr.find_last_of("V");
        fmt = vpos != string::npos ? ctxstr.substr(vpos + 1) : "1";
        if (fmt != "1")
          in_anns = true;
      } else { //< not a BEGIN line
        if (s.find("BEGIN ") != string::npos)
          throw ReadError("Unexpected BEGIN line in YODA format parsing before "
                          "ending current BEGIN..END block");
        // FINISHING THE CURRENT CONTEXT
        if (s.find("END ") != string::npos) {
          if (context == HISTO1D) {
            hmap["Histo1D"].insert({curpath, nbins});
          } else if (context == HISTO2D) {
            hmap["Histo2D"].insert({curpath, nbins});
          } else if (context == PROFILE1D) {
            hmap["Profile1D"].insert({curpath, nbins});
          } else if (context == PROFILE2D) {
            hmap["Profile2D"].insert({curpath, nbins});
          } else if (context == SCATTER1D) {
            hmap["Scatter1D"].insert({curpath, nbins});
          } else if (context == SCATTER2D) {
            hmap["Scatter2D"].insert({curpath, nbins});
          } else if (context == SCATTER3D) {
            hmap["Scatter3D"].insert({curpath, nbins});
          } else if (context == COUNTER) {
            hmap["Counter"].insert({curpath, nbins});
          }
          in_anns = false;
          context = NONE;
          continue;
        }
        // ANNOTATIONS PARSING
        if (fmt == "1") {
          // First convert to one-key-per-line YAML syntax
          const size_t ieq = s.find("=");
          if (ieq != string::npos)
            s.replace(ieq, 1, ": ");
          // Special-case treatment for syntax clashes
          const size_t icost = s.find(": *");
          if (icost != string::npos) {
            s.replace(icost, 1, ": '*");
            s += "'";
          }
          // Store reformatted annotation
          const size_t ico = s.find(":");
          if (ico != string::npos)
            continue;
        } else if (in_anns) {
          if (s == "---")
            in_anns = false;
          continue;
        }

        if ((context == HISTO1D) || (context == HISTO2D) ||
            (context == PROFILE1D) || (context == PROFILE2D)) {
          if (s.find("Total") != string::npos ||
              s.find("Underflow") != string::npos ||
              s.find("Overflow") != string::npos)
            continue;
          nbins++;
        } else if ((context == SCATTER1D) || (context == SCATTER2D) ||
                   (context == SCATTER3D)) {
          nbins++;
        }
      }
    }
    return Index(hmap);
  }
}
