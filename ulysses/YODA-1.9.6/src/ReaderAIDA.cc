// -*- C++ -*-
//
// This file is part of YODA -- Yet more Objects for Data Analysis
// Copyright (C) 2008-2021 The YODA collaboration (see AUTHORS for details)
//
#include "YODA/ReaderAIDA.h"
#include "YODA/Utils/StringUtils.h"
#include "YODA/Exceptions.h"
#include "tinyxml/tinyxml.h"

#include "YODA/Counter.h"
#include "YODA/Histo1D.h"
#include "YODA/Histo2D.h"
#include "YODA/Profile1D.h"
#include "YODA/Profile2D.h"
// #include "YODA/Scatter1D.h"
#include "YODA/Scatter2D.h"
// #include "YODA/Scatter3D.h"

#include <iostream>
#include <locale>
using namespace std;

namespace YODA {

  /// Singleton creation function
  Reader& ReaderAIDA::create() {
    static ReaderAIDA _instance;
    return _instance;
  }

  void ReaderAIDA::_readDoc(std::istream& stream, vector<AnalysisObject*>& aos) {
    TiXmlDocument doc;
    stream >> doc;
    if (doc.Error()) {
      string err = "Error in " + string(doc.Value());
      err += ": " + string(doc.ErrorDesc());
      cerr << err << endl;
      throw ReadError(err);
    }

    // Return value, to be populated
    vector<AnalysisObject*> rtn;
    try {
      // Walk down tree to get to the <dataPointSet> elements
      const TiXmlNode* aidaN = doc.FirstChild("aida");
      if (!aidaN) throw ReadError("Couldn't get <aida> root element");
      for (const TiXmlNode* dpsN = aidaN->FirstChild("dataPointSet"); dpsN; dpsN = dpsN->NextSibling("dataPointSet")) {
        const TiXmlElement* dpsE = dpsN->ToElement();
        if (dpsE == 0) continue;
        const string plotpath = dpsE->Attribute("path");
        const string plotname = dpsE->Attribute("name");

        // DPS to be stored
        string sep = "/";
        if (plotpath.rfind("/") == plotpath.size()-1 || plotname.find("/") == 0) sep = "";
        /// @todo Clarify the memory management resulting from this... need shared_ptr?
        Scatter2D* dps = new Scatter2D(plotpath + sep + plotname);

        /// @todo This code crashes when there are annotations in the AIDA file: fix
        //// Read in annotations
        //for (const TiXmlNode* annN = dpsN->FirstChild("annotation"); annN; annN = annN->NextSibling()) {
        //  for (const TiXmlNode* itN = annN->FirstChild("item"); itN; itN = itN->NextSibling()) {
        //    dps->setAnnotation(itN->ToElement()->Attribute("key"), itN->ToElement()->Attribute("value"));
        //  }
        //}

        size_t ipt = 0;
        for (const TiXmlNode* dpN = dpsN->FirstChild("dataPoint"); dpN; dpN = dpN->NextSibling("dataPoint")) {
          ipt += 1;
          const TiXmlNode* xMeasN = dpN->FirstChild("measurement");
          if (!xMeasN) {
            cerr << "Couldn't get any <measurement> tag in DPS " << dpsE->Attribute("name") << " point #" << ipt << endl; continue;
          }
          const TiXmlNode* yMeasN = xMeasN->NextSibling("measurement");
          if (!yMeasN) {
            cerr << "Couldn't get y-axis <measurement> tag in DPS " << dpsE->Attribute("name") << " point #" << ipt << endl; continue;
          }
          const TiXmlElement* xMeasE = xMeasN->ToElement();
          const TiXmlElement* yMeasE = yMeasN->ToElement();
          const string xcentreStr   = xMeasE->Attribute("value");
          const string xerrplusStr  = xMeasE->Attribute("errorPlus");
          const string xerrminusStr = xMeasE->Attribute("errorMinus");
          const string ycentreStr   = yMeasE->Attribute("value");
          const string yerrplusStr  = yMeasE->Attribute("errorPlus");
          const string yerrminusStr = yMeasE->Attribute("errorMinus");
          //if (!centreStr) throw ReadError("Couldn't get a valid bin centre");
          //if (!errplusStr) throw ReadError("Couldn't get a valid bin err+");
          //if (!errminusStr) throw ReadError("Couldn't get a valid bin err-");
          istringstream xssC(xcentreStr);
          xssC.imbue(std::locale::classic()); // Interpret numbers in the "C" locale
          istringstream xssP(xerrplusStr);
          xssP.imbue(std::locale::classic()); // Interpret numbers in the "C" locale
          istringstream xssM(xerrminusStr);
          xssM.imbue(std::locale::classic()); // Interpret numbers in the "C" locale
          istringstream yssC(ycentreStr);
          yssC.imbue(std::locale::classic()); // Interpret numbers in the "C" locale
          istringstream yssP(yerrplusStr);
          yssP.imbue(std::locale::classic()); // Interpret numbers in the "C" locale
          istringstream yssM(yerrminusStr);
          yssM.imbue(std::locale::classic()); // Interpret numbers in the "C" locale
          double xcentre, xerrplus, xerrminus, ycentre, yerrplus, yerrminus;
          xssC >> xcentre; xssP >> xerrplus; xssM >> xerrminus;
          yssC >> ycentre; yssP >> yerrplus; yssM >> yerrminus;
          dps->addPoint(xcentre, ycentre, xerrminus, xerrplus, yerrminus, yerrplus);
        }
        aos.push_back(dps);

      }
    } catch (std::exception& e) {
      cerr << e.what() << endl;
      throw;
    }
  }



  // void ReaderAIDA::readGenericAO(std::ostream& os, const Histo1D& h) {
  // }


  // void ReaderAIDA::readProfile1D(std::ostream& os, const Profile1D& p) {
  // }


  // void ReaderAIDA::readScatter2D(std::ostream& os, const Scatter2D& s) {
  // }


}
