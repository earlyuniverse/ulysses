#include "YODA/Scatter1D.h"
#include "YODA/Counter.h"
#include <sstream>
#include "yaml-cpp/yaml.h"
#ifdef YAML_NAMESPACE
#define YAML YAML_NAMESPACE
#endif

namespace YODA {


  /// Make a Scatter1D representation of a Histo1D
  Scatter1D mkScatter(const Counter& c) {
    Scatter1D rtn;
    for (const std::string& a : c.annotations())
      rtn.setAnnotation(a, c.annotation(a));
    rtn.setAnnotation("Type", c.type()); // might override the copied ones
    Point1D pt(c.val(), c.err());
    pt.setParent(&rtn);
    rtn.addPoint(pt);
    return rtn;
  }


  void Scatter1D::updateTotalUncertainty(){
    for (size_t thisPointIndex = 0; thisPointIndex < this->numPoints(); ++thisPointIndex) {
      Point1D& thisPoint = this->_points[thisPointIndex];
      thisPoint.updateTotalUncertainty();
    }
  }


  void Scatter1D::parseVariations() {
    if (this->_variationsParsed) { return;}
    if (!(this->hasAnnotation("ErrorBreakdown"))) { return; }
    YAML::Node errorBreakdown;
    errorBreakdown = YAML::Load(this->annotation("ErrorBreakdown"));

    if (errorBreakdown.size()) {
      for (size_t thisPointIndex = 0; thisPointIndex < this->numPoints(); ++thisPointIndex) {
        Point1D &thispoint = this->_points[thisPointIndex];
        YAML::Node variations = errorBreakdown[thisPointIndex];
        for (const auto& variation : variations) {
          const std::string variationName = variation.first.as<std::string>();
          // The empty-name variation is the total and should not be signed
          double eyp = 0, eym = 0;
          try {
            eyp = variation.second["up"].as<double>();
            //if (variationName.empty()) eyp = fabs(eyp);
          } catch (...) {
            eyp = std::numeric_limits<double>::quiet_NaN();
          }
          try {
            eym = variation.second["dn"].as<double>();
            if (variationName.empty()) eym = fabs(eym);
          } catch (...) {
            eym = std::numeric_limits<double>::quiet_NaN();
          }
          thispoint.setXErrs(eym, eyp, variationName);
        }
      }
      this->_variationsParsed =true;
    }
  }


  // Prepare the variations to be written
  void Scatter1D::writeVariationsToAnnotations(){
    YAML::Emitter em;
    em.SetMapFormat(YAML::Flow);
    em << YAML::BeginMap;
    for (size_t thisPointIndex = 0; thisPointIndex < this->numPoints(); ++thisPointIndex) {
      Point1D& thisPoint = this->_points[thisPointIndex];
      em << YAML::Key << thisPointIndex;
      em << YAML::Value << YAML::BeginMap;
      for (const auto& variation : this->variations()) {
        em << YAML::Key << variation;
        em << YAML::Value << YAML::BeginMap;
        em << YAML::Key <<  "up";
        em << YAML::Value <<  thisPoint.xErrPlus(variation);
        em << YAML::Key <<  "dn";
        em << YAML::Value <<  thisPoint.xErrMinus(variation);
        em << YAML::EndMap;
      }
      em << YAML::EndMap;
    }
    em << YAML::EndMap;
    const std::string val = em.c_str();
    this->setAnnotation("ErrorBreakdown", val);
  }


  /// @todo Reduce duplication between Scatter types
  std::vector<std::string> Scatter1D::variations() const  {
    /// @todo Auto-run parseVariations? Why expose the machinery to the user?
    std::vector<std::string> vecvariations;
    for (auto& point : this->_points) {
      for (auto& it : point.errMap()) {
        // if the variation is not already in the vector, add it!
        if (std::find(vecvariations.begin(), vecvariations.end(), it.first) == vecvariations.end()) {
          vecvariations.push_back(it.first);
        }
      }
    }
    return vecvariations;
  }


  void Scatter1D::rmVariations() {
    _variationsParsed = false;
    for (Point1D& point : this->_points) point.rmVariations();
  }


}
