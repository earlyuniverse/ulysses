#include "YODA/Scatter3D.h"
#include "YODA/Histo2D.h"
#include "YODA/Profile2D.h"
#include "YODA/Exceptions.h"
#include <sstream>
#include "yaml-cpp/yaml.h"
#ifdef YAML_NAMESPACE
#define YAML YAML_NAMESPACE
#endif

namespace YODA {


  Scatter3D mkScatter(const Histo2D& h, bool usefocus, bool binareadiv) {
    Scatter3D rtn;
    for (const std::string& a : h.annotations()) rtn.setAnnotation(a, h.annotation(a));
    rtn.setAnnotation("Type", h.type());

    for (size_t i = 0; i < h.numBins(); ++i) {
      const HistoBin2D& b = h.bin(i);

      /// SAME FOR ALL 2D BINS

      double x = b.xMid();
      if (usefocus) {
        try {
          x = b.xFocus();
        } catch (const LowStatsError& lse) {
          x = b.xMid();
        }
      }
      const double exminus = x - b.xMin();
      const double explus = b.xMax() - x;

      double y = b.yMid();
      if (usefocus) {
        try {
          y = b.yFocus();
        } catch (const LowStatsError& lse) {
          y = b.yMid();
        }
      }
      const double eyminus = y - b.yMin();
      const double eyplus = b.yMax() - y;

      /// END SAME FOR ALL 2D BINS

      double z;
      try {
        z = b.sumW();
      } catch (const Exception&) { // LowStatsError or WeightError
        z = std::numeric_limits<double>::quiet_NaN();
      }
      if (binareadiv) z /= b.xWidth()*b.yWidth();
      double ez;
      try {
        ez = b.relErr() * z;
      } catch (const Exception&) { // LowStatsError or WeightError
        ez = std::numeric_limits<double>::quiet_NaN();
      }

      Point3D pt(x, y, z, exminus, explus, eyminus, eyplus, ez, ez);
      pt.setParent(&rtn);
      rtn.addPoint(pt);
    }

    assert(h.numBins() == rtn.numPoints());
    return rtn;
  }


  Scatter3D mkScatter(const Profile2D& h, bool usefocus, bool usestddev) {
    Scatter3D rtn;
    for (const std::string& a : h.annotations())
      rtn.setAnnotation(a, h.annotation(a));
    rtn.setAnnotation("Type", h.type());
    for (size_t i = 0; i < h.numBins(); ++i) {
      const ProfileBin2D& b = h.bin(i);

      /// SAME FOR ALL 2D BINS

      double x = b.xMid();
      if (usefocus) {
        try {
          x = b.xFocus();
        } catch (const LowStatsError& lse) {
          x = b.xMid();
        }
      }
      const double exminus = x - b.xMin();
      const double explus = b.xMax() - x;

      double y = b.yMid();
      if (usefocus) {
        try {
          y = b.yFocus();
        } catch (const LowStatsError& lse) {
          y = b.yMid();
        }
      }
      const double eyminus = y - b.yMin();
      const double eyplus = b.yMax() - y;

      /// END SAME FOR ALL 2D BINS

      double z;
      try {
        z = b.mean();
      } catch (const LowStatsError& lse) {
        z = std::numeric_limits<double>::quiet_NaN();
      }
      double ez;
      try {
        ez = usestddev ? b.stdDev() : b.stdErr(); ///< Control z-error scheme via usestddev arg
      } catch (const LowStatsError& lse) {
        ez = std::numeric_limits<double>::quiet_NaN();
      }

      rtn.addPoint(x, y, z, exminus, explus, eyminus, eyplus, ez, ez);
    }

    return rtn;
  }


  // Prepare the variations to be written
  void Scatter3D::writeVariationsToAnnotations(){
    YAML::Emitter em;
    em.SetMapFormat(YAML::Flow);
    em << YAML::BeginMap;
    for (size_t thisPointIndex = 0; thisPointIndex < this->numPoints(); ++thisPointIndex) {
      Point3D& thisPoint = this->_points[thisPointIndex];
      em << YAML::Key << thisPointIndex;
      em << YAML::Value << YAML::BeginMap;
      for (const auto& variation : this->variations()) {
        em << YAML::Key << variation;
        em << YAML::Value << YAML::BeginMap;
        em << YAML::Key <<  "up";
        em << YAML::Value <<  thisPoint.zErrPlus(variation);
        em << YAML::Key <<  "dn";
        em << YAML::Value <<  thisPoint.zErrMinus(variation);
        em << YAML::EndMap;
      }
      em << YAML::EndMap;
    }
    em << YAML::EndMap;
    const std::string val = em.c_str();
    this->setAnnotation("ErrorBreakdown", val);
  }


  void Scatter3D::updateTotalUncertainty() {
    for (size_t thisPointIndex = 0; thisPointIndex < this->numPoints(); ++thisPointIndex) {
      Point3D& thisPoint = this->_points[thisPointIndex];
      thisPoint.updateTotalUncertainty();
    }
  }


  void Scatter3D::parseVariations() {
    if (this->_variationsParsed) return;
    if (!(this->hasAnnotation("ErrorBreakdown"))) return;
    YAML::Node errorBreakdown;
    errorBreakdown = YAML::Load(this->annotation("ErrorBreakdown"));
    if (errorBreakdown.size()) {
      for (size_t thisPointIndex = 0; thisPointIndex < this->numPoints(); ++thisPointIndex) {
        Point3D& thispoint = this->_points[thisPointIndex];
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
          thispoint.setZErrs(eym, eyp, variationName);
        }
      }
      this->_variationsParsed = true;
    }
  }


  /// @todo Reduce duplication between Scatter types
  std::vector<std::string> Scatter3D::variations() const  {
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


  void Scatter3D::rmVariations() {
    _variationsParsed = false;
    for (Point3D& point : this->_points) point.rmVariations();
  }


}
