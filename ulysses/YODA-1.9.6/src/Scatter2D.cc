#include "YODA/Scatter2D.h"
#include "YODA/Histo1D.h"
#include "YODA/Profile1D.h"
#include <sstream>
#include "yaml-cpp/yaml.h"
#ifdef YAML_NAMESPACE
#define YAML YAML_NAMESPACE
#endif

namespace YODA {

  using namespace std;


  /// Make a Scatter2D representation of a Histo1D
  Scatter2D mkScatter(const Histo1D& h, bool usefocus, bool binwidthdiv,
                      double uflow_binwidth, double oflow_binwidth) {

    Scatter2D rtn;
    for (const std::string& a : h.annotations()) rtn.setAnnotation(a, h.annotation(a));
    rtn.setAnnotation("Type", h.type()); // might override the copied ones

    // Underflow point
    if (uflow_binwidth > 0) {
      const double ex = uflow_binwidth/2;
      const double x = h.xMin() - ex;
      const Dbn1D& hu = h.underflow();

      double y;
      try {
        y = hu.sumW();
      } catch (const Exception&) { // LowStatsError or WeightError
        y = std::numeric_limits<double>::quiet_NaN();
      }
      if (binwidthdiv) y /= uflow_binwidth;
      double ey;
      try {
        ey = hu.relErrW() * y;
      } catch (const Exception&) { // LowStatsError or WeightError
        ey = std::numeric_limits<double>::quiet_NaN();
      }

      Point2D pt(x, y, ex, ex, ey, ey);
      pt.setParent(&rtn);
      rtn.addPoint(pt);
    }

    // In-range points
    for (const HistoBin1D& b : h.bins()) {
      const double x = usefocus ? b.xFocus() : b.xMid();
      const double ex_m = x - b.xMin();
      const double ex_p = b.xMax() - x;

      double y;
      try {
        y = b.sumW();
      } catch (const Exception&) { // LowStatsError or WeightError
        y = std::numeric_limits<double>::quiet_NaN();
      }
      if (binwidthdiv) y /= b.xWidth();
      double ey;
      try {
        ey = b.relErr() * y;
      } catch (const Exception&) { // LowStatsError or WeightError
        ey = std::numeric_limits<double>::quiet_NaN();
      }

      // Attach the point to its parent
      Point2D pt(x, y, ex_m, ex_p, ey, ey);
      pt.setParent(&rtn);
      rtn.addPoint(pt);
    }

    // Overflow point
    if (oflow_binwidth > 0) {
      const double ex = oflow_binwidth/2;
      const double x = h.xMin() - ex;
      const Dbn1D& ho = h.overflow();

      double y;
      try {
        y = ho.sumW();
      } catch (const Exception&) { // LowStatsError or WeightError
        y = std::numeric_limits<double>::quiet_NaN();
      }
      if (binwidthdiv) y /= oflow_binwidth;
      double ey;
      try {
        ey = ho.relErrW() * y;
      } catch (const Exception&) { // LowStatsError or WeightError
        ey = std::numeric_limits<double>::quiet_NaN();
      }

      Point2D pt(x, y, ex, ex, ey, ey);
      pt.setParent(&rtn);
      rtn.addPoint(pt);
    }

    return rtn;
  }


  /// Make a Scatter2D representation of a Profile1D
  Scatter2D mkScatter(const Profile1D& p, bool usefocus, bool usestddev,
                      double uflow_binwidth, double oflow_binwidth) {

    Scatter2D rtn;
    for (const std::string& a : p.annotations()) rtn.setAnnotation(a, p.annotation(a));
    rtn.setAnnotation("Type", p.type());

    // Underflow point
    if (uflow_binwidth > 0) {
      const double ex = uflow_binwidth/2;
      const double x = p.xMin() - ex;
      const Dbn2D& hu = p.underflow();

      double y;
      try {
        y = hu.yMean();
      } catch (const Exception&) { // LowStatsError or WeightError
        y = std::numeric_limits<double>::quiet_NaN();
      }
      double ey;
      try {
        ey = usestddev ? hu.yStdDev() : hu.yStdErr(); ///< Control y-error scheme via usestddev arg
      } catch (const Exception&) { // LowStatsError or WeightError
        ey = std::numeric_limits<double>::quiet_NaN();
      }

      Point2D pt(x, y, ex, ex, ey, ey);
      pt.setParent(&rtn);
      rtn.addPoint(pt);
    }

    // In-range bins
    for (const ProfileBin1D& b : p.bins()) {
      const double x = usefocus ? b.xFocus() : b.xMid();
      const double ex_m = x - b.xMin();
      const double ex_p = b.xMax() - x;

      double y;
      try {
        y = b.mean();
      } catch (const Exception&) { // LowStatsError or WeightError
        y = std::numeric_limits<double>::quiet_NaN();
      }
      double ey;
      try {
        ey = usestddev ? b.stdDev() : b.stdErr(); ///< Control y-error scheme via usestddev arg
      } catch (const Exception&) { // LowStatsError or WeightError
        ey = std::numeric_limits<double>::quiet_NaN();
      }

      //const Point2D pt(x, y, ex_m, ex_p, ey, ey);
      Point2D pt(x, y, ex_m, ex_p, ey, ey);
      pt.setParent(&rtn);
      rtn.addPoint(pt);
    }

    // Overflow point
    if (oflow_binwidth > 0) {
      const double ex = oflow_binwidth/2;
      const double x = p.xMin() - ex;
      const Dbn2D& ho = p.overflow();

      double y;
      try {
        y = ho.yMean();
      } catch (const Exception&) { // LowStatsError or WeightError
        y = std::numeric_limits<double>::quiet_NaN();
      }
      double ey;
      try {
        ey = usestddev ? ho.yStdDev() : ho.yStdErr(); ///< Control y-error scheme via usestddev arg
      } catch (const Exception&) { // LowStatsError or WeightError
        ey = std::numeric_limits<double>::quiet_NaN();
      }

      Point2D pt(x, y, ex, ex, ey, ey);
      pt.setParent(&rtn);
      rtn.addPoint(pt);
    }

    return rtn;
  }


  void Scatter2D::updateTotalUncertainty() {
    for (size_t thisPointIndex = 0; thisPointIndex < this->numPoints(); ++thisPointIndex) {
      Point2D& thisPoint = this->_points[thisPointIndex];
      thisPoint.updateTotalUncertainty();
    }
  }


  // Prepare the variations to be written
  void Scatter2D::writeVariationsToAnnotations() {
    // If there are no variations to write, exit early
    if (variations().empty()) return;
    // There *are* some variations to encode...
    YAML::Emitter em;
    em.SetMapFormat(YAML::Flow);
    em << YAML::BeginMap;
    for (size_t thisPointIndex = 0; thisPointIndex < this->numPoints(); ++thisPointIndex) {
      const Point2D& thisPoint = this->_points[thisPointIndex];
      em << YAML::Key << thisPointIndex;
      em << YAML::Value << YAML::BeginMap;
      for (const auto& variation : this->variations()) {
        em << YAML::Key << variation;
        em << YAML::Value << YAML::BeginMap;
        em << YAML::Key <<  "up";
        em << YAML::Value <<  thisPoint.yErrPlus(variation);
        em << YAML::Key <<  "dn";
        em << YAML::Value <<  thisPoint.yErrMinus(variation);
        em << YAML::EndMap;
      }
      em << YAML::EndMap;
    }
    em << YAML::EndMap;
    const std::string val = em.c_str();
    this->setAnnotation("ErrorBreakdown", val);
  }


  // Retrieve variations from annotation, parse them as YAML, and update the points
  void Scatter2D::parseVariations() {
    if (this->_variationsParsed) return;
    if (!(this->hasAnnotation("ErrorBreakdown"))) return;
    YAML::Node errorBreakdown;
    errorBreakdown = YAML::Load(this->annotation("ErrorBreakdown"));

    if (errorBreakdown.size()) {
      for (size_t thisPointIndex = 0; thisPointIndex < this->numPoints(); ++thisPointIndex) {
        Point2D& thispoint = this->_points[thisPointIndex];
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
          thispoint.setYErrs(eym, eyp, variationName);
        }
      }
      this->_variationsParsed =true;
    }
  }


  /// @todo Reduce duplication between Scatter types
  std::vector<std::string> Scatter2D::variations() const {
    std::vector<std::string> vecVariations;
    for (auto& point : this->_points) {
      for (auto& it : point.errMap()) {
        // if the variation is not already in the vector, add it!
        if (std::find(vecVariations.begin(), vecVariations.end(), it.first) == vecVariations.end()) {
          vecVariations.push_back(it.first);
        }
      }
    }
    return vecVariations;
  }


  void Scatter2D::rmVariations() {
    _variationsParsed = false;
    for (Point2D& point : this->_points) point.rmVariations();
  }


  std::vector<std::vector<double> > Scatter2D::covarianceMatrix(bool ignoreOffDiagonalTerms) {
    int nPoints = this->numPoints();
    //double covM[nPoints][nPoints] = {};
    std::vector<std::vector<double> > covM;


    // initialise cov matrix to be the right shape
    for (int i = 0; i < nPoints; i++) {
      std::vector<double> row;
      row.resize(nPoints);
      covM.push_back(row);
    }

    // case where only have nominal, ie total uncertainty, labelled "" (empty string)
    if (this->variations().size() == 1) {
      for (int i = 0; i < nPoints; i++) {
        covM[i][i] = sqr(((this->_points[i].yErrs().first+this->_points[i].yErrs().second)/2));
        if (covM[i][i] == 0 ) covM[i][i] = 1;
      }
      return covM;
    }
    // more interesting case where we actually have some uncertainty breakdown!
    auto  systList= this->variations();
    for (auto sname : systList){
      if (sname.length() == 0) continue;
      std::vector<double> systErrs;
      systErrs.resize(nPoints);
      for (int i = 0; i < nPoints; i++) {
        auto point = this->_points[i];
        try {
          auto variations = point.errMap().at(sname);
          // up/dn are symmetrized since this method can't handle asymmetric errors
          systErrs[i] = (fabs(variations.first)+fabs(variations.second))*0.5;
        } catch (const std::exception& e) { // missing bin
          systErrs[i] = 0.0;
        }
      }
      if (ignoreOffDiagonalTerms || sname.find("stat") != std::string::npos || sname.find("uncor") != std::string::npos) {
        for (int i = 0; i < nPoints; i++) {
          covM[i][i] += systErrs[i]*systErrs[i]; // just the diagonal; bins are considered uncorrelated
        }
      } else {
        for (int i = 0; i < nPoints; i++) {
          for (int j = 0; j < nPoints; j++) {
            covM[i][j] += systErrs[i]*systErrs[j];
          }
        }
      }
    }
    return covM;
  }


}
