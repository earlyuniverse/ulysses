#include "YODA/Point2D.h"
#include "YODA/Scatter2D.h"

namespace YODA {


  /// Get error map for direction @a i
  const std::map<std::string, std::pair<double,double>>& Point2D::errMap() const {
    getVariationsFromParent();
    return _ey;
  }


  void Point2D::getVariationsFromParent() const {
    if (this->getParent()) {
      this->getParent<Scatter2D>()->parseVariations();
    }
  }


}
