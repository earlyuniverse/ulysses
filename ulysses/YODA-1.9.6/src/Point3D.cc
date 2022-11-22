#include "YODA/Point3D.h"
#include "YODA/Scatter3D.h"

namespace YODA {


  void Point3D::getVariationsFromParent() const {
    if (this->getParent()) this->getParent<Scatter3D>()->parseVariations();
  }


}
