#include "YODA/Point1D.h"
#include "YODA/Scatter1D.h"

namespace YODA {


  void Point1D::getVariationsFromParent() const {
    if (this->getParent()) this->getParent<Scatter1D>()->parseVariations();
  }


}
