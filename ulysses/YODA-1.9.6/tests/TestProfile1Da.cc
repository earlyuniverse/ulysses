#include "YODA/Profile1D.h"
#include "YODA/Utils/Formatting.h"
#include <cmath>
#include <iostream>
#include <unistd.h>

using namespace std;
using namespace YODA;


int main() {

  Profile1D h(20, 0.0, 1.0);
  for (size_t n = 0; n < 10000; ++n) {
    const double x = rand()/static_cast<double>(RAND_MAX);
    const double y = rand()/static_cast<double>(RAND_MAX) * 20 * x;
    h.fill(x, y, 2);
  }

  for (int i = 0; i < 4; ++i) {
    if (i > 0) h.rebin(2);
    MSG("Profile (rebinning #" << i << ", num bins = " << h.numBins() << ")");
    for (const ProfileBin1D& b : h.bins()) {
      MSG(b.xMin() << "-" << b.xMax() << ": "
          << RED(b.mean()) << ", " << BLUE(b.stdDev()) << ", " << RED(b.stdErr()));
    }
  }

  return EXIT_SUCCESS;
}
