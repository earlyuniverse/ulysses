#include "YODA/Histo1D.h"
#include "YODA/Utils/Formatting.h"
#include <cmath>
#include <iostream>

using namespace std;
using namespace YODA;


int main() {

  Histo1D h(20, 0.0, 1.0);
  for (size_t n = 0; n < 1000; ++n) {
    const double num = rand()/static_cast<double>(RAND_MAX);
    h.fill(num);
  }

  MSG("Path = " << h.path());
  MSG("Mean value = " << h.xMean() << " +- " << h.xStdErr());
  MSG("Total area = " << h.integral());

  auto compareHeight = [](const HistoBin1D& a, const HistoBin1D& b) { return a.height() < b.height(); };
  const HistoBin1D& highestBin = *( max_element(h.bins().begin(), h.bins().end(), compareHeight) );
  const double maxHeight = highestBin.height();

  for (int i = 0; i < 4; ++i) {
    if (i > 0) h.rebin(2);
    MSG("Histo (rebinning #" << i << ", num bins = " << h.numBins() << ")");
    for (const HistoBin1D& b : h.bins()) {
      const int numElements = static_cast<int>(round(20 * b.height()/maxHeight));
      MSG(string().insert(0, numElements, '=') << "  " << RED(b.height()));
    }
  }

  return EXIT_SUCCESS;
}
