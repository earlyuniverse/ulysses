#include "YODA/Histo1D.h"
#include "YODA/WriterYODA.h"
#include <cmath>
#include <vector>
#include <iostream>
#include <memory>
#include <type_traits>

using namespace std;
using namespace YODA;


int main() {

  //shared_ptr<AnalysisObject>
  auto h = make_shared<Histo1D>(20, 0.0, 1.0);
  for (size_t n = 0; n < 1000; ++n) {
    const double num = rand()/static_cast<double>(RAND_MAX);
    h->fill(num);
  }
  WriterYODA::create().write("testwriter1.yoda", h);
  WriterYODA::write("testwriter1.yoda", h);

  vector< shared_ptr<Histo1D> > hs;
  hs.push_back(h);
  WriterYODA::write("testwriter2.yoda", hs);

  vector< shared_ptr<AnalysisObject> > aos;
  aos.push_back( static_pointer_cast<AnalysisObject>(h) );
  WriterYODA::write("testwriter2.yoda", aos);
  WriterYODA::write("testwriter2.yoda.gz", aos);

  return EXIT_SUCCESS;
}
