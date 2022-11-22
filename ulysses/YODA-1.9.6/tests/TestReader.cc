#include "YODA/Histo1D.h"
#include "YODA/ReaderYODA.h"
#include "YODA/IO.h"
#include <iostream>

using namespace std;
using namespace YODA;

int main() {

  vector<AnalysisObject*> aos1 = YODA::read("testwriter2.yoda");

  Index idx = YODA::ReaderYODA::create().mkIndex("testwriter2.yoda");
  cout << "idx 1:" << idx.toString() << endl;
  Index idx2= YODA::mkIndex("testwriter2.yoda");
  cout << "idx 2:" << idx2.toString() << endl;

  #ifdef WITH_ZLIB
  vector<AnalysisObject*> aos2 = YODA::read("testwriter2.yoda.gz");
  #endif

  return EXIT_SUCCESS;
}
