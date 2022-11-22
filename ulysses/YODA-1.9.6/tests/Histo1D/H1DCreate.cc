#include "YODA/Histo1D.h"
#include "YODA/Utils/Formatting.h"

using namespace YODA;
using namespace std;


int main() {
  MSG_BLUE("Testing Histo1D constructors: ");

  MSG_(PAD(70) << "The most basic, linear constructor:");
  Histo1D h(100, 0, 100);
  if (h.numBins() != 100) {
    MSG_RED("FAIL: Wrong number of bins was created!");
    return -1;
  }
  if (h.xMin() != 0) {
    MSG_RED("FAIL: Low edge wasn't properly set!");
    return -1;
  }
  if (h.xMax() != 100) {
    MSG_RED("FAIL: High edge wasn't properly set!");
    return -1;
  }
  if (h.integral() != 0) {
    MSG_RED("FAIL: The constructor is setting some statistics!");
    return -1;
  }
  MSG_GREEN("PASS");


  MSG_(PAD(70) << "Explicit bin edges constructor: ");
  vector<double> edges;
  for (int i = 0; i < 101; ++i) edges.push_back(i);
  Histo1D h1(edges);
  MSG_GREEN("FOO1");
  if (h1.numBins() != 100) {
    MSG_RED("FAIL: Wrong number of bins was created!");
    return -1;
  }
  MSG_GREEN("FOO2");
  if (h1.xMin() != 0) {
    MSG_RED("FAIL: Low edge wasn't properly set!");
    return -1;
  }
  MSG_GREEN("FOO3");
  if (h1.xMax() != 100) {
    MSG_RED("FAIL: High edge wasn't properly set!");
    return -1;
  }
  MSG_GREEN("FOO4");
  if (h1.integral() != 0) {
    MSG_RED("FAIL: The constructor is setting some statistics!");
    return -1;
  }
  MSG_GREEN("PASS");


  MSG_(PAD(70) << "Copy constructor: ");
  Histo1D h2(h);
  if (h2.numBins() != 100) {
    MSG_RED("FAIL: Wrong number of bins was created!");
    return -1;
  }
  if (h2.xMin() != 0) {
    MSG_RED("FAIL: Low edge wasn't properly set!");
    return -1;
  }
  if (h2.xMax() != 100) {
    MSG_RED("FAIL: High edge wasn't properly set!");
    return -1;
  }
  if (h2.integral() != 0) {
    MSG_RED("FAIL: The constructor is setting some statistics!");
    return -1;
  }
  MSG_GREEN("PASS");

  return EXIT_SUCCESS;
}
