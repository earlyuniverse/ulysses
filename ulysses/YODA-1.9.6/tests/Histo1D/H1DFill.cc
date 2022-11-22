#include "YODA/Histo1D.h"
#include "YODA/Utils/Formatting.h"
#include <cmath>
#include <limits>

using namespace YODA;
using namespace std;

int main() {
  MSG_BLUE("Testing Histo1D filling: ");

  MSG_(PAD(70) << "Setting up a 100 bin histo from 0-100: ");
  Histo1D h(100, 0, 100);
  MSG_GREEN("PASS");

  MSG_(PAD(70) << "Trying to fill the sample histogram: ");
  h.fill(0, 2);
  MSG_GREEN("PASS");

  MSG_(PAD(70) << "Checking sumW: ");
  if (h.sumW() != 2) {
    MSG_RED("FAIL");
    return -1;
  }
  MSG_GREEN("PASS");

  MSG_(PAD(70) << "Checking sumW2: ");
  if (h.sumW2() != 4) {
    MSG_RED("FAIL");
    return -1;
  }
  MSG_GREEN("PASS");

  MSG_(PAD(70) << "Trying to fill again: ");
  h.fill(10, 2);
  MSG_GREEN("PASS");

  MSG_(PAD(70) << "Checking mean: ");
  if (!fuzzyEquals(5, h.xMean(false))) {
    MSG_RED("FAIL");
    return -1;
  }
  MSG_GREEN("PASS");

  MSG_(PAD(70) << "Checking variance (total): ");
  if (!fuzzyEquals(50, h.xVariance())) { //< Note the effective / (N-1) for unbiasing
    MSG_RED("FAIL");
    return -1;
  }
  MSG_GREEN("PASS");

  MSG_(PAD(70) << "Checking variance (sum): ");
  if (!fuzzyEquals(50, h.xVariance(false))) { //< Note the effective / (N-1) for unbiasing
    MSG_RED("FAIL");
    return -1;
  }
  MSG_GREEN("PASS");

  MSG_(PAD(70) << "Checking standard deviation (total): ");
  if (!fuzzyEquals(sqrt(50), h.xStdDev())) {
    MSG_RED("FAIL");
    return -1;
  }
  MSG_GREEN("PASS");

  MSG_(PAD(70) << "Checking standard deviation (sum): ");
  if (!fuzzyEquals(sqrt(50), h.xStdDev(false))) {
    MSG_RED("FAIL");
    return -1;
  }
  MSG_GREEN("PASS");

  MSG_(PAD(70) << "Checking bin index lookup: ");
  const int i = h.binIndexAt(74.5);
  if (i != 74) {
    MSG_RED("FAIL");
    return -1;
  }
  MSG_GREEN("PASS");

  MSG_(PAD(70) << "Checking integral: ");
  const double integral = h.integral((size_t) i);
  MSG_BLUE(integral);
  if (!fuzzyEquals(integral, 4)) {
    MSG_RED("FAIL");
    return -1;
  }
  MSG_GREEN("PASS");

  MSG_(PAD(70) << "Trying to fill the underflow: ");
  h.fill(-10,1);
  MSG_GREEN("PASS");

  MSG_(PAD(70) << "Checking if stats were set correctly: ");
  if (!fuzzyEquals(h.underflow().xMean(), -10)) {
    MSG_RED("FAIL");
    return -1;
  }
  MSG_GREEN("PASS");

  MSG_(PAD(70) << "Trying to fill the overflow: ");
  h.fill(110,1);
  MSG_GREEN("PASS");

  MSG_(PAD(70) << "Checking if stats were set correctly: ");
  if (!fuzzyEquals(h.overflow().xMean(), 110)) {
    MSG_RED("FAIL");
    return -1;
  }
  MSG_GREEN("PASS");


  // MSG_(PAD(70) << "Trying to fill with inf: ");
  // h.fill(numeric_limits<double>::infinity(),1);
  // MSG_GREEN("PASS");

  // MSG_(PAD(70) << "Checking if stats were set correctly: ");
  // if (isinf(h.xMean())) {
  //   MSG_RED("FAIL");
  //   return -1;
  // }
  // if (isinf(h.xMean(false))) {
  //   MSG_RED("FAIL");
  //   return -1;
  // }
  // MSG_GREEN("PASS");


  // MSG_(PAD(70) << "Trying to fill with NaN: ");
  // h.fill(numeric_limits<double>::quiet_NaN(),1);
  // MSG_GREEN("PASS");

  // MSG_(PAD(70) << "Checking if stats were set correctly: ");
  // if (isnan(h.xMean())) {
  //   MSG_RED("FAIL");
  //   return -1;
  // }
  // MSG_GREEN("PASS");

  return EXIT_SUCCESS;
}
