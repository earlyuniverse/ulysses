#include "YODA/Profile1D.h"
#include "YODA/Utils/Formatting.h"

using namespace YODA;
using namespace std;

int main() {
  MSG_BLUE("Testing Profile1D filling: ");

  MSG_(PAD(70) << "Setting up 100-bin profile histo ");
  Profile1D p(100, 0, 100);
  MSG_GREEN("PASS");

  MSG_(PAD(70) << "Testing fill operator: ");
  p.fill(1,1,2);
  if (p.sumW() != 2 || p.sumW2() != 4) {
    MSG_RED("FAIL");
    return -1;
  }
  MSG_GREEN("PASS");


  MSG_(PAD(70) << "Testing the fill of the underflow: ");
  p.fill(-10,2,3);
  if (p.underflow().xMean() != -10) {
    MSG_RED("FAIL");
    return -1;
  }
  if (p.underflow().yMean() != 2) {
    MSG_RED("FAIL");
    return -1;
  }
  if (p.underflow().sumW() != 3 || p.underflow().sumW2() != 9) {
    MSG_RED("FAIL");
    return -1;
  }
  MSG_GREEN("PASS");


  MSG_(PAD(70) << "Testing the fill of the underflow: ");
  p.fill(110,2,3);
  if (p.overflow().xMean() != 110) {
    MSG_RED("FAIL");
    return -1;
  }
  if (p.overflow().yMean() != 2) {
    MSG_RED("FAIL");
    return -1;
  }
  if (p.overflow().sumW() != 3 || p.overflow().sumW2() != 9) {
    MSG_RED("FAIL");
    return -1;
  }
  MSG_GREEN("PASS");

  return EXIT_SUCCESS;
}
