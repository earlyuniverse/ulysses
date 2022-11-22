#include "YODA/Histo2D.h"
#include "YODA/Utils/Formatting.h"

#include <iostream>
using namespace std;
using namespace YODA;

int main() {
  MSG_BLUE("Testing Histo2D creation: ");

  MSG_(PAD(70) << "Setting up two 10x10 bin histos from 0-100 x 0-100: ");
  Histo2D first(10, 0, 100, 10, 0, 100);
  Histo2D second(10, 0, 100, 10, 0, 100);
  MSG_GREEN("PASS");

  MSG_(PAD(70) << "Filling the histos: ");
  first.fill(1, 1, 1);
  first.fill(23, 1, 1);
  second.fill(1, 1, 1);
  MSG_GREEN("PASS");

  MSG_(PAD(70) << "Testing the copy constructor: ");
  Histo2D copyTest(first);
  MSG_GREEN("PASS");

  MSG_(PAD(70) << "Adding histos: ");
  Histo2D added(first + second);
  MSG_GREEN("PASS");

  MSG_(PAD(70) << "Subtracting histos: ");
  Histo2D subtracted(first - second);
  MSG_GREEN("PASS");

  // MSG_(PAD(70) << "Dividing histos: ");
  // Scatter3D divided(first / second);
  // MSG_GREEN("PASS");

  return EXIT_SUCCESS;
}
