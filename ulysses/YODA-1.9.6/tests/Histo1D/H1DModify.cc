#include "YODA/Histo1D.h"
#include "YODA/Utils/Formatting.h"

using namespace YODA;
using namespace std;

int main() {
  MSG_BLUE("Testing Histo1D modifiers: ");

  MSG_(PAD(70) << "Setting up a 100 bin histo from 0-100: ");
  Histo1D h(100, 0, 100);
  MSG_GREEN("PASS");

  MSG_(PAD(70) << "Trying to merge bins: ");
  h.mergeBins(0, 10);
  MSG_GREEN("PASS");

  MSG_(PAD(70) << "Testing if bin number was updated: ");
  if (h.numBins() != 90) {
    MSG_RED("FAIL");
    return -1;
  }
  MSG_GREEN("PASS");

  MSG_(PAD(70) << "Checking if bin size was updated: ");
  if (!fuzzyEquals(h.bin(0).xMin(), 0) || !fuzzyEquals(h.bin(0).xMax(), 11)) {
    MSG_RED("FAIL");
    return -1;
  }
  MSG_GREEN("PASS");

  MSG_(PAD(70) << "Checking if bin removal works: ");
  h.eraseBin(0);
  MSG_GREEN("PASS");

  MSG_(PAD(70) << "Was the bin number updated properly? ");
  if (h.numBins() != 89) {
    MSG_RED("FAIL");
    return -1;
  }
  MSG_GREEN("PASS");

  MSG_(PAD(70) << "Was the right bin removed? ");
  if (fuzzyEquals(h.bin(0).xMin(), 0) && fuzzyEquals(h.bin(0).xMax(), 11)) {
    MSG_RED("FAIL");
    return -1;
  }
  MSG_GREEN("PASS");

  MSG_(PAD(70) << "Checking if it is possible to add a bin: ");
  h.addBin(0, 11);
  MSG_GREEN("PASS");

  MSG_(PAD(70) << "Checking if it was added properly: ");
  if (!fuzzyEquals(h.bin(0).xMin(), 0) || !fuzzyEquals(h.bin(0).xMax(),11)) {
    MSG_RED("FAIL");
    return -1;
  }
  MSG_GREEN("PASS");

  MSG_(PAD(70) << "Checking the reset function: ");
  h.reset();
  if (h.integral() != 0) {
    MSG_RED("FAIL");
    return -1;
  }
  MSG_GREEN("PASS");

  return EXIT_SUCCESS;
}
