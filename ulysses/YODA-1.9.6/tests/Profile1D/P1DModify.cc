#include "YODA/Profile1D.h"
#include "YODA/Utils/Formatting.h"

using namespace YODA;
using namespace std;

int main() {
  MSG_BLUE("Testing Profile1D modifiers: ");

  MSG_(PAD(70) << "Creating the Profile1D: ");
  Profile1D p(100,0,100);
  p.fill(1,1,2);
  MSG_GREEN("PASS");

  MSG_(PAD(70) << "Scaling the height: ");
  p.scaleW(3);
  if (p.sumW() != 6 || p.sumW2() != 36) {
    MSG_RED("FAIL");
    return -1;
  }
  MSG_GREEN("PASS");

  MSG_(PAD(70) << "Resetting the profile: ");
  p.reset();
  if (p.sumW() != 0 || p.sumW2() != 0){
    MSG_RED("FAIL");
    return -1;
  }
  MSG_GREEN("PASS");

  MSG_(PAD(70) << "Merging the bins: ");
  p.mergeBins(0,10);
  if (p.bin(0).xMin() != 0 || p.bin(0).xMax() != 11){
    MSG_RED("FAIL");
    return -1;
  }
  if (p.numBins() != 90){
    MSG_RED("FAIL");
    return -1;
  }
  MSG_GREEN("PASS");

  MSG_(PAD(70) << "Testing rebinning: ");
  p.rebin(2);
  for (size_t i = 1; i < p.bins().size() - 1; ++i){
    if (2 != p.bin(i).xWidth()){
      MSG_RED("FAIL");
      return -1;
    }
  }
  if (p.numBins() != 45) {
    MSG_RED("FAIL");
    return -1;
  }
  MSG_GREEN("PASS");

  MSG_(PAD(70) << "Trying to add a bin (first method): ");
  p.addBin(110, 120);
  if (p.numBins() != 46){
    MSG_RED("FAIL");
    return -1;
  }
  MSG_GREEN("PASS");

  MSG_(PAD(70) << "Trying to add a bin (second method): ");
  vector<double> test;
  test.push_back(120); test.push_back( 140); test.push_back(145);
  p.addBins(test);
  if (p.numBins() != 48){
    MSG_RED("FAIL");
    return -1;
  }
  MSG_GREEN("PASS");

  // MSG_(PAD(70) << "Trying to add a bin (third method): ");
  // vector<pair<double,double> > test2;
  // test2.push_back(make_pair(180,190));
  // p.addBins(test2);
  // if(p.numBins() != 49){
  //     MSG_RED("FAIL");
  //   return -1;
  // }
  // MSG_GREEN("PASS");

  return EXIT_SUCCESS;
}
