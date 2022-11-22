#include "YODA/Scatter2D.h"
#include "YODA/Utils/Formatting.h"

using namespace YODA;
using namespace std;

int main() {
  MSG_BLUE("Testing Scatter2D construction: ");


  MSG_(PAD(70) << "Constructing a scatter (empty const): ");
  Scatter2D s1;
  MSG_GREEN("PASS");

  MSG_(PAD(70) << "Constructing a scatter (vector of points)");
  vector<Point2D> points;
  Point2D apoint(0,0,0);
  points.push_back(apoint);
  Scatter2D s2(points);
  MSG_GREEN("PASS");

  MSG_(PAD(70) << "Constructing a scatter (values, no errs) ");
  vector<double> values;
  values.push_back(0);
  Scatter2D s3(values, values);
  MSG_GREEN("PASS");

  MSG_(PAD(70) << "Constructing a scatter (values, sym errs)");
  Scatter2D s4(values, values, values, values);
  MSG_GREEN("PASS");

  vector<pair<double, double> > valuesS;
  valuesS.push_back(make_pair(0,0));

  // MSG_(PAD(70) << "Constructing a scatter (sym err x,asym y)");
  // Scatter2D s5(values, values, values, valuesS);
  // MSG_GREEN("PASS");

  // MSG_(PAD(70) << "Constructing a scatter (asym x, sym y) ");
  // Scatter2D s6(values, values, valuesS, values);
  // MSG_GREEN("PASS");

  MSG_(PAD(70) << "Constructing a scatter (asym x, asym y) ");
  Scatter2D s7(values, values, valuesS, valuesS);
  MSG_GREEN("PASS");

  MSG_(PAD(70) << "Constructing a scatter (explicit asym) ");
  Scatter2D s8(values, values, values, values, values, values);
  MSG_GREEN("PASS");

  MSG_(PAD(70) << "Testing a copy operator: ");
  Scatter2D s9(s8);
  MSG_GREEN("PASS");

  MSG_(PAD(70) << "Testing an assignment operator: ");
  Scatter2D s10(s9);
  s10 = s7;
  MSG_GREEN("PASS");

  return EXIT_SUCCESS;
}
