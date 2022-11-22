#include "YODA/Scatter2D.h"
#include "YODA/Utils/Formatting.h"

using namespace YODA;
using namespace std;

int main() {
  MSG_BLUE("Testing Scatter2D modifiers: ");


  MSG_(PAD(70) << "Constructing a scatter: ");
  vector<double> coords;
  coords.push_back(0); coords.push_back(1);
  Scatter2D s1(coords, coords);
  MSG_GREEN("PASS");


  MSG_(PAD(70) << "Scaling the scatter: ");
  s1.scaleXY(2, 3);
  if (s1.point(1).x() != 2 || s1.point(1).y() != 3) {
    MSG_RED("FAIL");
    return -1;
  }
  MSG_GREEN("PASS");


  MSG_(PAD(70) << "Adding a point (first method): ");
  Point2D point(1, 1);
  s1.addPoint(point);
  if (s1.numPoints() != 3) {
    MSG_RED("FAIL");
    return -1;
  }
  if (s1.point(1).x() != 1 || s1.point(1).y() != 1) {
    MSG_RED("FAIL");
    return -1;
  }
  MSG_GREEN("PASS");


  MSG_(PAD(70) << "Adding a point (second method): ");
  s1.addPoint(-1, -1);
  if (s1.numPoints() != 4) {
    MSG_RED("FAIL");
    return -1;
  }
  if (s1.point(0).x() != -1 || s1.point(0).y() != -1) {
    MSG_RED("FAIL");
    return -1;
  }
  MSG_GREEN("PASS");


  MSG_(PAD(70) << "Adding a point (third method): ");
  s1.addPoint(5, 4, 6, 3);
  if (s1.numPoints() != 5) {
    MSG_RED("FAIL");
    return -1;
  }
  if (s1.point(4).x() != 5 || s1.point(4).y() != 4 ||
      s1.point(4).xErrMinus() != 6 || s1.point(4).xErrPlus() != 6 ||
      s1.point(4).yErrMinus() != 3 || s1.point(4).yErrPlus() != 3){
    MSG_RED("FAIL");
    return -1;
  }
  MSG_GREEN("PASS");

  MSG_(PAD(70) << "Adding a point (fourth method): ");
  pair<double, double> errs = make_pair(1.0, 2.0);
  s1.addPoint(10, 11, errs, errs);
  if (s1.numPoints() != 6) {
    MSG_RED("FAIL");
    return -1;
  }
  if (s1.point(5).x() != 10 || s1.point(5).y() != 11 ||
      s1.point(5).xErrMinus() != 1 || s1.point(5).xErrPlus() != 2 ||
      s1.point(5).yErrMinus() != 1 || s1.point(5).yErrPlus() != 2) {
    MSG_RED("FAIL");
    return -1;
  }
  MSG_GREEN("PASS");


  MSG_(PAD(70) << "Adding a point (fifth method): ");
  s1.addPoint(12, 14, make_pair(6, 5), errs);
  if (s1.numPoints() != 7) {
    MSG_RED("FAIL");
    return -1;
  }
  if (s1.point(6).x() != 12 || s1.point(6).y() != 14 ||
      s1.point(6).xErrMinus() != 6 || s1.point(6).xErrPlus() != 5 ||
      s1.point(6).yErrMinus() != 1 || s1.point(6).yErrPlus() != 2){
    MSG_RED("FAIL");
    return -1;
  }
  MSG_GREEN("PASS");


  MSG_(PAD(70) << "Adding a point (sixth method): ");
  pair<double, double> errsy;
  errsy = make_pair(100, 200);
  s1.addPoint(300, 400, errs, errsy);
  if (s1.numPoints() != 8) {
    MSG_RED("FAIL");
    return -1;
  }
  if (s1.point(7).x() != 300 || s1.point(7).y() != 400 ||
      s1.point(7).xErrMinus() != 1 || s1.point(7).xErrPlus() != 2 ||
      s1.point(7).yErrMinus() != 100 || s1.point(7).yErrPlus() != 200){
    MSG_RED("FAIL");
    return -1;
  }
  MSG_GREEN("PASS");


  MSG_(PAD(70) << "Adding a point (seventh method): ");
  s1.addPoint(300, 400, 1, 2, 3, 4);
  if (s1.numPoints() != 9) {
    MSG_RED("FAIL");
    return -1;
  }

  if (s1.point(8).x() != 300 || s1.point(8).y() != 400 ||
      s1.point(8).xErrMinus() != 1 || s1.point(8).xErrPlus() != 2 ||
      s1.point(8).yErrMinus() != 3 || s1.point(8).yErrPlus() != 4){
    MSG_RED("FAIL");
    return -1;
  }
  MSG_GREEN("PASS");


  MSG_(PAD(70) << "Adding a point (eighth method): ");
  Point2D p2(800, 900); Point2D p3(1000, 1000);
  vector<Point2D> p5; p5.push_back(p2); p5.push_back(p3);
  s1.addPoints(p5);
  if (s1.numPoints() != 11) {
    MSG_RED("FAIL");
    return -1;
  }
  if (s1.point(9).x() != 800 || s1.point(9).y() != 900 ||
      s1.point(9).xErrMinus() != 0 || s1.point(9).xErrPlus() != 0 ||
      s1.point(9).yErrMinus() != 0 || s1.point(9).yErrPlus() != 0 ||
      s1.point(10).x() != 1000 || s1.point(10).y() != 1000){
    MSG_RED("FAIL");
    return -1;
  }
  MSG_GREEN("PASS");


  MSG_(PAD(70) << "Trying to reset the scatter: ");
  s1.reset();
  if (s1.numPoints() != 0){
    MSG_RED("FAIL");
    return -1;
  }
  MSG_GREEN("PASS");

  return EXIT_SUCCESS;
}
