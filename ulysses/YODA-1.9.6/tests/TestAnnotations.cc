#include "YODA/Counter.h"
#include <iostream>

using namespace std;
using namespace YODA;

int main() {

  Counter c("/blah", "Foo bar baz");
  c.setAnnotation("Foo", "1.234");
  c.setAnnotation("Bar", 5.678);

  if (c.title() != "Foo bar baz") return 1;

  if (c.annotation("Foo") != "1.234") return 2;
  if (c.annotation<int>("Foo") != 1) return 3;
  if (c.annotation<double>("Foo") != 1.234) return 4;

  // if (c.annotation("Bar") != "5.678") return 5;
  if (c.annotation<int>("Bar") != 5) return 6;
  if (c.annotation<double>("Bar") != 5.678) return 7;

  return 0;
}
