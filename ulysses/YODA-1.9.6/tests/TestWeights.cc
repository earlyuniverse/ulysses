#include "YODA/Weights.h"
#include <iostream>

using namespace std;
using namespace YODA;


int foo(const Weights& w) {
  return w.size();
}


int main() {

  vector<string> keys;
  keys.push_back("FOO");
  keys.push_back("BAR");
  keys.push_back("BAZ");

  Weights w1(keys);
  w1["FOO"] += 1;
  w1["BAR"] += 2;
  w1["BAZ"] += 3.4;
  w1 *= 3;
  cout << "W1 = " << w1 << endl;

  /// @todo Allow construction from a literal C array of key names
  Weights w2(keys, 1.0);
  w2 += w1;
  cout << "W2 = " << w2 << endl;

  Weights w3 = w1 / w2;
  cout << "W3 = " << w3 << endl;

  Weights w4 = 5;
  w4 += 4.5;
  w4 += -1.5;
  w4 += -1.5;
  foo(w1);
  foo(w2);
  foo(w3);
  foo(w4);
  foo(5);
  foo(Weights(6));

  return 0;
}
