#include "YODA/Utils/sortedvector.h"
#include <iostream>
#include <cstdlib>
#include <set>

using namespace std;
using namespace YODA;

int main() {
  Utils::sortedvector<int> sv;
  sv.insert(1);
  sv.insert(5);
  sv.insert(3);
  sv.insert(2);
  sv.insert(4);

  for (Utils::sortedvector<int>::const_iterator it = sv.begin(); it != sv.end(); ++it) {
    cout << *it << endl;
  }

  for (int i = 0; i < 5; ++i) {
    if (sv[i] != i+1) return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
