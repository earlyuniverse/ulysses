#include "YODA/Utils/indexedset.h"
#include <iostream>
#include <cstdlib>
#include <set>

using namespace std;
using namespace YODA;

int main() {
  Utils::indexedset<int> iset;
  iset.insert(1);
  iset.insert(5);
  iset.insert(3);
  iset.insert(2);
  iset.insert(4);

  // for (utils::indexedset<int>::const_iterator it = iset.begin(); it != iset.end(); ++it) {
  //   cout << *it << endl;
  // }

  cout << "iset[3]: " << iset[3] << " == 4: "
       << boolalpha << (iset[3] == 4)
       << endl;

  return (iset[3] == 4) ? EXIT_SUCCESS : EXIT_FAILURE;
}
