#include<vector>
#include<iostream>
#include<memory>
using namespace std;
#if __cplusplus <= 199711L
  #error This library needs at least a C++11 compliant compiler
#endif
int main() {
  vector<double> vv(5,5.0);
  for (auto & v: vv) {
    cout << v << endl;
  }
}
