#include "yaml-cpp/null.h"

namespace YODA_YAML {
_Null Null;

bool IsNullString(const std::string& str) {
  return str.empty() || str == "~" || str == "null" || str == "Null" ||
         str == "NULL";
}
}
