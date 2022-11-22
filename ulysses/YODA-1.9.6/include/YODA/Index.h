#ifndef YODA_INDEX_H
#define YODA_INDEX_H

#include <sstream>
#include <unordered_map>

namespace YODA {

  /// @class Index Index.h "YODA/Index.h"
  /// @brief Index of a file
  /// 
  /// This class holds index of a file in form of nested unordered_map.
  class Index {
  public:

    /// @brief Alias for nested unordered_map for index.
    ///
    /// 
    /// Structure:
    /// {
    ///     "Scatter2D": {
    ///         "/ATLAS_2017_I1514251/d09-x01-y01": 3
    ///     }
    /// }
    using AOIndex =
      std::unordered_map<std::string, std::unordered_map<std::string, int>>;

    Index() = default;

    Index(AOIndex&& idx) noexcept : _index(idx) {}

    Index(const AOIndex& idx) : _index(idx) {}

    Index(Index&& other) : _index(std::move(other._index)) {}

    Index& operator=(const Index& rhs) {
      _index = rhs._index;
      return *this;
    }

    Index& operator=(Index&& rhs) {
      _index = std::move(rhs._index);
      return *this;
    }

    ///@brief Implicitly cast to a nested index map.
    operator AOIndex() const {return _index; }

    ///@brief Get nested index map.
    AOIndex getAOIndex() const { return _index; }

    /// @brief Get string representation of index.
    std::string toString() const {
      std::ostringstream indexStr;
      for (const auto& kv : this->_index) {
        indexStr << "OBJECT TYPE: " << kv.first << "\n";

        for (const auto& path_to_bincount : kv.second) {
          indexStr << "    ----------\n";
          indexStr << "    "
                   << "PATH:      " << path_to_bincount.first << "\n";
          indexStr << "    "
                   << "BIN COUNT: " << path_to_bincount.second << "\n";
          indexStr << "    ----------\n";
        }
      }
      return indexStr.str();
    }

  private:
    /// @brief Holds index.
    AOIndex _index;
  };

}

#endif