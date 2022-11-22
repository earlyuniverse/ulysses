// -*- C++ -*-
//
// This file is part of YODA -- Yet more Objects for Data Analysis
// Copyright (C) 2008-2021 The YODA collaboration (see AUTHORS for details)
//
#ifndef YODA_Weights_h
#define YODA_Weights_h

#include "YODA/Exceptions.h"
#include <vector>
#include <string>
#include <map>
#include <ostream>

namespace YODA {


  /// @brief A named, vectorised generalisation of an event weight.
  ///
  /// @todo Accept general Boost.Ranges as constructor args... but start with literal arrays for convenience
  /// @todo Autogenerate numerical names if not given
  class Weights {
  public:

    /// @name Constructors
    /// @{

    Weights(const Weights& other)
      : _values(other._values)
    {  }

    /// Convenience auto-constructor from a single double, since that's the commonest use case.
    Weights(double value) {
      _values["0"] = value;
    }

    /// Constructor from a vector of key/value pairs
    Weights(const std::vector<std::pair<std::string, double> >& keys_values) {
      for (std::vector<std::pair<std::string, double> >::const_iterator i = keys_values.begin(); i != keys_values.end(); ++i) {
        _values[i->first] = i->second;
      }
    }

    /// Constructor from vectors of keys and values
    Weights(const std::vector<std::string>& keys, const std::vector<double>& values) {
      if (keys.size() != values.size()) {
        throw WeightError("Mismatch in lengths of keys and values vectors in Weights constructor");
      }
      for (size_t i = 0; i < keys.size(); ++i) {
        _values[keys[i]] = values[i];
      }
    }

    /// Constructor from vectors of keys and a single value, defaulting to 0.0
    Weights(const std::vector<std::string>& keys, double value=0.0) {
      for (std::vector<std::string>::const_iterator i = keys.begin(); i != keys.end(); ++i) {
        _values[*i] = value;
      }
    }

    /// @}

  public:


    /// @todo Sorted iterators of pair<std::string, double>

    /// @name Accessors
    /// @{

    typedef std::map<std::string, double>::iterator iterator;
    typedef std::map<std::string, double>::const_iterator const_iterator;
    iterator begin() { return _values.begin(); }
    const_iterator begin() const { return _values.begin(); }
    iterator end() { return _values.end(); }
    const_iterator end() const { return _values.end(); }

    double& operator [] (const std::string& key) {
      if (_values.find(key) == _values.end()) {
        throw WeightError("No weight found with supplied name");
      }
      return _values[key];
    }
    const double& operator [] (const std::string& key) const {
      const_iterator rtn = _values.find(key);
      if (rtn == _values.end()) {
        throw WeightError("No weight found with supplied name");
      }
      return rtn->second;
    }

    double& operator [] (size_t index) {
      if (index >= size()) {
        throw WeightError("Requested weight index is larger than the weights collection");
      }
      return _values[keys()[index]];
    }
    const double& operator [] (size_t index) const {
      if (index >= size()) {
        throw WeightError("Requested weight index is larger than the weights collection");
      }
      return _values.find(keys()[index])->second;
    }

    /// Number of weights keys
    unsigned int size() const {
      return _values.size();
    }

    /// Sorted list of weight keys
    std::vector<std::string> keys() const {
      std::vector<std::string> rtn;
      rtn.reserve(size());
      for (const_iterator i = begin(); i != end(); ++i) {
        rtn.push_back(i->first);
      }
      return rtn;
    }

    /// List of weight values, in the order of the sorted keys
    std::vector<double> values() const {
      std::vector<double> rtn;
      rtn.reserve(size());
      for (const_iterator i = begin(); i != end(); ++i) {
        rtn.push_back(i->second);
      }
      return rtn;
    }

    /// @}


    /// @name Arithmetic operators as members
    /// @{

    /// Add another weights to this
    Weights& operator += (const Weights& toAdd) {
      if (keys().empty()) _initToMatch(toAdd);
      if (keys() != toAdd.keys()) {
        throw WeightError("Mismatch in args to Weights += operator");
      }
      for (size_t i = 0; i < size(); ++i) {
        _values[keys()[i]] += toAdd[keys()[i]];
      }
      return *this;
    }

    /// Subtract another weights from this
    Weights& operator -= (const Weights& toSubtract) {
      if (keys().empty()) _initToMatch(toSubtract);
      if (keys() != toSubtract.keys()) {
        throw WeightError("Mismatch in args to Weights -= operator");
      }
      for (size_t i = 0; i < size(); ++i) {
        _values[keys()[i]] -= toSubtract[keys()[i]];
      }
      return *this;
    }

    /// Multiply by another weights
    Weights& operator *= (const Weights& toMultiplyBy) {
      if (keys().empty()) _initToMatch(toMultiplyBy);
      if (keys() != toMultiplyBy.keys()) {
        throw WeightError("Mismatch in args to Weights *= operator");
      }
      for (size_t i = 0; i < size(); ++i) {
        _values[keys()[i]] *= toMultiplyBy[keys()[i]];
      }
      return *this;
    }

    /// Divide by another weights
    Weights& operator /= (const Weights& toDivideBy) {
      if (keys().empty()) _initToMatch(toDivideBy);
      if (keys() != toDivideBy.keys()) {
        throw WeightError("Mismatch in args to Weights /= operator");
      }
      for (size_t i = 0; i < size(); ++i) {
        _values[keys()[i]] /= toDivideBy[keys()[i]];
      }
      return *this;
    }

    /// Multiply by a double
    Weights& operator *= (double toMultiplyBy) {
      for (size_t i = 0; i < size(); ++i) {
        _values[keys()[i]] *= toMultiplyBy;
      }
      return *this;
    }

    /// Divide by a double
    Weights& operator /= (double toDivideBy) {
      for (size_t i = 0; i < size(); ++i) {
        _values[keys()[i]] /= toDivideBy;
      }
      return *this;
    }

    /// Negate
    /// @todo Can/should this modify itself and return a reference?
    Weights operator - () const {
      Weights rtn = *this;
      rtn *= -1;
      return rtn;
    }

    /// @}


    /// @name Comparison operators
    /// @{

    /// Equals
    bool operator == (const Weights& other) const {
      return this->_values == other._values;
    }

    /// Not equals
    bool operator != (const Weights& other) const {
      return !(*this == other);
    }

    /// @}


    /// @todo Allow implicit casting to double, if single-entried? Or too dangerous and not useful enough?
    // double operator (double) () {}

  private:

    /// Initialise an empty list of weights keys to match those of another Weights object
    void _initToMatch(const Weights& other) {
      if (keys().empty()) {
        throw LogicError("Weights::_initToMatch shouldn't ever be called if there are already defined weights keys");
      }
      for (size_t i = 0; i < other.size(); ++i) {
        _values[other.keys()[i]] = 0;
      }
    }


  private:

    std::map<std::string, double> _values;

  };


  /// @name Combining weights: global operators
  /// @{

  /// Add two weights
  inline Weights operator + (const Weights& first, const Weights& second) {
    Weights tmp = first;
    tmp += second;
    return tmp;
  }

  /// Subtract two weights
  inline Weights operator - (const Weights& first, const Weights& second) {
    Weights tmp = first;
    tmp -= second;
    return tmp;
  }

  /// Multiply two weights
  inline Weights operator * (const Weights& first, const Weights& second) {
    Weights tmp = first;
    tmp *= second;
    return tmp;
  }

  /// Divide two weights
  inline Weights operator / (const Weights& numer, const Weights& denom) {
    Weights tmp = numer;
    tmp /= denom;
    return tmp;
  }


  /// Multiply by a double
  inline Weights operator * (double a, const Weights& w) {
    Weights tmp = w;
    tmp *= a;
    return tmp;
  }

  /// Multiply by a double
  inline Weights operator * (const Weights& w, double a) {
    return a * w;
  }

  /// Divide by a double
  inline Weights operator / (const Weights& w, double a) {
    Weights tmp = w;
    tmp /= a;
    return tmp;
  }

  /// Divide a double by a Weights
  /// @todo Is this really needed?
  inline Weights operator / (double a, const Weights& w) {
    Weights tmp(w.keys(), a);
    tmp /= w;
    return tmp;
  }

  /// @}


  /// Standard text representaion
  inline std::ostream& operator<<(std::ostream& out, const Weights& w) {
    out << "{ ";
    for (Weights::const_iterator i = w.begin(); i != w.end(); ++i) {
      if (i != w.begin()) out << ", ";
      out << i->first << ": " << i->second;
    }
    out << "}";
    return out;
  }


}

#endif
