// -*- C++ -*-
//
// This file is part of YODA -- Yet more Objects for Data Analysis
// Copyright (C) 2008-2021 The YODA collaboration (see AUTHORS for details)
//
#ifndef YODA_Exception_h
#define YODA_Exception_h

#include <string>
#include <exception>
#include <stdexcept>

namespace YODA {


  /// @brief Generic unspecialised YODA runtime error.
  ///
  /// NB. We don't use "Error" because that's a useful stats word to have
  /// available!
  class Exception : public std::runtime_error {
  public:
    Exception(const std::string& what);
  };


  /// Error for general binning problems.
  class BinningError : public Exception {
  public:
    BinningError(const std::string& what);
  };


  /// Error for e.g. use of invalid bin ranges.
  class RangeError : public Exception {
  public:
    RangeError(const std::string& what);
  };


  /// Error for modification of a data object where filling has already begun.
  class LockError : public Exception {
  public:
    LockError(const std::string& what);
  };


  /// Error to throw when a slicing is requested on a non-slicable state of an object.
  class GridError : public Exception {
  public:
    GridError(const std::string& what);
  };


  /// Error for places where it should not have been possible to get to!
  class LogicError : public Exception {
  public:
    LogicError(const std::string& what);
  };


  /// @brief Errors relating to event/bin weights
  ///
  /// Arises in computing statistical quantities because e.g. the bin
  /// weight is zero or negative.
  class WeightError : public Exception {
  public:
    WeightError(const std::string& what);
  };


  /// Errors relating to insufficient (effective) statistics
  class LowStatsError : public Exception {
  public:
    LowStatsError(const std::string& what);
  };


  /// Error for unfound or broken AnalysisObject annotations
  class AnnotationError : public Exception {
  public:
    AnnotationError(const std::string& what);
  };


  /// Error for file reading errors
  class ReadError : public Exception {
  public:
    ReadError(const std::string& what);
  };


  /// Error for file writing errors
  class WriteError : public Exception {
  public:
    WriteError(const std::string& what);
  };


  /// Error for problems introduced outside YODA, to put it nicely.
  class UserError : public Exception {
  public:
    UserError(const std::string& what);
  };


}

#endif
