// -*- C++ -*-
//
// This file is part of YODA -- Yet more Objects for Data Analysis
// Copyright (C) 2008-2017 The YODA collaboration (see AUTHORS for details)
//
#ifndef YODA_TRAITS_H
#define YODA_TRAITS_H

#include <type_traits>

namespace YODA {


  namespace SFINAE {
    /// @todo Replace by C++11 std::false/true_type
    typedef char yes[1]; typedef char no[2];

    /// C++11 equivalent of C++17 std::void_t
    template <typename ...>
    using void_t = void;
  }


  /// SFINAE definition of dereferenceability trait, cf. Boost has_dereference
  template <typename T, typename=void>
  struct Derefable : std::false_type {};
  //
  template <typename T>
  struct Derefable<T, SFINAE::void_t< decltype(*std::declval<T>())> > : std::true_type {};


  /// SFINAE struct to check for dereferencing to an AnalysisObject (reference) at compile time
  template <typename T, typename=AnalysisObject>
  struct DerefableToAO : std::false_type {};
  //
  template <typename T>
  struct DerefableToAO<T, typename std::conditional<std::is_base_of<AnalysisObject, typename std::decay< decltype(*std::declval<T>()) >::type>::value,
                                                    AnalysisObject, void>::type> : std::true_type {};
  // The following would not be enough since it doesn't work with e.g. Histo1D*:
  // struct DerefableToAO<T, typename std::decay< decltype(*std::declval<T>()) >::type > : std::true_type {};


  // SFINAE struct to check for iterator concept at compile time
  /// @todo Replace with C++11 std stuff
  template<typename T>
  struct Iterable {
    template <typename C> static SFINAE::yes& test(typename C::iterator* c);
    template <typename> static SFINAE::no& test(...);
    static const bool value = sizeof(test<T>(0)) == sizeof(SFINAE::yes);
  };


  // SFINAE struct to check for const_iterator concept at compile time
  /// @todo Replace with C++11 std stuff
  template<typename T>
  struct CIterable {
    template <typename C> static SFINAE::yes& test(typename C::const_iterator* c);
    template <typename> static SFINAE::no& test(...);
    static const bool value = sizeof(test<T>(0)) == sizeof(SFINAE::yes);
  };


  // SFINAE struct to check for pushing concept at compile time
  /// @todo Replace with C++11 std stuff
  template<typename T,typename VAL>
  struct Pushable {
    template <typename SIG,SIG> struct has_sig;
    template<typename C> static SFINAE::yes& test(has_sig<void (C::*) (const VAL&),&C::push_back>*);
    template<typename C> static SFINAE::no& test(...);
    static const bool value = sizeof(test<T>(0)) == sizeof(SFINAE::yes);
  };


}

#endif
