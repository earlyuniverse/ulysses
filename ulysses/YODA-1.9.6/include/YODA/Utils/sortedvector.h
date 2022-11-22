// -*- C++ -*-
//
// This file is part of YODA -- Yet more Objects for Data Analysis
// Copyright (C) 2008-2017 The YODA collaboration (see AUTHORS for details)
//
#ifndef YODA_SORTEDVECTOR_H
#define YODA_SORTEDVECTOR_H

#include <vector>
#include <algorithm>
#include <stdexcept>
#include <iostream>

namespace YODA {
  namespace Utils {


    /// @brief Specialisation of std::vector to allow indexed access to ordered elements
    ///
    /// @warning This still scales as n^2 log n or so, if inserting n
    /// elements. Prefer to populate a whole std::vector first, then add, so the
    /// sorting/searching only needs to be done once.
    ///
    /// @todo Need to template on the value-comparison definition?
    /// @todo Generalise the type of source container for constructor argument
    template <typename T>
    class sortedvector : public std::vector<T> {
    public:

      /// Default constructor
      sortedvector() {}

      /// Conversion from std::vector
      sortedvector(const std::vector<T> & vec)
        : std::vector<T>(vec) {
        std::sort(this->begin(), this->end());
      }

      /// Insertion operator (push_back should not be used!)
      void insert(const T& val) {
        // Dumb way:
        //   std::vector<T>::push_back(val);
        //   std::sort(this->begin(), this->end());
        //
        // A bit better:
        std::vector<T>::insert(std::upper_bound(std::vector<T>::begin(), std::vector<T>::end(), val), val);
      }


    private:

      /// Hiding push_back from the base class
      void push_back();

    };


  }
}

#endif
