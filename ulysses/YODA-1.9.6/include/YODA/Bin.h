// -*- C++ -*-
//
// This file is part of YODA -- Yet more Objects for Data Analysis
// Copyright (C) 2008-2021 The YODA collaboration (see AUTHORS for details)
//
#ifndef YODA_Bin_h
#define YODA_Bin_h

#include <string>
#include <utility>

namespace YODA {


  /// @brief Base class for bins in 1D and 2D histograms.
  ///
  /// This base class only provides very basic functionality for fill
  /// weight statistics access, as 1D/2D and basic/profile histos have
  /// quite difference implementations.
  class Bin {
  public:

    /// Virtual destructor for inheritance
    virtual ~Bin() { }

    /// @name Miscellaneous
    /// @{

    /// Reset this bin
    virtual void reset() = 0;

    /// @}


    /// @name Dimensions
    /// @{

    /// Dimension of the fill space
    ///
    /// @todo Convert to be the total dimension
    virtual size_t dim() = 0;

    /// Dimension of the fill space
    virtual size_t fillDim() = 0;

    /// @}


    /// @name Fill statistics
    /// @{

    /// The number of entries (fractional fills are possible)
    virtual double numEntries() const = 0;

    /// The effective number of entries
    virtual double effNumEntries() const = 0;

    /// The sum of weights
    virtual double sumW() const = 0;

    /// The sum of weights squared
    virtual double sumW2() const = 0;

    /// @}


    /// @todo Add integer ID access to axis quantities (i.e. min, max, mid, focus)


  };


}



#endif
