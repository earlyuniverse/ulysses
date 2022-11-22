// -*- C++ -*-
//
// This file is part of YODA -- Yet more Objects for Data Analysis
// Copyright (C) 2008-2020 The YODA collaboration (see AUTHORS for details)
//
#ifndef YODA_Binned_h
#define YODA_Binned_h

#include "YODA/AnalysisObject.h"
#include "YODA/Exceptions.h"

namespace YODA {


  /// A base class for all binned objects
  class Binned {
  public:

    /// @name Constructors
    //@{

    /// Virtual destructor for inheritance
    virtual ~Binned() = default;


    /// @name Modifiers
    //@{

    /// Reset
    virtual void reset() = 0;

    // /// Fill histo bin @a i with the given weight, optionally as a fractional fill
    // virtual void fillBin(size_t i, double weight=1.0, double fraction=1.0) = 0;


    /// @todo These make sense in 1D; less clear about other dimensionalities
    ///
    // /// Merge every group of n bins, starting from the LHS
    // void rebinBy(unsigned int n, size_t begin=0, size_t end=UINT_MAX) {
    //   _axis.rebinBy(n, begin, end);
    // }
    // /// Overloaded alias for rebinBy
    // void rebin(unsigned int n, size_t begin=0, size_t end=UINT_MAX) {
    //   rebinBy(n, begin, end);
    // }

    // /// Rebin to the given list of bin edges
    // void rebinTo(const std::vector<double>& newedges) {
    //   _axis.rebinTo(newedges);
    // }
    // /// Overloaded alias for rebinTo
    // void rebin(const std::vector<double>& newedges) {
    //   rebinTo(newedges);
    // }


    /// @name Bin accessors
    /// @{

    /// Number of bins (not counting under/overflow)
    virtual size_t numBins() const = 0;

    /// Check if binning is the same as different Histo1D
    ///
    /// @todo Needs RTTI for inequivalent types
    // bool sameBinning(const Binned& bd) {
    //   return _axis == h1._axis;
    // }

    /// @}


    /// @name Bin adding and removal
    /// @{

    // /// Add a new bin specifying its lower and upper bound
    // void addBin(double from, double to) { _axis.addBin(from, to); }

    // /// Add new bins by specifying a vector of edges
    // void addBins(std::vector<double> edges) { _axis.addBins(edges); }

    // // /// Add new bins specifying a beginning and end of each of them
    // // void addBins(std::vector<std::pair<double,double> > edges) {
    // //   _axis.addBins(edges);
    // // }

    // /// Add a new bin, perhaps already populated: CAREFUL!
    // void addBin(const HistoBin1D& b) { _axis.addBin(b); }

    // /// @brief Bins addition operator
    // ///
    // /// Add multiple bins without resetting
    // void addBins(const Bins& bins) {
    //   _axis.addBins(bins);
    // }

    /// Remove a bin
    virtual void rmBin(size_t index) = 0;

    /// Remove several bins
    ///
    /// @todo Implement, cf. Scatter::rmPoints()
    // virtual void rmBins(size_t index) = 0;

    /// @deprecated Use rmBin
    void eraseBin(size_t index) { rmBin(index); }

    //@}

  };


}

#endif
