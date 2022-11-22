// -*- C++ -*-
//
// This file is part of YODA -- Yet more Objects for Data Analysis
// Copyright (C) 2008-2021 The YODA collaboration (see AUTHORS for details)
//
#ifndef YODA_Profile1D_h
#define YODA_Profile1D_h

#include "YODA/AnalysisObject.h"
#include "YODA/Fillable.h"
#include "YODA/Binned.h"
#include "YODA/ProfileBin1D.h"
#include "YODA/Scatter2D.h"
#include "YODA/Dbn2D.h"
#include "YODA/Axis1D.h"
#include "YODA/Exceptions.h"
#include <vector>
#include <string>
#include <map>
#include <tuple>

namespace YODA {


  // Forward declarations
  class Histo1D;
  class Scatter2D;


  /// Convenience typedef
  typedef Axis1D<ProfileBin1D, Dbn2D> Profile1DAxis;


  /// A one-dimensional profile histogram
  class Profile1D : public AnalysisObject, public Fillable, public Binned {
  public:

    /// Convenience typedefs
    typedef Profile1DAxis Axis;
    typedef Axis::Bins Bins;
    typedef ProfileBin1D Bin;

    typedef std::tuple<double, double> FillType;
    typedef double BinType;
    typedef std::shared_ptr<Profile1D> Ptr;


    /// @name Constructors
    /// @{

    /// Default constructor
    Profile1D(const std::string& path="", const std::string& title="")
      : AnalysisObject("Profile1D", path, title),
        _axis()
    { }


    /// Constructor giving range and number of bins
    Profile1D(size_t nxbins, double xlower, double xupper,
              const std::string& path="", const std::string& title="")
      : AnalysisObject("Profile1D", path, title),
        _axis(nxbins, xlower, xupper)
    { }


    /// Constructor giving explicit bin edges
    ///
    /// For n bins, binedges.size() == n+1, the last one being the upper bound
    /// of the last bin
    Profile1D(const std::vector<double>& xbinedges,
              const std::string& path="", const std::string& title="")
      : AnalysisObject("Profile1D", path, title),
        _axis(xbinedges)
    { }


    /// Copy constructor with optional new path
    /// @todo Also allow title setting from the constructor?
    Profile1D(const Profile1D& p, const std::string& path="");


    /// Constructor from a Scatter2D's binning, with optional new path
    /// @todo Also allow title setting from the constructor?
    Profile1D(const Scatter2D& s, const std::string& path="");


    /// Constructor from a Histo1D's binning, with optional new path
    /// @todo Also allow title setting from the constructor?
    Profile1D(const Histo1D& h, const std::string& path="");


    /// @brief State-setting constructor.
    ///
    /// Intended principally for internal persistency use.
    Profile1D(const std::vector<ProfileBin1D>& bins,
              const Dbn2D& dbn_tot, const Dbn2D& dbn_uflow, const Dbn2D& dbn_oflow,
              const std::string& path="", const std::string& title="")
      : AnalysisObject("Profile1D", path, title),
        _axis(bins, dbn_tot, dbn_uflow, dbn_oflow)
    { }


    /// Assignment operator
    Profile1D& operator = (const Profile1D& p1) {
      AnalysisObject::operator = (p1); //< AO treatment of paths etc.
      _axis = p1._axis;
      return *this;
    }

    /// Make a copy on the stack
    Profile1D clone() const {
      return Profile1D(*this);
    }

    /// Make a copy on the heap, via 'new'
    Profile1D* newclone() const {
      return new Profile1D(*this);
    }

    /// @}


    /// @brief Fill dimension of this data object
    ///
    /// @todo Change this to the total dimension (in v2)
    size_t dim() const { return 1; }

    /// Fill dimension of this data object
    size_t fillDim() const { return 1; }


    /// @name Modifiers
    /// @{

    /// Fill histo by value and weight
    virtual void fill(double x, double y, double weight=1.0, double fraction=1.0);
    void fill(const FillType & xs, double weight=1.0, double fraction=1.0) {
        fill(std::get<0>(xs), std::get<1>(xs), weight, fraction);
    }

    /// Fill histo x bin i with the given y value and weight
    virtual void fillBin(size_t i, double y, double weight=1.0, double fraction=1.0);


    /// @brief Reset the histogram
    ///
    /// Keep the binning but set all bin contents and related quantities to zero
    void reset() {
      _axis.reset();
    }


    /// Rescale as if all fill weights had been different by factor @a scalefactor.
    void scaleW(double scalefactor) {
      _axis.scaleW(scalefactor);
    }


    /// Rescale as if all y values had been different by factor @a scalefactor.
    void scaleY(double scalefactor) {
      _axis.totalDbn().scaleY(scalefactor);
      _axis.overflow().scaleY(scalefactor);
      _axis.underflow().scaleY(scalefactor);
      for (size_t i = 0; i < bins().size(); ++i)
        bin(i).scaleY(scalefactor);
    }


    /// Merge together the bin range with indices from @a from to @a to, inclusive
    void mergeBins(size_t from, size_t to) {
      _axis.mergeBins(from, to);
    }

    /// check if binning is the same as different Profile1D
    bool sameBinning(const Profile1D& p1) {
      return _axis == p1._axis;
    }

    /// Merge every group of n bins, starting from the LHS
    void rebinBy(unsigned int n, size_t begin=0, size_t end=UINT_MAX) {
      _axis.rebinBy(n, begin, end);
    }
    /// Overloaded alias for rebinBy
    void rebin(unsigned int n, size_t begin=0, size_t end=UINT_MAX) {
      rebinBy(n, begin, end);
    }

    /// Rebin to the given list of bin edges
    void rebinTo(const std::vector<double>& newedges) {
      _axis.rebinTo(newedges);
    }
    /// Overloaded alias for rebinTo
    void rebin(const std::vector<double>& newedges) {
      rebinTo(newedges);
    }

    /// @}


    /// @name Bin adding and removing
    /// @{

    /// Bin addition operator
    void addBin(double xlow, double xhigh) {
      _axis.addBin(xlow, xhigh);
    }

    /// Bin addition operator
    void addBins(const std::vector<double> binedges) {
      _axis.addBins(binedges);
    }

    // /// Bin addition operator
    // void addBins(const std::vector<std::pair<double,double> > edges) {
    //   _axis.addBins(edges);
    // }

    /// Add a new bin, perhaps already populated: CAREFUL!
    void addBin(const ProfileBin1D& b) { _axis.addBin(b); }

    /// @brief Bins addition operator
    ///
    /// Add multiple bins without resetting
    void addBins(const Bins& bins) {
      _axis.addBins(bins);
    }

    void rmBin(size_t index) {
      _axis.eraseBin(index);
    }

    /// @}


    /// @name Bin accessors
    /// @{

    /// Number of bins (not counting under/overflow)
    size_t numBins() const { return bins().size(); }

    /// Number of bins on the x-axis (not counting under/overflow)
    size_t numBinsX() const { return numBins(); }

    /// Low edge of this histo's axis
    double xMin() const { return _axis.xMin(); }

    /// High edge of this histo's axis
    double xMax() const { return _axis.xMax(); }

    /// All bin edges on this histo's axis
    ///
    /// @note This only returns the finite edges, i.e. -inf and +inf are removed
    /// @todo Make the +-inf stripping controllable by a default-valued bool arg
    std::vector<double> xEdges() const { return _axis.xEdges(); }

    /// All bin widths on this histo's axis
    ///
    /// @note This only returns the finite edges, i.e. -inf and +inf are removed
    /// @todo Make the +-inf stripping controllable by a default-valued bool arg
    std::vector<double> xWidths() const { return _axis.xWidths(); }

    /// Access the bin vector
    std::vector<YODA::ProfileBin1D>& bins() { return _axis.bins(); }

    /// Access the bin vector
    const std::vector<YODA::ProfileBin1D>& bins() const { return _axis.bins(); }


    /// Access a bin by index (non-const version)
    ProfileBin1D& bin(size_t index) { return _axis.bins()[index]; }
    /// Access a bin by index (const version)
    const ProfileBin1D& bin(size_t index) const { return _axis.bins()[index]; }

    /// Access a bin index by x-coordinate.
    int binIndexAt(double x) { return _axis.binIndexAt(x); }

    /// Access a bin by x-coordinate (const version)
    const ProfileBin1D& binAt(double x) const { return _axis.binAt(x); }


    /// Access summary distribution, including gaps and overflows (non-const version)
    Dbn2D& totalDbn() { return _axis.totalDbn(); }
    /// Access summary distribution, including gaps and overflows (const version)
    const Dbn2D& totalDbn() const { return _axis.totalDbn(); }
    /// Set summary distribution, mainly for persistency: CAREFUL!
    void setTotalDbn(const Dbn2D& dbn) { _axis.setTotalDbn(dbn); }


    /// Access underflow (non-const version)
    Dbn2D& underflow() { return _axis.underflow(); }
    /// Access underflow (const version)
    const Dbn2D& underflow() const { return _axis.underflow(); }
    /// Set underflow distribution, mainly for persistency: CAREFUL!
    void setUnderflow(const Dbn2D& dbn) { _axis.setUnderflow(dbn); }


    /// Access overflow (non-const version)
    Dbn2D& overflow() { return _axis.overflow(); }
    /// Access overflow (const version)
    const Dbn2D& overflow() const { return _axis.overflow(); }
    /// Set overflow distribution, mainly for persistency: CAREFUL!
    void setOverflow(const Dbn2D& dbn) { _axis.setOverflow(dbn); }

    /// @}


    /// @name Whole histo data
    /// @{

    /// @todo Add integrals? Or are they too ambiguous to make a core function?

    /// Get the number of fills (fractional fills are possible)
    double numEntries(bool includeoverflows=true) const;

    /// Get the effective number of fills
    double effNumEntries(bool includeoverflows=true) const;

    /// Get sum of weights in histo.
    double sumW(bool includeoverflows=true) const;

    /// Get sum of squared weights in histo.
    double sumW2(bool includeoverflows=true) const;

    /// Get the mean x
    double xMean(bool includeoverflows=true) const;

    /// Get the variance in x
    double xVariance(bool includeoverflows=true) const;

    /// Get the standard deviation in x
    double xStdDev(bool includeoverflows=true) const {
      return std::sqrt(xVariance(includeoverflows));
    }

    /// Get the standard error on the mean x
    double xStdErr(bool includeoverflows=true) const;

    /// Get the RMS in x
    double xRMS(bool includeoverflows=true) const;

    /// @}


    /// @name Adding and subtracting histograms
    /// @{

    /// Add another profile to this one
    Profile1D& operator += (const Profile1D& toAdd) {
      if (hasAnnotation("ScaledBy")) rmAnnotation("ScaledBy");
      _axis += toAdd._axis;
      return *this;
    }

    /// Subtract another profile from this one
    Profile1D& operator -= (const Profile1D& toSubtract) {
      if (hasAnnotation("ScaledBy")) rmAnnotation("ScaledBy");
      _axis -= toSubtract._axis;
      return *this;
    }

    inline bool operator == (const Profile1D& other){
      return _axis == other._axis;
    }

    inline bool operator != (const Profile1D& other){
      return ! operator == (other);
    }
    /// @}


  protected:

    /// Access a bin by x-coordinate (non-const version)
    ProfileBin1D& _binAt(double x) { return _axis.binAt(x); }


  private:

    /// @name Bin data
    /// @{

    /// The bins contained in this profile histogram
    Axis1D<ProfileBin1D, Dbn2D> _axis;

    /// @}

  };


  /// Convenience typedef
  typedef Profile1D P1D;


  /// @name Combining profile histos: global operators
  /// @{


  /// Add two profile histograms
  inline Profile1D add(const Profile1D& first, const Profile1D& second) {
    Profile1D tmp = first;
    if (first.path() != second.path()) tmp.setPath("");
    tmp += second;
    return tmp;
  }


  /// Add two profile histograms
  inline Profile1D operator + (const Profile1D& first, const Profile1D& second) {
    return add(first, second);
  }


  /// Subtract two profile histograms
  inline Profile1D subtract(const Profile1D& first, const Profile1D& second) {
    Profile1D tmp = first;
    if (first.path() != second.path()) tmp.setPath("");
    tmp -= second;
    return tmp;
  }


  /// Subtract two profile histograms
  inline Profile1D operator - (const Profile1D& first, const Profile1D& second) {
    return subtract(first, second);
  }


  /// Divide two profile histograms
  Scatter2D divide(const Profile1D& numer, const Profile1D& denom);


  /// Divide two profile histograms
  inline Scatter2D operator / (const Profile1D& numer, const Profile1D& denom) {
    return divide(numer, denom);
  }

  /// @}


}

#endif
