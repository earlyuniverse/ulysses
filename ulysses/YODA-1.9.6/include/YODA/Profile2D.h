// -*- C++ -*-
//
// This file is part of YODA -- Yet more Objects for Data Analysis
// Copyright (C) 2008-2021 The YODA collaboration (see AUTHORS for details)
//
#ifndef YODA_Profile2D_h
#define YODA_Profile2D_h

#include "YODA/AnalysisObject.h"
#include "YODA/Fillable.h"
#include "YODA/Binned.h"
#include "YODA/ProfileBin2D.h"
#include "YODA/Dbn3D.h"
#include "YODA/Axis2D.h"
#include "YODA/Scatter3D.h"
#include "YODA/Exceptions.h"
#include <vector>
#include <tuple>

namespace YODA {


  // Forward declarations
  class Histo2D;
  class Scatter3D;

  /// Convenience typedef
  typedef Axis2D<ProfileBin2D, Dbn3D> Profile2DAxis;


  /// A two-dimensional profile histogram
  class Profile2D : public AnalysisObject, public Fillable, public Binned {
  public:

    /// Convenience typedefs
    typedef Profile2DAxis Axis;
    typedef Axis::Bins Bins;
    typedef ProfileBin2D Bin;
    typedef Axis::Outflows Outflows;

    typedef std::tuple<double, double, double> FillType;
    typedef std::tuple<double, double> BinType;
    typedef std::shared_ptr<Profile2D> Ptr;


    /// @name Constructors
    /// @{

    /// Default constructor
    Profile2D(const std::string& path="", const std::string& title="")
      : AnalysisObject("Profile2D", path, title),
        _axis()
    { }


    /// Constructor giving range and number of bins
    Profile2D(size_t nbinsX, double lowerX, double upperX,
              size_t nbinsY, double lowerY, double upperY,
              const std::string& path="", const std::string& title="")
      : AnalysisObject("Profile2D", path, title),
        _axis(nbinsX, std::make_pair(lowerX, upperX),  nbinsY, std::make_pair(lowerY, upperY))
    { }


    /// Constructor giving explicit bin edges in the direction of X and Y
    Profile2D(const std::vector<double>& xedges, const std::vector<double>& yedges,
              const std::string& path="", const std::string& title="")
      : AnalysisObject("Profile2D", path, title),
        _axis(xedges, yedges)
    { }


    /// Constructor accepting an explicit collection of bins.
    Profile2D(const std::vector<Bin>& bins,
              const std::string& path="", const std::string& title="")
      : AnalysisObject("Profile2D", path, title),
        _axis(bins)
    { }


    /// A copy constructor with optional new path
    /// @todo Also allow title setting from the constructor?
    Profile2D(const Profile2D& p, const std::string& path="");

    /// A constructor from a Scatter3D's binning, with optional new path
    /// @todo Also allow title setting from the constructor?
    Profile2D(const Scatter3D& s, const std::string& path="");

    /// Constructor from a Histo2D's binning, with optional new path
    /// @todo Also allow title setting from the constructor?
    Profile2D(const Histo2D& h, const std::string& path="");

    /// @brief State-setting constructor
    ///
    /// Mainly intended for internal persistency use.
    Profile2D(const std::vector<ProfileBin2D>& bins,
              const Dbn3D& totalDbn,
              const Outflows& outflows,
              const std::string& path="", const std::string& title="")
      : AnalysisObject("Profile2D", path, title),
        _axis(bins, totalDbn, outflows)
    { }


    /// Assignment operator
    Profile2D& operator = (const Profile2D& p2) {
      AnalysisObject::operator = (p2); //< AO treatment of paths etc.
      _axis = p2._axis;
      return *this;
    }


    /// Make a copy on the stack
    Profile2D clone() const {
      return Profile2D(*this);
    }

    /// Make a copy on the heap, via 'new'
    Profile2D* newclone() const {
      return new Profile2D(*this);
    }

    /// @}


    /// @brief Fill dimension of this data object
    ///
    /// @todo Change this to the total dimension (in v2)
    size_t dim() const { return 2; }

    /// Fill dimension of this data object
    size_t fillDim() const { return 2; }


    /// @name Modifiers
    /// @{

    /// Fill histo by value and weight
    virtual void fill(double x, double y, double z, double weight=1.0, double fraction=1.0);
    virtual void fill(const FillType & xs, double weight=1.0, double fraction=1.0) {
        fill(std::get<0>(xs), std::get<1>(xs), std::get<2>(xs), weight, fraction);
    }

    /// Fill histo x-y bin i with the given z value and weight
    virtual void fillBin(size_t i, double z, double weight=1.0, double fraction=1.0);


    /// @brief Reset the histogram
    ///
    /// Keep the binning but reset the statistics
    void reset() {
      _axis.reset();
    }

    /// Rescale as if all fill weights had been different by a @a scalefactor
    void scaleW(double scalefactor) {
      /// @todo Is this ScaledBy annotation needed?
      setAnnotation("ScaledBy", annotation<double>("ScaledBy", 1.0) * scalefactor);
      _axis.scaleW(scalefactor);
    }

    /// Rescale as if all z values had been different by factor @a scalefactor.
    void scaleZ(double scalefactor) {
      _axis.totalDbn().scaleZ(scalefactor);
      /// @todo Need to rescale overflows too, when they exist.
      // _axis.overflow().scaleZ(scalefactor);
      // _axis.underflow().scaleZ(scalefactor);
      for (size_t i = 0; i < bins().size(); ++i)
        bin(i).scaleZ(scalefactor);
    }


    /// @todo TODO
    // /// Merge together the bin range with indices from @a from to @a to, inclusive
    // void mergeBins(size_t from, size_t to) {
    //   _axis.mergeBins(from, to);
    // }

    /// @todo TODO
    // /// Merge every group of n bins, starting from the LHS
    // void rebin(size_t n) {
    //   throw "IMPLEMENT!";
    //   //_axis.rebin(n);
    // }

    /// @}


    /// @name Bin adding and removing
    /// @{

    // /// @brief Bin addition operator
    // ///
    // /// Add a bin to the axis, described by its x and y ranges.
    void addBin(Axis::EdgePair1D xrange, Axis::EdgePair1D yrange) {
       _axis.addBin(xrange, yrange);
    }

    // /// @brief Bin addition operator
    // ///
    // /// Add a bin to the axis, possibly pre-populated
    void addBin(const Bin& bin) {
       _axis.addBin(bin);
    }

    /// @brief Bins addition operator
    ///
    /// Add multiple bins from edge cuts without resetting
    void addBins(const Axis::Edges& xcuts, const Axis::Edges& ycuts) {
      _axis.addBins(xcuts, ycuts);
    }


    /// @brief Bins addition operator
    ///
    /// Add multiple bins without resetting
    void addBins(const Bins& bins) {
      _axis.addBins(bins);
    }

    /// Check if binning is the same as different Profile2D
    bool sameBinning(const Profile2D& p2) {
      return _axis == p2._axis;
    }


    /// @todo TODO
    // /// @brief Bin addition operator
    // ///
    // /// Add a set of bins delimiting coordinates of which are contained
    // /// in binLimits vector.
    // void addBin(const std::vector<Segment>& binLimits) {
    //   _axis.addBin(binLimits);
    // }

    void rmBin(size_t index) {
      _axis.eraseBin(index);
    }

    /// @}


    /// @name Bin accessors
    /// @{

    /// All bin edges on this histo's x axis
    ///
    /// @note This only returns the finite edges, i.e. -inf and +inf are removed
    /// @todo Make the +-inf stripping controllable by a default-valued bool arg
    std::vector<double> xEdges() const { return _axis.xEdges(); }

    /// All bin edges on this histo's y axis
    ///
    /// @note This only returns the finite edges, i.e. -inf and +inf are removed
    /// @todo Make the +-inf stripping controllable by a default-valued bool arg
    std::vector<double> yEdges() const { return _axis.yEdges(); }

    /// All bin widths on this histo's x axis
    ///
    /// @note This only returns the finite edges, i.e. -inf and +inf are removed
    /// @todo Make the +-inf stripping controllable by a default-valued bool arg
    std::vector<double> xWidths() const { return _axis.xWidths(); }

    /// All bin widths on this histo's y axis
    ///
    /// @note This only returns the finite edges, i.e. -inf and +inf are removed
    /// @todo Make the +-inf stripping controllable by a default-valued bool arg
    std::vector<double> yWidths() const { return _axis.yWidths(); }

    /// @todo Add xMins, xMaxs, xMids, xFoci, and y-versions

    /// Low x edge of this histo's axis
    double xMin() const { return _axis.xMin(); }

    /// High x edge of this histo's axis
    double xMax() const { return _axis.xMax(); }


    /// Low y edge of this histo's axis
    double yMin() const { return _axis.yMin(); }

    /// High y edge of this histo's axis
    double yMax() const { return _axis.yMax(); }


    /// Access the bin vector (non-const)
    std::vector<YODA::ProfileBin2D>& bins() { return _axis.bins(); }

    /// Access the bin vector (const)
    const std::vector<YODA::ProfileBin2D>& bins() const { return _axis.bins(); }


    /// Access a bin by index (non-const)
    ProfileBin2D& bin(size_t index) { return _axis.bins()[index]; }

    /// Access a bin by index (const)
    const ProfileBin2D& bin(size_t index) const { return _axis.bins()[index]; }

    /// Access a bin index by coordinate
    int binIndexAt(double x, double y) { return _axis.binIndexAt(x, y); }
    int binIndexAt(const BinType& t) { return _axis.binIndexAt(std::get<0>(t), std::get<1>(t)); }

    /// Access a bin by coordinate (const)
    const ProfileBin2D& binAt(double x, double y) const { return _axis.binAt(x, y); }

    const ProfileBin2D& binAt(const BinType& t) const { return _axis.binAt(std::get<0>(t), std::get<1>(t)); }


    /// Number of bins of this axis (not counting under/over flow)
    size_t numBins() const { return _axis.bins().size(); }

    /// Number of bins along the x axis
    size_t numBinsX() const { return _axis.numBinsX(); }

    /// Number of bins along the y axis
    size_t numBinsY() const { return _axis.numBinsY(); }


    /// Access summary distribution, including gaps and overflows (non-const version)
    Dbn3D& totalDbn() { return _axis.totalDbn(); }

    /// Access summary distribution, including gaps and overflows (const version)
    const Dbn3D& totalDbn() const { return _axis.totalDbn(); }

    /// Set summary distribution, including gaps and overflows
    void setTotalDbn(const Dbn3D& dbn) { _axis.setTotalDbn(dbn); }


    // /// @brief Access an outflow (non-const)
    // ///
    // /// Two indices are used, for x and y: -1 = underflow, 0 = in-range, and +1 = overflow.
    // /// (0,0) is not a valid overflow index pair, since it is in range for both x and y.
    // Dbn3D& outflow(int ix, int iy) {
    //   return _axis.outflow(ix, iy);
    // }

    // /// @brief Access an outflow (const)
    // ///
    // /// Two indices are used, for x and y: -1 = underflow, 0 = in-range, and +1 = overflow.
    // /// (0,0) is not a valid overflow index pair, since it is in range for both x and y.
    // const Dbn3D& outflow(int ix, int iy) const {
    //   return _axis.outflow(ix, iy);
    // }

    /// @}


    /// @name Whole histo data
    /// @{

    /// Get the number of fills (fractional fills are possible)
    double numEntries(bool includeoverflows=true) const;

    /// Get the effective number of fills
    double effNumEntries(bool includeoverflows=true) const;

    /// Get sum of weights in histo
    double sumW(bool includeoverflows=true) const;

    /// Get the sum of squared weights in histo
    double sumW2(bool includeoverflows=true) const;

    /// Get the mean x
    double xMean(bool includeoverflows=true) const;

    /// Get the mean y
    double yMean(bool includeoverflows=true) const;

    /// Get the variance in x
    double xVariance(bool includeoverflows=true) const;

    /// Get the variance in y
    double yVariance(bool includeoverflows=true) const;

    /// Get the standard deviation in x
    double xStdDev(bool includeoverflows=true) const {
      return std::sqrt(xVariance(includeoverflows));
    }

    /// Get the standard deviation in y
    double yStdDev(bool includeoverflows=true) const {
      return std::sqrt(yVariance(includeoverflows));
    }

    /// Get the standard error on the mean x
    double xStdErr(bool includeoverflows=true) const;

    /// Get the standard error on the mean y
    double yStdErr(bool includeoverflows=true) const;

    /// Get the RMS in x
    double xRMS(bool includeoverflows=true) const;

    /// Get the RMS in y
    double yRMS(bool includeoverflows=true) const;

    /// @}


    /// @name Adding and subtracting histograms
    /// @{

    /// Add another profile to this one
    Profile2D& operator += (const Profile2D& toAdd) {
      if (hasAnnotation("ScaledBy")) rmAnnotation("ScaledBy");
      _axis += toAdd._axis;
      return *this;
    }

    /// Subtract another profile from this one
    Profile2D& operator -= (const Profile2D& toSubtract) {
      if (hasAnnotation("ScaledBy")) rmAnnotation("ScaledBy");
      _axis -= toSubtract._axis;
      return *this;
    }

    inline bool operator == (const Profile2D& other){
      return _axis == other._axis;
    }

    inline bool operator != (const Profile2D& other){
      return ! operator == (other);
    }
    /// @}-


  protected:

    /// Access a bin by coordinate (non-const)
    ProfileBin2D& _binAt(double x, double y) { return _axis.binAt(x, y); }


  private:

    /// @name Bin data
    /// @{

    /// The bins contained in this profile histogram
    Axis2D<ProfileBin2D, Dbn3D> _axis;

    /// @}
  };


  /// Convenience typedef
  typedef Profile2D P2D;


  /// @name Combining profile histos: global operators
  /// @{

  /// Add two profile histograms
  inline Profile2D add(const Profile2D& first, const Profile2D& second) {
    Profile2D tmp = first;
    if (first.path() != second.path()) tmp.setPath("");
    tmp += second;
    return tmp;
  }

  /// Add two profile histograms
  inline Profile2D operator + (const Profile2D& first, const Profile2D& second) {
    return add(first,second);
  }

  /// Subtract two profile histograms
  inline Profile2D subtract(const Profile2D& first, const Profile2D& second) {
    Profile2D tmp = first;
    if (first.path() != second.path()) tmp.setPath("");
    tmp -= second;
    return tmp;
  }

  /// Subtract two profile histograms
  inline Profile2D operator - (const Profile2D& first, const Profile2D& second) {
    return subtract(first,second);
  }

  /// Divide two profile histograms
  Scatter3D divide(const Profile2D& numer, const Profile2D& denom);

  /// Divide two profile histograms
  inline Scatter3D operator / (const Profile2D& numer, const Profile2D& denom) {
    return divide(numer, denom);
  }

  /// @}

}

#endif
