// -*- C++ -*-
//
// This file is part of YODA -- Yet more Objects for Data Analysis
// Copyright (C) 2008-2021 The YODA collaboration (see AUTHORS for details)
//
#ifndef YODA_Histo2D_h
#define YODA_Histo2D_h

#include "YODA/AnalysisObject.h"
#include "YODA/Fillable.h"
#include "YODA/Binned.h"
#include "YODA/HistoBin2D.h"
#include "YODA/Dbn2D.h"
#include "YODA/Axis2D.h"
#include "YODA/Scatter3D.h"
#include "YODA/Exceptions.h"
#include <vector>
#include <tuple>

namespace YODA {


  // Forward declaration
  class Profile2D;
  class Scatter3D;

  /// Convenience typedef
  typedef Axis2D<HistoBin2D, Dbn2D> Histo2DAxis;


  /// A two-dimensional histogram
  class Histo2D : public AnalysisObject, public Fillable, public Binned {
  public:

    /// Convenience typedefs
    typedef Histo2DAxis Axis;
    typedef Axis::Bins Bins;
    typedef HistoBin2D Bin;
    typedef Axis::Outflows Outflows;

    typedef std::tuple<double, double> FillType;
    typedef FillType BinType;
    typedef std::shared_ptr<Histo2D> Ptr;


    /// @name Constructors
    /// @{

    /// Default constructor
    Histo2D(const std::string& path="", const std::string& title="")
      : AnalysisObject("Histo2D", path, title),
        _axis()
    { }


    /// Constructor giving range and number of bins.
    Histo2D(size_t nbinsX, double lowerX, double upperX,
            size_t nbinsY, double lowerY, double upperY,
            const std::string& path="", const std::string& title="")
      : AnalysisObject("Histo2D", path, title),
        _axis(nbinsX, std::make_pair(lowerX, upperX),  nbinsY, std::make_pair(lowerY, upperY))
    { }


    /// Constructor accepting the bin edges on X and Y axis.
    Histo2D(const std::vector<double>& xedges, const std::vector<double>& yedges,
            const std::string& path="", const std::string& title="")
            : AnalysisObject("Histo2D", path, title),
            _axis(xedges, yedges)
    { }


    /// Constructor accepting an explicit collection of bins.
    Histo2D(const std::vector<Bin>& bins,
            const std::string& path="", const std::string& title="")
            : AnalysisObject("Histo2D", path, title),
            _axis(bins)
    { }


    /// Copy constructor with optional new path
    /// @todo Also allow title setting from the constructor?
    Histo2D(const Histo2D& h, const std::string& path="");

    /// A constructor from a Scatter3D's binning, with optional new path
    /// @todo Also allow title setting from the constructor?
    Histo2D(const Scatter3D& s, const std::string& path="");

    /// Constructor from a Profile2D's binning, with optional new path
    /// @todo Also allow title setting from the constructor?
    Histo2D(const Profile2D& h, const std::string& path="");

    /// @brief State-setting constructor
    ///
    /// Mainly intended for internal persistency use.
    Histo2D(const std::vector<HistoBin2D>& bins,
            const Dbn2D& totalDbn,
            const Outflows& outflows,
            const std::string& path="", const std::string& title="")
      : AnalysisObject("Histo2D", path, title),
        _axis(bins, totalDbn, outflows)
    { }


    /// Assignment operator
    Histo2D& operator = (const Histo2D& h2) {
      AnalysisObject::operator = (h2); //< AO treatment of paths etc.
      _axis = h2._axis;
      return *this;
    }


    /// Make a copy on the stack
    Histo2D clone() const {
      return Histo2D(*this);
    }

    /// Make a copy on the heap, via 'new'
    Histo2D* newclone() const {
      return new Histo2D(*this);
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

    /// Fill histo with weight at (x,y)
    virtual void fill(double x, double y, double weight=1.0, double fraction=1.0);

    ///
    virtual void fill(const FillType & xs, double weight=1.0, double fraction=1.0) {
        fill(std::get<0>(xs), std::get<1>(xs), weight, fraction);
    }


    /// Fill histo x-y bin i with the given weight
    virtual void fillBin(size_t i, double weight=1.0, double fraction=1.0);


    /// @brief Reset the histogram.
    ///
    /// Keep the binning but set all bin contents and related quantities to zero
    void reset() {
      _axis.reset();
    }

    /// Rescale as if all fill weights had been different by factor @a scalefactor.
    void scaleW(double scalefactor) {
      setAnnotation("ScaledBy", annotation<double>("ScaledBy", 1.0) * scalefactor);
      _axis.scaleW(scalefactor);
    }


    /// Normalize the (visible) histo "volume" to the @a normto value.
    ///
    /// If @a includeoverflows is true, the original normalisation is computed with
    /// the overflow bins included, so that the resulting visible normalisation can
    /// be less than @a normto. This is probably what you want.
    void normalize(double normto=1.0, bool includeoverflows=true) {
      const double oldintegral = integral(includeoverflows);
      if (oldintegral == 0) throw WeightError("Attempted to normalize a histogram with null area");
      scaleW(normto / oldintegral);
    }


    /// Scale the dimensions
    void scaleXY(double scaleX = 1.0, double scaleY = 1.0) {
      _axis.scaleXY(scaleX, scaleY);
    }

    /// @}


    /// @name Bin adding and removing
    /// @{

    /// @brief Bin addition operator
    ///
    /// Add a bin to an axis described by its x and y ranges.
    void addBin(Axis::EdgePair1D xrange, Axis::EdgePair1D yrange) {
       _axis.addBin(xrange, yrange);
    }

    /// @brief Bin addition operator
    ///
    /// Add a bin, possibly already populated
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


    // /// Adding bins
    /// @todo TODO
    // void addBin(const std::vector<std::pair<std::pair<double,double>, std::pair<double,double> > > coords) {
    //     _axis.addBin(coords);
    // }

    // /// Adding bins which is not so eloquent
    /// @todo TODO
    // void addBin(double lowX, double lowY, double highX, double highY)   {
    //     _axis.addBin(lowX, lowY, highX, highY);
    // }

    // /// Merge the bins
    /// @todo TODO
    // void mergeBins(size_t from, size_t to) {
    //   _axis.mergeBins(from, to);
    // }

    /// Rebin the whole histo by a @a factorX in the X direction and
    /// @a factorY in the Y direction
    /// @todo TODO
    // void rebin(size_t factorX, size_t factorY){
    //   _axis.rebin(factorX, factorY);
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


    /// check if binning is the same as different Histo2D
    bool sameBinning(const Histo2D& h2) {
      return _axis == h2._axis;
    }

    /// Access the bin vector (non-const version)
    std::vector<YODA::HistoBin2D>& bins() { return _axis.bins(); }
    /// Access the bin vector (const version)
    const std::vector<YODA::HistoBin2D>& bins() const { return _axis.bins(); }


    /// Access a bin by index (non-const version)
    HistoBin2D& bin(size_t index) { return _axis.bin(index); }
    /// Access a bin by index (const version)
    const HistoBin2D& bin(size_t index) const { return _axis.bin(index); }


    /// Access a bin index by coordinate
    int binIndexAt(double x, double y) { return _axis.binIndexAt(x, y); }

    int binIndexAt(const BinType& t) { return _axis.binIndexAt(std::get<0>(t), std::get<1>(t)); }

    /// Access a bin by coordinate (const version)
    const HistoBin2D& binAt(double x, double y) const { return _axis.binAt(x, y); }

    const HistoBin2D& binAt(const BinType& t) { return _axis.binAt(std::get<0>(t), std::get<1>(t)); }


    /// Number of bins
    size_t numBins() const { return _axis.numBins(); }

    /// Number of bins along the x axis
    size_t numBinsX() const { return _axis.numBinsX(); }

    /// Number of bins along the y axis
    size_t numBinsY() const { return _axis.numBinsY(); }


    /// Access summary distribution, including gaps and overflows (non-const version)
    Dbn2D& totalDbn() { return _axis.totalDbn(); }
    /// Access summary distribution, including gaps and overflows (const version)
    const Dbn2D& totalDbn() const { return _axis.totalDbn(); }
    /// Set summary distribution, including gaps and overflows
    void setTotalDbn(const Dbn2D& dbn) { _axis.setTotalDbn(dbn); }


    // /// @brief Access an outflow (non-const)
    // ///
    // /// Two indices are used, for x and y: -1 = underflow, 0 = in-range, and +1 = overflow.
    // /// (0,0) is not a valid overflow index pair, since it is in range for both x and y.
    // Dbn2D& outflow(int ix, int iy) {
    //   std::cout << "Histo2D::outflow\n";
    //   return _axis.outflow(ix, iy);
    // }

    // /// @brief Access an outflow (const)
    // ///
    // /// Two indices are used, for x and y: -1 = underflow, 0 = in-range, and +1 = overflow.
    // /// (0,0) is not a valid overflow index pair, since it is in range for both x and y.
    // const Dbn2D& outflow(int ix, int iy) const {
    //   return _axis.outflow(ix, iy);
    // }

    /// @}


    /// @name Whole histo data
    /// @{

    /// Get the total volume of the histogram
    double integral(bool includeoverflows=true) const { return sumW(includeoverflows); }

    /// Get the number of fills (fractional fills are possible)
    double numEntries(bool includeoverflows=true) const;

    /// Get the effective number of fills
    double effNumEntries(bool includeoverflows=true) const;

    /// Get the sum of weights in histo
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

    /// Get the standard error in x
    double xStdErr(bool includeoverflows=true) const;

    /// Get the standard error in y
    double yStdErr(bool includeoverflows=true) const;

    /// Get the RMS in x
    double xRMS(bool includeoverflows=true) const;

    /// Get the RMS in y
    double yRMS(bool includeoverflows=true) const;

    /// @}


    /// @name Adding and subtracting histograms
    /// @{

    /// @brief Add another histogram to this one
    ///
    /// @note Adding histograms will unset any ScaledBy attribute from prevous calls to scaleW or normalize.
    Histo2D& operator += (const Histo2D& toAdd) {
      if (hasAnnotation("ScaledBy")) rmAnnotation("ScaledBy");
      _axis += toAdd._axis;
      return *this;
    }

    /// @brief Subtract another histogram from this one
    ///
    /// @note Subtracting histograms will unset any ScaledBy attribute from prevous calls to scaleW or normalize.
    Histo2D& operator -= (const Histo2D& toSubtract) {
      if (hasAnnotation("ScaledBy")) rmAnnotation("ScaledBy");
      _axis -= toSubtract._axis;
      return *this;
    }

    bool operator == (const Histo2D& other) const {
      return _axis == other._axis;
    }

    bool operator != (const Histo2D& other) const {
        return ! operator == (other);
    }

    /// @}


    // /// @name Slicing operators
    // /// @{

    // /// @brief Create a Histo2D for the bin slice parallel to the x axis at the specified y coordinate
    // ///
    // /// Note that the created histogram will not have correctly filled underflow and overflow bins.
    // /// @todo It's not really *at* the specified y coord: it's for the corresponding bin row.
    // /// @todo Change the name!
    // Histo2D cutterX(double atY, const std::string& path="", const std::string& title="");


    // /// @brief Create a Histo2D for the bin slice parallel to the y axis at the specified x coordinate
    // ///
    // /// Note that the created histogram will not have correctly filled underflow and overflow bins.
    // /// @todo It's not really *at* the specified x coord: it's for the corresponding bin row.
    // /// @todo Change the name!
    // Histo2D cutterY(double atX, const std::string& path="", const std::string& title="");


    // /// X-wise Profile1D creator from Histo2D
    // Profile1D mkProfileX();

    // /// Y-wise Profile1D creator from Histo2D
    // Profile1D mkProfileY();
    // /// @}


  protected:

    /// Access a bin by coordinate (non-const version)
    HistoBin2D& _binAt(double x, double y) { return _axis.binAt(x, y); }


  private:

    /// @name Bin data
    /// @{

    /// Definition of bin edges and contents
    Axis2D<HistoBin2D, Dbn2D> _axis;

    /// @}

  };


  /// Convenience typedef
  typedef Histo2D H2D;


  /// @name Combining histos: global operators
  /// @{

  /// Add two histograms
  inline Histo2D add(const Histo2D& first, const Histo2D& second) {
    Histo2D tmp = first;
    if (first.path() != second.path()) tmp.setPath("");
    tmp += second;
    return tmp;
  }


  /// Add two histograms
  inline Histo2D operator + (const Histo2D& first, const Histo2D& second) {
    return add(first, second);
  }


  /// Subtract two histograms
  inline Histo2D subtract(const Histo2D& first, const Histo2D& second) {
    Histo2D tmp = first;
    if (first.path() != second.path()) tmp.setPath("");
    tmp -= second;
    return tmp;
  }


  /// Subtract two histograms
  inline Histo2D operator - (const Histo2D& first, const Histo2D& second) {
    return subtract(first, second);
  }


  /// @todo Add multiply(H2, H2) -> Scatter3D?


  /// @brief Divide two histograms, with an uncorrelated error treatment
  ///
  /// @todo Wouldn't it be nice to be able to supply a correlation matrix or function as optional arg?
  ///
  /// @note The two histos must have _exactly_ the same binning.
  Scatter3D divide(const Histo2D& numer, const Histo2D& denom);


  /// Divide two histograms, with an uncorrelated error treatment
  ///
  /// @note The two histos must have _exactly_ the same binning.
  inline Scatter3D operator / (const Histo2D& numer, const Histo2D& denom) {
    return divide(numer, denom);
  }


  /// @brief Calculate a histogrammed efficiency ratio of two histograms
  ///
  /// @note The two histos must have _exactly_ the same binning.
  ///
  /// @note An efficiency is not the same thing as a standard division of two
  /// histograms: the errors are treated as correlated via binomial statistics.
  Scatter3D efficiency(const Histo2D& accepted, const Histo2D& total);


  /// @brief Calculate the asymmetry (a-b)/(a+b) of two histograms
  ///
  /// @note The two histos must have _exactly_ the same binning.
  inline Scatter3D asymm(const Histo2D& a, const Histo2D& b) {
    return (a-b) / (a+b);
  }

  /// @}


}

#endif
