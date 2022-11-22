// -*- C++ -*-
//
// This file is part of YODA -- Yet more Objects for Data Analysis
// Copyright (C) 2008-2021 The YODA collaboration (see AUTHORS for details)
//
#ifndef YODA_Histo1D_h
#define YODA_Histo1D_h

#include "YODA/AnalysisObject.h"
#include "YODA/Binned.h"
#include "YODA/Fillable.h"
#include "YODA/HistoBin1D.h"
#include "YODA/Dbn1D.h"
#include "YODA/Scatter2D.h"
#include "YODA/Axis1D.h"
#include "YODA/Exceptions.h"
#include <vector>
#include <string>
#include <map>

namespace YODA {


  /// Convenience typedef
  typedef Axis1D<HistoBin1D, Dbn1D> Histo1DAxis;

  /// A one-dimensional histogram
  class Histo1D : public AnalysisObject, public Binned, public Fillable {
  public:

    /// Convenience typedefs
    typedef Histo1DAxis Axis;
    typedef Axis::Bins Bins;
    typedef HistoBin1D Bin;

    typedef double FillType;
    typedef FillType BinType;
    typedef std::shared_ptr<Histo1D> Ptr;

    /// @name Constructors
    /// @{

    /// Default constructor
    Histo1D(const std::string& path="", const std::string& title="")
      : AnalysisObject("Histo1D", path, title),
        _axis()
    { }


    /// Constructor giving range and number of bins.
    Histo1D(size_t nbins, double lower, double upper,
            const std::string& path="", const std::string& title="")
      : AnalysisObject("Histo1D", path, title),
        _axis(nbins, lower, upper)
    { }


    /// @brief Constructor giving explicit bin edges.
    ///
    /// For n bins, binedges.size() == n+1, the last one being the upper bound
    /// of the last bin
    Histo1D(const std::vector<double>& binedges,
            const std::string& path="", const std::string& title="")
      : AnalysisObject("Histo1D", path, title),
        _axis(binedges)
    { }


    /// Constructor accepting an explicit collection of bins.
    Histo1D(const std::vector<Bin>& bins,
            const std::string& path="", const std::string& title="")
            : AnalysisObject("Histo1D", path, title),
            _axis(bins)
    { }


    /// Copy constructor with optional new path
    /// @todo Also allow title setting from the constructor?
    Histo1D(const Histo1D& h, const std::string& path="");


    /// Constructor from a Scatter2D's binning, with optional new path
    /// @todo Also allow title setting from the constructor?
    Histo1D(const Scatter2D& s, const std::string& path="");


    /// Constructor from a Profile1D's binning, with optional new path
    /// @todo Also allow title setting from the constructor?
    Histo1D(const Profile1D& p, const std::string& path="");


    /// @brief State-setting constructor
    ///
    /// Intended principally for internal persistency use.
    Histo1D(const std::vector<HistoBin1D>& bins,
            const Dbn1D& dbn_tot, const Dbn1D& dbn_uflow, const Dbn1D& dbn_oflow,
            const std::string& path="", const std::string& title="")
      : AnalysisObject("Histo1D", path, title),
        _axis(bins, dbn_tot, dbn_uflow, dbn_oflow)
    { }


    /// Assignment operator
    Histo1D& operator = (const Histo1D& h1) {
      AnalysisObject::operator = (h1); //< AO treatment of paths etc.
      _axis = h1._axis;
      return *this;
    }

    /// Make a copy on the stack
    Histo1D clone() const {
      return Histo1D(*this);
    }

    /// Make a copy on the heap, via 'new'
    Histo1D* newclone() const {
      return new Histo1D(*this);
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

    /// @brief Reset the histogram.
    ///
    /// Keep the binning but set all bin contents and related quantities to zero
    virtual void reset() {
      _axis.reset();
    }

    /// Fill histo by value and weight, optionally as a fractional fill
    virtual void fill(double x, double weight=1.0, double fraction=1.0);

    /// Fill histo bin i with the given weight, optionally as a fractional fill
    virtual void fillBin(size_t i, double weight=1.0, double fraction=1.0);


    /// Rescale as if all fill weights had been different by factor @a scalefactor.
    void scaleW(double scalefactor) {
      setAnnotation("ScaledBy", annotation<double>("ScaledBy", 1.0) * scalefactor);
      _axis.scaleW(scalefactor);
    }


    /// Normalize the (visible) histo area to the @a normto value.
    ///
    /// If @a includeoverflows is true, the original normalisation is computed with
    /// the overflow bins included, so that the resulting visible normalisation can
    /// be less than @a normto. This is probably what you want.
    void normalize(double normto=1.0, bool includeoverflows=true) {
      const double oldintegral = integral(includeoverflows);
      if (oldintegral == 0) throw WeightError("Attempted to normalize a histogram with null area");
      /// @todo Check that this is the desired behaviour
      scaleW(normto / oldintegral);
    }


    /// Merge together the bin range with indices from @a from to @a to, inclusive
    void mergeBins(size_t from, size_t to) {
      _axis.mergeBins(from, to);
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


    /// @name Bin accessors
    /// @{

    /// Number of bins (not counting under/overflow)
    size_t numBins() const { return bins().size(); }

    /// Number of bins on the x (only) axis (not counting under/overflow)
    size_t numBinsX() const { return numBins(); }

    /// Low edge of this histo's axis
    double xMin() const { return _axis.xMin(); }

    /// High edge of this histo's axis
    double xMax() const { return _axis.xMax(); }

    /// check if binning is the same as different Histo1D
    bool sameBinning(const Histo1D& h1) {
      return _axis == h1._axis;
    }

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
    std::vector<YODA::HistoBin1D>& bins() { return _axis.bins(); }
    /// Access the bin vector (const version)
    const std::vector<YODA::HistoBin1D>& bins() const { return _axis.bins(); }


    /// Access a bin by index (non-const version)
    HistoBin1D& bin(size_t index) { return _axis.bins()[index]; }
    /// Access a bin by index (const version)
    const HistoBin1D& bin(size_t index) const { return _axis.bins()[index]; }


    /// Access a bin index by coordinate
    /// @todo Convert to ssize_t?
    int binIndexAt(double x) {
      return _axis.binIndexAt(x);
    }

    /// Access a bin by coordinate (const version)
    const HistoBin1D& binAt(double x) const { return _axis.binAt(x); }


    /// Access summary distribution, including gaps and overflows (non-const version)
    Dbn1D& totalDbn() { return _axis.totalDbn(); }
    /// Access summary distribution, including gaps and overflows (const version)
    const Dbn1D& totalDbn() const { return _axis.totalDbn(); }
    /// Set summary distribution, mainly for persistency: CAREFUL!
    void setTotalDbn(const Dbn1D& dbn) { _axis.setTotalDbn(dbn); }


    /// Access underflow (non-const version)
    Dbn1D& underflow() { return _axis.underflow(); }
    /// Access underflow (const version)
    const Dbn1D& underflow() const { return _axis.underflow(); }
    /// Set underflow distribution, mainly for persistency: CAREFUL!
    void setUnderflow(const Dbn1D& dbn) { _axis.setUnderflow(dbn); }


    /// Access overflow (non-const version)
    Dbn1D& overflow() { return _axis.overflow(); }
    /// Access overflow (const version)
    const Dbn1D& overflow() const { return _axis.overflow(); }
    /// Set overflow distribution, mainly for persistency: CAREFUL!
    void setOverflow(const Dbn1D& dbn) { _axis.setOverflow(dbn); }

    /// @}


    /// @name Bin adding and removing
    /// @{

    /// Add a new bin specifying its lower and upper bound
    void addBin(double from, double to) { _axis.addBin(from, to); }

    /// Add new bins by specifying a vector of edges
    void addBins(std::vector<double> edges) { _axis.addBins(edges); }

    // /// Add new bins specifying a beginning and end of each of them
    // void addBins(std::vector<std::pair<double,double> > edges) {
    //   _axis.addBins(edges);
    // }

    /// Add a new bin, perhaps already populated: CAREFUL!
    void addBin(const HistoBin1D& b) { _axis.addBin(b); }

    /// @brief Bins addition operator
    ///
    /// Add multiple bins without resetting
    void addBins(const Bins& bins) {
      _axis.addBins(bins);
    }

    /// Remove a bin
    void rmBin(size_t index) { _axis.eraseBin(index); }

    /// @}


    /// @name Whole histo data
    /// @{

    /// Get the total area (sumW) of the histogram
    double integral(bool includeoverflows=true) const { return sumW(includeoverflows); }

    /// @brief Get the integrated area of the histogram between bins @a binindex1 and @a binindex2.
    ///
    /// @note The area of bin @a binindex2 _is_ included in the returned
    /// value. To include the underflow and overflow areas, you should add them
    /// explicitly with the underflow() and overflow() methods.
    ///
    /// @todo Allow int bin index args for type compatibility with binIndexAt()?
    double integralRange(size_t binindex1, size_t binindex2) const {
      assert(binindex2 >= binindex1);
      if (binindex1 >= numBins()) throw RangeError("binindex1 is out of range");
      if (binindex2 >= numBins()) throw RangeError("binindex2 is out of range");
      double rtn = 0;
      for (size_t i = binindex1; i <= binindex2; ++i) {
        rtn += bin(i).sumW();
      }
      return rtn;
    }

    /// @brief Get the integrated area of the histogram up to bin @a binindex.
    ///
    /// @note The area of bin @a binindex _is_ included in the returned
    /// value. To not include the underflow, set includeunderflow=false.
    ///
    /// @todo Allow int bin index args for type compatibility with binIndexAt()?
    double integralTo(size_t binindex, bool includeunderflow=true) const {
      double rtn = includeunderflow ? underflow().sumW() : 0;
      rtn += integralRange(0, binindex);
      return rtn;
    }

    /// Get the number of fills
    double numEntries(bool includeoverflows=true) const;

    /// Get the effective number of fills
    double effNumEntries(bool includeoverflows=true) const;

    /// Get sum of weights in histo
    double sumW(bool includeoverflows=true) const;

    /// Get sum of squared weights in histo
    double sumW2(bool includeoverflows=true) const;

    /// Get the mean in x
    double xMean(bool includeoverflows=true) const;

    /// Get the variance in x
    double xVariance(bool includeoverflows=true) const;

    /// Get the standard deviation in x
    double xStdDev(bool includeoverflows=true) const {
      if (includeoverflows) return _axis.totalDbn().xStdDev();
      return std::sqrt(xVariance(includeoverflows));
    }

    /// Get the standard error in x
    double xStdErr(bool includeoverflows=true) const;

    /// Get the RMS in x
    double xRMS(bool includeoverflows=true) const;

    /// @}


    /// @name Adding and subtracting histograms
    /// @{

    /// @brief Add another histogram to this one
    ///
    /// @note Adding histograms will unset any ScaledBy attribute from prevous calls to scaleW or normalize.
    Histo1D& operator += (const Histo1D& toAdd) {
      if (hasAnnotation("ScaledBy")) rmAnnotation("ScaledBy");
      _axis += toAdd._axis;
      return *this;

      // if (!hasAnnotation("ScaledBy") && !toAdd.hasAnnotation("ScaledBy")) {
      // _axis += toAdd._axis;
      // } else {
      //   // Undo scaling of both histograms
      //   double scaledBy = annotation<double>("ScaledBy", 1.0);
      //   _axis.scaleW(1.0/scaledBy);

      //   double toAddScaledBy = toAdd.annotation<double>("ScaledBy", 1.0);
      //   Axis1D<HistoBin1D, Dbn1D> toAddAxis = toAdd._axis;
      //   toAddAxis.scaleW(1.0/toAddScaledBy);

      //   _axis += toAddAxis;

      //   // Re-apply combined scaling
      //   double newScaledBy = scaledBy*toAddScaledBy/(scaledBy+toAddScaledBy);
      //   _axis.scaleW(newScaledBy);
      //   setAnnotation("ScaledBy", newScaledBy);
      // }
      /// @todo What about if one histo sets ScaledBy, and the other doesn't?!? Aaaargh
      // return *this;
    }

    /// @brief Subtract another histogram from this one
    ///
    /// @note Subtracting histograms will unset any ScaledBy attribute from prevous calls to scaleW or normalize.
    Histo1D& operator -= (const Histo1D& toSubtract) {
      if (hasAnnotation("ScaledBy")) rmAnnotation("ScaledBy");
      _axis -= toSubtract._axis;
      return *this;
    }

    /// @}


  protected:

    /// Access a bin by coordinate (non-const version)
    HistoBin1D& _binAt(double x) { return _axis.binAt(x); }


  private:

    /// @name Bin data
    /// @{

    /// Definition of bin edges and contents
    Axis1D<HistoBin1D, Dbn1D> _axis;

    /// @}

  };


  /// Convenience typedef
  typedef Histo1D H1D;


  /// @name Combining histos: global operators
  /// @{

  /// Add two histograms
  inline Histo1D add(const Histo1D& first, const Histo1D& second) {
    Histo1D tmp = first;
    if (first.path() != second.path()) tmp.setPath("");
    tmp += second;
    return tmp;
  }


  /// Add two histograms
  inline Histo1D operator + (const Histo1D& first, const Histo1D& second) {
    return add(first, second);
  }


  /// Subtract two histograms
  inline Histo1D subtract(const Histo1D& first, const Histo1D& second) {
    Histo1D tmp = first;
    if (first.path() != second.path()) tmp.setPath("");
    tmp -= second;
    return tmp;
  }


  /// Subtract two histograms
  inline Histo1D operator - (const Histo1D& first, const Histo1D& second) {
    return subtract(first, second);
  }


  /// @todo Add multiply(H1, H1) -> Scatter2D?


  /// Divide two histograms, with an uncorrelated error treatment
  ///
  /// @todo Wouldn't it be nice to be able to supply a correlation matrix or function as optional arg?
  ///
  /// @note The two histos must have _exactly_ the same binning.
  Scatter2D divide(const Histo1D& numer, const Histo1D& denom);


  /// Divide two histograms, with an uncorrelated error treatment
  ///
  /// @note The two histos must have _exactly_ the same binning.
  inline Scatter2D operator / (const Histo1D& numer, const Histo1D& denom) {
    return divide(numer, denom);
  }


  /// Add histogram and scatter
  Scatter2D add(const Histo1D& histo, const Scatter2D& scatt);

  inline Scatter2D add(const Scatter2D& scatt, const Histo1D& histo) {
    return add(histo, scatt);
  }

  /// Subtract scatter from histogram
  Scatter2D subtract(const Histo1D& histo, const Scatter2D& scatt);

  /// Subtract histogram from scatter
  Scatter2D subtract(const Scatter2D& scatt, const Histo1D& histo);

  /// Multiply histogram with scatter
  Scatter2D multiply(const Histo1D& histo, const Scatter2D& scatt);

  /// Multiply scatter with histogram
  inline Scatter2D multiply(const Scatter2D& scatt, const Histo1D& histo) {
    return multiply(histo, scatt);
  }

  /// Divide histogram by scatter
  Scatter2D divide(const Histo1D& numer, const Scatter2D& denom);

  /// Divide scatter by histogram
  Scatter2D divide(const Scatter2D& numer, const Histo1D& denom);

  /// @todo Add functions/operators on pointers



  /// @brief Calculate a histogrammed efficiency ratio of two histograms
  ///
  /// @note The two histos must have _exactly_ the same binning.
  ///
  /// @note An efficiency is not the same thing as a standard division of two
  /// histograms: the errors are treated as correlated via binomial statistics.
  Scatter2D efficiency(const Histo1D& accepted, const Histo1D& total);


  /// @brief Calculate the asymmetry (a-b)/(a+b) of two histograms
  ///
  /// @note The two histos must have _exactly_ the same binning.
  inline Scatter2D asymm(const Histo1D& a, const Histo1D& b) {
    return (a-b) / (a+b);
  }


  /// @brief Convert a Histo1D to a Scatter2D representing the integral of the histogram
  ///
  /// @note The integral histo errors are calculated as sqrt(binvalue), as if they
  /// are uncorrelated. This is not in general true for integral histograms, so if you
  /// need accurate errors you should explicitly monitor bin-to-bin correlations.
  ///
  /// The includeunderflow param chooses whether the underflow bin is included
  /// in the integral numbers as an offset.
  ///
  /// @todo Rename/alias as mkIntegral
  Scatter2D toIntegralHisto(const Histo1D& h, bool includeunderflow=true);


  /// @brief Convert a Histo1D to a Scatter2D where each bin is a fraction of the total
  ///
  /// @note This sounds weird: let's explain a bit more! Sometimes we want to
  /// take a histo h, make an integral histogram H from it, and then divide H by
  /// the total integral of h, such that every bin in H represents the
  /// cumulative efficiency of that bin as a fraction of the total. I.e. an
  /// integral histo, scaled by 1/total_integral and with binomial errors.
  ///
  /// The includeunderflow param behaves as for toIntegral, and applies to both
  /// the initial integration and the integral used for the scaling. The
  /// includeoverflow param applies only to obtaining the scaling factor.
  ///
  /// @todo Rename/alias as mkIntegralEff
  Scatter2D toIntegralEfficiencyHisto(const Histo1D& h, bool includeunderflow=true, bool includeoverflow=true);

  /// @}


}

#endif
