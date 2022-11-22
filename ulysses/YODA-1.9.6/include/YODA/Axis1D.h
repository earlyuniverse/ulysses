#ifndef YODA_Axis1D_h
#define YODA_Axis1D_h

#include "YODA/AnalysisObject.h"
#include "YODA/Exceptions.h"
#include "YODA/Bin.h"
#include "YODA/Utils/MathUtils.h"
#include "YODA/Utils/BinSearcher.h"
#include <limits>
#include <string>

namespace YODA {


  /// @brief 1D bin container
  ///
  /// This class handles most of the low-level operations on an axis of bins
  /// arranged in a 1D line (including gaps).
  template <typename BIN1D, typename DBN>
  class Axis1D {
  public:

    /// Typedefs
    /// @{

    /// Bin type
    typedef BIN1D Bin;

    /// A vector containing 1D bins. Not used for searching.
    typedef typename std::vector<Bin> Bins;

    /// @}

    /// @name Constructors
    /// @{

    /// Empty constructor
    Axis1D()
      : _locked(false)
    { }


    /// Constructor accepting a list of bin edges
    Axis1D(const std::vector<double>& binedges)
      : _locked(false)
    {
      addBins(binedges);
    }


    /// @brief Constructor accepting a vector of bins.
    ///
    /// Note that not only dimensions of these bins will be set,
    /// all the contents of the bins will be copied across, including
    /// the statistics
    Axis1D(const std::vector<BIN1D>& bins)
      : _locked(false)
    {
      addBins(bins);
    }


    /// Constructor with the number of bins and the axis limits
    /// @todo Rewrite interface to use a pair for the low/high
    Axis1D(size_t nbins, double lower, double upper)
      : _locked(false)
    {
      addBins(linspace(nbins, lower, upper));
    }


    /// @brief Constructor accepting a vector of bins.
    ///
    /// Note that not only dimensions of these bins will be set,
    /// all the contents of the bins will be copied across, including
    /// the statistics
    Axis1D(const Bins& bins, const DBN& dbn_tot, const DBN& dbn_uflow, const DBN& dbn_oflow)
      : _dbn(dbn_tot), _underflow(dbn_uflow), _overflow(dbn_oflow), _locked(false)
    {
      addBins(bins);
    }

    /// @name Statistics accessor functions
    /// @{

    /// Get the number of bins on the axis
    size_t numBins() const {
      return bins().size();
    }

    /// Return a vector of bins (const)
    const Bins& bins() const {
      return _bins;
    }

    /// Return a vector of bins (non-const)
    Bins& bins() {
      return _bins;
    }

    /// Return the lowest-value bin edge on the axis
    double xMin() const {
      if (numBins() == 0) throw RangeError("This axis contains no bins and so has no defined range");
      return bins().front().xMin();
    }

    /// Return the highest-value bin edge on the axis
    double xMax() const {
      if (numBins() == 0) throw RangeError("This axis contains no bins and so has no defined range");
      return bins().back().xMax();
    }

    /// Return all the Nbin+1 bin edges on the axis
    ///
    /// @note This only returns the finite edges, i.e. -inf and +inf are removed
    /// @todo Make the +-inf stripping controllable by a default-valued bool arg
    std::vector<double> xEdges() const {
      std::vector<double> rtn(_binsearcher.edges().begin()+1, _binsearcher.edges().end()-1);
      return rtn;
    }

    /// Return all the Nbin bin widths on the axis
    ///
    /// @note This only returns the finite edges, i.e. -inf and +inf are removed
    /// @todo Make the +-inf stripping controllable by a default-valued bool arg
    std::vector<double> xWidths() const {
      // std::vector<double> rtn; rtn.reserve(numBins());
      // for (const Bin& b : bins()) rtn.push_back( b.xWidth() );
      std::vector<double> rtn = xEdges();
      for (size_t ib = 0; ib < rtn.size()-1; ++ib) rtn[ib] = rtn[ib+1] - rtn[ib];
      rtn.pop_back();
      return rtn;
    }

    /// Return a bin at a given index (non-const)
    BIN1D& bin(size_t index) {
      if (index >= numBins()) throw RangeError("YODA::Histo1D: index out of range!");
      return _bins[index];
    }

    /// Return a bin at a given index (const)
    const BIN1D& bin(size_t index) const {
      if (index >= numBins()) throw RangeError("YODA::Histo1D: index out of range!");
      return _bins[index];
    }

    /// Returns an index of a bin at a given coord, -1 if no bin matches
    ssize_t binIndexAt(double coord) const {
      // Yes, this is robust even with an empty axis: there's always at least one outflow
      return _indexes[_binsearcher.index(coord)];
    }

    /// Return a bin at a given coordinate (non-const)
    BIN1D& binAt(double x) {
      const ssize_t index = binIndexAt(x);
      if (index == -1) throw RangeError("There is no bin at the specified x");
      return bin(index);
    }

    /// Return a bin at a given coordinate (const)
    const BIN1D& binAt(double x) const {
      const ssize_t index = binIndexAt(x);
      if (index == -1) throw RangeError("There is no bin at the specified x");
      return bin(index);
    }


    /// Return the total distribution (const)
    const DBN& totalDbn() const {
      return _dbn;
    }
    /// Return the total distribution (non-const)
    DBN& totalDbn() {
      return _dbn;
    }
    /// Set the total distribution: CAREFUL!
    void setTotalDbn(const DBN& dbn) {
      _dbn = dbn;
    }

    /// Return underflow (const)
    const DBN& underflow() const {
      return _underflow;
    }
    /// Return underflow (non-const)
    DBN& underflow() {
      return _underflow;
    }
    /// Set the underflow distribution: CAREFUL!
    void setUnderflow(const DBN& dbn) {
      _underflow = dbn;
    }

    /// Return overflow (const)
    const DBN& overflow() const {
      return _overflow;
    }
    /// Return overflow (non-const)
    DBN& overflow() {
      return _overflow;
    }
    /// Set the overflow distribution: CAREFUL!
    void setOverflow(const DBN& dbn) {
      _overflow = dbn;
    }

    /// @}


    /// @name Modifiers and helpers
    /// @{

    /// Reset all the bin statistics on the axis
    void reset() {
      _dbn.reset();
      _underflow.reset();
      _overflow.reset();
      for (Bin& bin : _bins) bin.reset();
      _locked = false;
    }


    /// Set the axis lock state
    void _setLock(bool locked) {
      _locked = locked;
    }


    /// Get lists of Nbins+1 edges, and Nbins+2 indices (-1 for gaps) from the current binning
    /// @todo Push this into an abstract Binning object in v2?
    /// @todo Should be a const ref arg?
    std::pair< std::vector<double>, std::vector<long> > _mk_edges_indexes(Bins& bins) const {
      std::vector<double> edges; edges.reserve(bins.size()+1); // Nbins+1 edges
      std::vector<long> indexes; edges.reserve(bins.size()+2); // Nbins + 2*outflows

      // Sort the bins
      std::sort(bins.begin(), bins.end());

      // Keep a note of the last high edge
      double last_high = -std::numeric_limits<double>::infinity();

      // Check for overlaps
      for (size_t i = 0; i < bins.size(); ++i) {
        Bin& currentBin = bins[i];
        const double new_low  = currentBin.xMin();
        const double reldiff = (new_low - last_high) / currentBin.xWidth();
        if (reldiff < -1e-3) { //< @note If there is a "large" negative gap (i.e. overlap), throw an exception
          std::stringstream ss;
          ss << "Bin edges overlap: " << last_high << " -> " << new_low;
          throw RangeError(ss.str());
        } else if (reldiff > 1e-3) { //< @note If there is a "large" positive gap, create a bin gap
          indexes.push_back(-1); // First index will be for underflow
          edges.push_back(new_low); // Will include first edge
        }

        // Bins check that they are not zero or negative width. It's okay for
        // them to throw an exception here, as we haven't changed anything yet.
        indexes.push_back(i);
        edges.push_back(currentBin.xMax());

        last_high = currentBin.xMax();
      }
      indexes.push_back(-1); // Overflow

      return std::make_pair(edges, indexes);
    }


    /// Merge together the bin range with indices from @a from to @a to, inclusive.
    /// Merge a series of bins, between the bins identified by indices @a from and @a to
    void mergeBins(size_t from, size_t to) {
      // Correctness checking
      if (from >= numBins())
        throw RangeError("Initial merge index is out of range");
      if (to >= numBins())
        throw RangeError("Final merge index is out of range");
      if (from > to)
        throw RangeError("Final bin must be greater than or equal to initial bin");
      if (_gapInRange(from, to))
        throw RangeError("Bin ranges containing gaps cannot be merged");
      if (from == to)
        return; // nothing to be done

      Bin& b = bin(from);
      for (size_t i = from + 1; i <= to; ++i)
        b.merge(_bins[i]);
      eraseBins(from+1, to);
    }


    /// @brief Merge every group of @a n bins, from start to end inclusive
    ///
    /// If the number of bins is not a multiple of @a n, the last @a m < @a n
    /// bins on the RHS will also be merged, as the closest possible approach to
    /// factor @n rebinning everywhere.
    void rebinBy(unsigned int n, size_t begin=0, size_t end=UINT_MAX) {
      if (n < 1) throw UserError("Rebinning requested in groups of 0!");
      for (size_t m = begin; m < end; ++m) {
        if (m > numBins()) break;
        const size_t myend = (m+n-1 < numBins()) ? m+n-1 : numBins()-1;
        if (myend > m) {
          mergeBins(m, myend);
          end -= myend-m; //< reduce upper index by the number of removed bins
        }
      }
    }

    /// @brief Overloaded alias for rebinBy
    void rebin(unsigned int n, size_t begin=0, size_t end=UINT_MAX) { rebinBy(n, begin, end); }


    /// @brief Rebin to the given list of bin edges
    void rebinTo(const std::vector<double>& newedges) {
      if (newedges.size() < 2)
        throw UserError("Requested rebinning to an edge list which defines no bins");
      const Utils::BinSearcher newbs(newedges);
      const std::vector<double> eshared = newbs.shared_edges(_binsearcher);
      if (eshared.size() != newbs.size())
        throw BinningError("Requested rebinning to incompatible edges");
      // std::cout << "Before merging" << std::endl;
      // for (double x : _binsearcher.edges()) std::cout << x << std::endl;
      // If the new min finite edge isn't the same, merge it into the underflow
      // NB. Edge search match the *next* bin, so step back one unit... and note these are BinSearcher indices, i.e. i+1
      if (!fuzzyEquals(xMin(), newedges.front())) {
        const size_t kmatch = _binsearcher.index(newedges.front()) - 1;
        mergeBins(0, kmatch-1);
        _underflow += bin(0).dbn();
        eraseBin(0);
      }
      // std::cout << "Merged start bins" << std::endl;
      // for (double x : _binsearcher.edges()) std::cout << x << std::endl;
      // Now the same for the overflow
      if (!fuzzyEquals(xMax(), newedges.back())) {
        const size_t kmatch = _binsearcher.index(newedges.back()) - 1;
        // std::cout << newedges.back() << " -> " << kmatch << " .. " << _bins.size()-1 << " / " << numBins() << std::endl;
        mergeBins(kmatch, _bins.size()-1);
        _overflow += bin(_bins.size()-1).dbn();
        eraseBin(_bins.size()-1);
      }
      // std::cout << "Merged end bins" << std::endl;
      // for (double x : _binsearcher.edges()) std::cout << x << std::endl;
      // Now merge the in-range bins
      size_t jcurr = 0;
      for (size_t i = 1; i < newedges.size(); ++i) { //< we already know that i=0 matches (until we support merging into overflows)
        const size_t kmatch = _binsearcher.index(newedges.at(i)) - 1; //< Will match the *next* bin, so step back one unit... and note these are BinSearcher indices
        assert(kmatch >= jcurr+1);
        mergeBins(jcurr, kmatch-1);
        jcurr += 1; //< The next bin to be merged, in the new numbering
      }
      // std::cout << "After merging" << std::endl;
      // for (double x : _binsearcher.edges()) std::cout << x << std::endl;
    }

    /// @brief Overloaded alias for rebinTo
    void rebin(const std::vector<double>& newedges) { rebinTo(newedges); }


    /// Add a bin, passed explicitly
    void addBin(const Bin& b) {
      /// @todo Efficiency?
      Bins newBins(_bins);
      newBins.push_back(b);
      _updateAxis(newBins);
    }


    /// Add a bin, providing its low and high edge
    void addBin(double low, double high) {
      addBin(Bin(low, high));
    }


    /// Add a contiguous set of bins to an axis, via their list of edges
    void addBins(const std::vector<double>& binedges) {
      Bins newBins(_bins);
      if (binedges.size() == 0) return;

      double low = binedges.front();
      for (size_t i = 1; i < binedges.size(); ++i) {
        const double high = binedges[i];
        assert(high>low); // Make sure binedges are meaningful
        newBins.push_back(Bin(low, high));
        low = high;
      }

      _updateAxis(newBins);
    }


    /// Add a list of bins as pairs of lowEdge, highEdge
    void addBins(const std::vector<std::pair<double, double> >& binpairs) {
      // Make a copy of the current binning
      Bins newBins(_bins);

      // Iterate over given bins
      for (size_t i = 0; i < binpairs.size(); ++i) {
        std::pair<double, double> b = binpairs[i];
        newBins.push_back(Bin(b.first, b.second));
      }
      _updateAxis(newBins);
    }


    /// Add a list of Bin objects
    void addBins(const Bins& bins) {
      Bins newBins(_bins);
      for (const Bin& b : bins) newBins.push_back(b);
      _updateAxis(newBins);
    }


    /// Remove a bin
    void eraseBin(const size_t i) {
      // Might as well erase from the internal bins, as we can guarantee
      // consistency.
      if (i >= numBins())
        throw RangeError("Bin index is out of range");

      const bool wasLocked = _locked;
      _locked = false;
      _bins.erase(_bins.begin() + i);
      _updateAxis(_bins);
      _locked = wasLocked;
    }


    /// Remove a bin range
    void eraseBins(const size_t from, const size_t to) {
      if (from >= numBins())
        throw RangeError("Initial index out of range");
      if (to >= numBins())
        throw RangeError("Final index out of range");
      if (from > to)
        throw RangeError("Final index is less than initial index");

      const bool wasLocked = _locked;
      _locked = false;
      _bins.erase(_bins.begin() + from, _bins.begin() + to + 1);
      _updateAxis(_bins);
      _locked = wasLocked;
    }


    /// Scale the size of an axis by a factor
    // @todo What if somebody passes in negative scalefactor? (test idea)
    void scaleX(double scalefactor) {
      _dbn.scaleX(scalefactor);
      _underflow.scaleX(scalefactor);
      _overflow.scaleX(scalefactor);
      for (size_t i = 0; i < _bins.size(); ++i)
        _bins[i].scaleX(scalefactor);
      _updateAxis(_bins);
    }


    /// Scale the amount of fills by a factor
    void scaleW(double scalefactor) {
      _dbn.scaleW(scalefactor);
      _underflow.scaleW(scalefactor);
      _overflow.scaleW(scalefactor);
      for (size_t i = 0; i < _bins.size(); ++i) _bins[i].scaleW(scalefactor);
    }

    /// @}


    /// @name Comparisons to other Axis objects
    /// @{

    bool sameBinning(const Axis1D& other) const {
      if (numBins() != other.numBins()) return false;
      if (_indexes != other._indexes) return false;
      return _binsearcher.same_edges(other._binsearcher);
    }

    bool subsetBinning(const Axis1D& other) const {
      const int ndiff = numBins() - other.numBins();
      if (ndiff == 0) return sameBinning(other);
      /// @todo Do we require the finite axis begin/end to be the same?
      return !_binsearcher.shared_edges(other._binsearcher).empty();
    }

    /// @}


    /// @name Operators
    /// @{

    /// Check if two of the Axis have the same binning, within numeric tolerance
    bool operator == (const Axis1D& other) const {
      return sameBinning(other);
    }


    /// Check if the binning of two of the Axis is different
    bool operator != (const Axis1D& other) const {
      return ! operator == (other);
    }


    /// Add two axes together
    Axis1D<BIN1D,DBN>& operator += (const Axis1D<BIN1D,DBN>& toAdd) {
      if (*this != toAdd) throw LogicError("YODA::Histo1D: Cannot add axes with different binnings.");

      for (size_t i = 0; i < _bins.size(); ++i) {
        _bins[i] += toAdd.bins().at(i);
      }

      _dbn += toAdd._dbn;
      _underflow += toAdd._underflow;
      _overflow += toAdd._overflow;
      return *this;
    }


    /// Subtract two axes
    Axis1D<BIN1D,DBN>& operator -= (const Axis1D<BIN1D,DBN>& toSubtract) {
      if (*this != toSubtract) throw LogicError("YODA::Histo1D: Cannot add axes with different binnings.");

      for (size_t i = 0; i < _bins.size(); ++i) {
        _bins[i] -= toSubtract.bins().at(i);
      }

      _dbn -= toSubtract._dbn;
      _underflow -= toSubtract._underflow;
      _overflow -= toSubtract._overflow;
      return *this;
    }

    /// @}


  private:

    /// Sort the given bins vector, and regenerate the bin searcher
    //
    /// The bin searcher is purely for searching, and is generated from
    /// the bins list only.
    void _updateAxis(Bins& bins) {
      // Ensure that axis is not locked
      if (_locked) {
        throw LockError("Attempting to update a locked 1D axis");
      }

      // Get the new cuts and indexes (throws if overlaps), and set them on the searcher
      const std::pair< std::vector<double>, std::vector<long> > es_is = _mk_edges_indexes(bins);
      _binsearcher = Utils::BinSearcher(es_is.first);
      _indexes = es_is.second;
      _bins = bins;
    }


    /// Check if there are any gaps in the axis' binning between bin indices @a
    /// from and @a to, inclusive.
    bool _gapInRange(size_t ifrom, size_t ito) const {
      if (ifrom == ito) return false;
      assert(ito < numBins() && ifrom < ito);

      /// @todo Why do we need to re-find the bin indices?
      const size_t from_ix = _binsearcher.index(bin(ifrom).xMid());
      const size_t to_ix = _binsearcher.index(bin(ito).xMid());
      // std::cout << ifrom << " vs. " << from_ix << std::endl;
      // std::cout << ito << " vs. " << to_ix << std::endl;

      for (size_t i = from_ix; i <= to_ix; i++)
      // for (size_t i = ifrom; i <= ito; i++)
        if (_indexes[i] == -1) return true;
      return false;
    }


  private:

    /// @name Data structures
    /// @{

    /// Bins vector
    Bins _bins;

    /// Total distribution
    DBN _dbn;

    /// Under- and overflows
    DBN _underflow, _overflow;

    // Binsearcher, for searching bins
    Utils::BinSearcher _binsearcher;

    // Mapping from binsearcher indices to bin indices (allowing gaps)
    std::vector<long> _indexes;

    /// Whether modifying bin edges is permitted
    bool _locked;

    /// @}

  };


  /// Add the statistics on two axes.
  template <typename BIN1D, typename DBN>
  inline Axis1D<BIN1D,DBN> operator + (const Axis1D<BIN1D,DBN>& first, const Axis1D<BIN1D,DBN>& second) {
    Axis1D<BIN1D,DBN> tmp = first;
    tmp += second;
    return tmp;
  }

  /// Subtract the statistics on two axis.
  template <typename BIN1D, typename DBN>
  inline Axis1D<BIN1D,DBN> operator - (const Axis1D<BIN1D,DBN>& first, const Axis1D<BIN1D,DBN>& second) {
    Axis1D<BIN1D,DBN> tmp = first;
    tmp -= second;
    return tmp;
  }


}

#endif
