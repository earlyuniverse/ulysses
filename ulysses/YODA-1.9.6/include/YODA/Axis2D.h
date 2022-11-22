#ifndef YODA_Axis2D_h
#define YODA_Axis2D_h

#include "YODA/AnalysisObject.h"
#include "YODA/Exceptions.h"
#include "YODA/Bin.h"
#include "YODA/Utils/MathUtils.h"
#include "YODA/Utils/Predicates.h"
#include "YODA/Utils/BinSearcher.h"
#include <limits>
#include <string>

namespace YODA {


  /// @brief 2D bin container
  ///
  /// This class handles most of the low-level operations on an axis of bins
  /// arranged in a 2D grid (including gaps).
  template <typename BIN2D, typename DBN>
  class Axis2D {
  public:

    /// Typedefs
    /// @{

    /// Bin type
    typedef BIN2D Bin;

    /// A vector containing 2D bins. Not used for searching.
    typedef typename std::vector<Bin> Bins;

    // Distinguishing between single edges and edge pairs (and pairs of pairs) is useful
    typedef std::vector<double> Edges;
    typedef std::pair<double, double> EdgePair1D;
    typedef std::pair<EdgePair1D, EdgePair1D> EdgePair2D;
    typedef std::vector<EdgePair2D> EdgePair2Ds;

    // Outflow distribution lists: see outflow(int, int)
    typedef std::vector<DBN> Outflow;
    typedef std::vector<Outflow> Outflows;

    /// @}

    /// @name Constructors
    /// @{

    // Empty constructor
    Axis2D()
      : _locked(false)
    {
      reset();
    }

    /// A constructor with specified x and y axis bin cuts.
    Axis2D(const Edges& xedges, const Edges& yedges)
      : _locked(false)
    {
      addBins(xedges, yedges);
      reset();
    }

    /// Constructor accepting X/Y ranges and number of bins
    /// on each of the axis. Both axes are divided linearly.
    Axis2D(size_t nbinsX, const std::pair<double,double>& rangeX,
           size_t nbinsY, const std::pair<double,double>& rangeY)
      : _locked(false)
    {
      addBins(linspace(nbinsX, rangeX.first, rangeX.second),
              linspace(nbinsY, rangeY.first, rangeY.second));
      reset();
    }

    /// Constructor accepting a list of bins
    Axis2D(const Bins& bins)
      : _locked(false)
    {
      addBins(bins);
      reset();
    }

    /// State-setting constructor for persistency
    Axis2D(const Bins& bins,
           const DBN& totalDbn,
           const Outflows& outflows)
      : _dbn(totalDbn), _outflows(outflows),
        _locked(false) // Does this make sense?
    {
      if (_outflows.size() != 8) {
        throw Exception("Axis2D outflow containers must have exactly 8 elements");
      }
      addBins(bins);
    }


    void reset() {
      _dbn.reset();
      _outflows.assign(8, Outflow());
      for (Bin& bin : _bins) bin.reset();
      _locked = false;
    }


    /// Get the number of bins.
    size_t numBins() const {
      return _bins.size();
    }

    /// Get the number of bins on the x-axis. This is only sensible for
    /// perfectly regular gridded bins. For irregular binnings, this is
    /// the number of cuts that were necessary to grid the data.
    size_t numBinsX() const {
      return _nx-1;
    }

    /// Get the number of bins on the y-axis. This is only sensible for
    /// perfectly regular gridded bins. For irregular binnings, this is
    /// the number of cuts that were necessary to grid the data.
    size_t numBinsY() const {
      return _ny-1;
    }

    /// @}
    //
    /// @name Statistics accessor functions
    /// @{

    /// @brief Get the outflow by x-index and y-index (non-const version)
    ///
    /// Indices are -1 = below range, 0 = in range, +1 = above range, e.g. (+1,
    /// -1) is in the "bottom right" position by being greater than the greatest
    /// x-edge and less than the lowest y-edge.
    ///
    Outflow& outflow(int ix, int iy) {
      return _outflows[_outflowIndex(ix, iy)];
    }

    /// @brief Get the outflow by x-index and y-index (const version)
    ///
    /// Indices are -1 = below range, 0 = in range, +1 = above range, e.g. (+1,
    /// -1) is in the "bottom right" position by being greater than the greatest
    /// x-edge and less than the lowest y-edge.
    ///
    const Outflow& outflow(int ix, int iy) const {
      return _outflows[_outflowIndex(ix, iy)];
    }

    /// Scale each bin as if the entire x-axis had been scaled by this factor.
    void scaleX(double xscale) {
      scaleXY(xscale, 1.0);
    }

    /// Scale each bin as if the entire y-axis had been scaled by this factor.
    void scaleY(double yscale) {
      scaleXY(1.0, yscale);
    }

    /// Scale each bin as if the entire x and y-axes had been scaled by
    /// their respective factors.
    void scaleXY(double sx, double sy) {
      _dbn.scaleXY(sx, sy);
      for (Outflow& outflow : _outflows)
        for (DBN& dbn : outflow)
          dbn.scaleXY(sx, sy);
      for (Bin& b : _bins)
        b.scaleXY(sx, sy);
      _updateAxis(_bins);
    }


    /// Rescale as if all fill weights had been different by factor @a
    /// scalefactor.
    void scaleW(double scalefactor) {
      _dbn.scaleW(scalefactor);
      for (Outflow& outflow : _outflows)
        for (DBN& dbn : outflow)
          dbn.scaleW(scalefactor);
      for (Bin& bin : _bins)
        bin.scaleW(scalefactor);
      _updateAxis(_bins);
    }


    /// Remove the bin at the given index. If many bins need to be
    /// removed, prefer eraseBins(vector[size_t] &) over many calls to this,
    /// as recreating the binhash is comparatively expensive.
    void eraseBin(size_t i) {
      if (i >= numBins())
        throw RangeError("Bin index is out of range");

      // Temporarily unlock the axis during the update
      _bins.erase(_bins.begin() + i);
      _updateAxis(_bins);
    }


    /// Erase a rectangle of bins.
    void eraseBins(const size_t from, const size_t to) {
      if (from >= numBins())
        throw RangeError("Initial bin index is out of range");
      if (from >= numBins())
        throw RangeError("Final bin index is out of range");

      Bin& fromBin = bin(from);
      Bin& toBin = bin(to);

      eraseBins(std::make_pair(fromBin.xMin(), toBin.xMax()),
                std::make_pair(fromBin.yMin(), toBin.yMax()));
    }

    /// Erase bins in an x- and y-range. Any bins which lie entirely within the
    /// region are deleted. If any part of the bin lies outside this
    /// range, the bin remains, so this has similar behaviour to select
    /// tools in vector graphics GUI packages.
    ///
    /// @todo How to test this?
    void eraseBins(const std::pair<double, double>& xrange,
                   const std::pair<double, double>& yrange)
    {
      size_t xiLow = _binSearcherX.index(xrange.first) - 1;
      size_t xiHigh = _binSearcherX.index(xrange.second) - 1;

      size_t yiLow = _binSearcherY.index(yrange.first) - 1;
      size_t yiHigh = _binSearcherY.index(yrange.second) - 1;

      /// @todo Beware the specialisation problems with vector<bool>...
      std::vector<bool> deleteMask(numBins(), false);

      for (size_t yi = yiLow; yi < yiHigh; yi++) {
        for (size_t xi = xiLow; xi < xiHigh; xi++) {
          ssize_t i = _indexes[_index(_nx, xi, yi)];
          if (i == -1 || deleteMask[i]) continue;
          if (bin(i).fitsInside(xrange, yrange)) deleteMask[i] = true;
        }
      }

      // Now we just update
      eraseBins(deleteMask);
    }


    /// Erase using a vector<bool>, where true represents that a bin
    /// will be deleted, and false means it will be kept.
    void eraseBins(const std::vector<bool>& deleteMask) {
      Bins newBins;
      for (size_t i = 0; i < numBins(); i++)
        if (!deleteMask[i]) newBins.push_back(bins(i));
      _update(newBins);
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


    /// Rebin with the same rebinning factor @a n in x and y
    void rebin(unsigned int n) {
      rebinXY(n, n);
    }

    /// Rebin with separate rebinning factors @a nx, @a ny in x and y
    void rebinXY(unsigned int nx, unsigned int ny) {
      rebinX(nx);
      rebinY(ny);
    }

    /// Rebin in x by factor @a nx
    void rebinX(unsigned int nx) {
      /// @todo WRITE THIS!
    }

    /// Rebin in y by factor @a ny
    void rebinY(unsigned int ny) {
      /// @todo WRITE THIS!
    }


    /// Set the axis lock state
    /// @todo Remove? Should not be public
    void _setLock(bool locked) { _locked = locked; }


    /// @todo Add xMins, xMaxs, xMids, xFoci, and y-versions


    /// Return the lowest-valued bin edge along the x-axis
    double xMin() const { return _xRange.first; }

    /// Return the highest-valued bin edge along the x-axis
    double xMax() const { return _xRange.second; }


    /// Return the lowest-valued bin edge along the y-axis
    double yMin() const { return _yRange.first; }

    /// Return the highest-valued bin edge along the y-axis
    double yMax() const { return _yRange.second; }


    /// Return all the NbinX+1 bin edges on the x-axis
    ///
    /// @note This only returns the finite edges, i.e. -inf and +inf are removed
    /// @todo Make the +-inf stripping controllable by a default-valued bool arg
    std::vector<double> xEdges() const {
      std::vector<double> rtn(_binSearcherX.edges().begin()+1, _binSearcherX.edges().end()-1);
      return rtn;
    }

    /// Return all the NbinY+1 bin edges on the y-axis
    ///
    /// @note This only returns the finite edges, i.e. -inf and +inf are removed
    /// @todo Make the +-inf stripping controllable by a default-valued bool arg
    std::vector<double> yEdges() const {
      std::vector<double> rtn(_binSearcherY.edges().begin()+1, _binSearcherY.edges().end()-1);
      return rtn;
    }

    /// Return all the NbinX bin widths on the x-axis
    ///
    /// @note This only returns the finite edges, i.e. -inf and +inf are removed
    /// @todo Make the +-inf stripping controllable by a default-valued bool arg
    std::vector<double> xWidths() const {
      std::vector<double> rtn = xEdges();
      for (size_t ib = 0; ib < rtn.size()-1; ++ib) rtn[ib] = rtn[ib+1] - rtn[ib];
      rtn.pop_back();
      return rtn;
    }

    /// Return all the NbinY bin widths on the y-axis
    ///
    /// @note This only returns the finite edges, i.e. -inf and +inf are removed
    /// @todo Make the +-inf stripping controllable by a default-valued bool arg
    std::vector<double> yWidths() const {
      std::vector<double> rtn = yEdges();
      for (size_t ib = 0; ib < rtn.size()-1; ++ib) rtn[ib] = rtn[ib+1] - rtn[ib];
      rtn.pop_back();
      return rtn;
    }


    /// Add a bin, providing its x- and y- edge ranges
    void addBin(EdgePair1D xrange, EdgePair1D yrange) {
      _checkUnlocked();
      Bins newBins = _bins;
      newBins.push_back(Bin(xrange, yrange));
      _updateAxis(newBins);
    }

    /// Add a pre-made bin
    void addBin(const Bin& bin) {
      _checkUnlocked();
      Bins newBins = _bins;
      newBins.push_back(bin);
      _updateAxis(newBins);
    }

    /// Add a vector of pre-made bins
    void addBins(const Bins& bins) {
      if (bins.size() == 0) return;
      _checkUnlocked();
      Bins newBins = _bins;
      for (const Bin& b : bins)
        newBins.push_back(b);
      _updateAxis(newBins);
    }

    /// Add a contiguous set of bins to an axis, via their list of edges
    void addBins(const std::vector<double>& xedges, const std::vector<double>& yedges) {
      if (xedges.size() == 0) return;
      if (yedges.size() == 0) return;
      _checkUnlocked();

      Bins newBins = _bins;
      for (size_t xi = 0; xi < xedges.size()-1; xi++) {
        for (size_t yi = 0; yi < yedges.size()-1; yi++) {
          const std::pair<double,double> xx = std::make_pair(xedges[xi], xedges[xi+1]);
          const std::pair<double,double> yy = std::make_pair(yedges[yi], yedges[yi+1]);
          // std::cout << "New bin with edges: [(" << xx.first << "," << xx.second << "), " << yy.first << "," << yy.second << ")]" << std::endl;
          newBins.push_back(Bin(xx, yy));
        }
      }

      _updateAxis(newBins);
    }


    /// Access bin by index
    Bin& bin(size_t i) {
      return _bins[i];
    }

    /// Access bin by index (const)
    const Bin& bin(size_t i) const {
      return _bins[i];
    }

    /// Get the bin index of the bin containing point (x, y).
    int binIndexAt(double x, double y) const {
      size_t xi = _binSearcherX.index(x) - 1;
      size_t yi = _binSearcherY.index(y) - 1;
      if (xi > _nx) return -1;
      if (yi > _ny) return -1;

      return _indexes[_index(_nx, xi, yi)];
    }

    /// Get the bin containing point (x, y).
    Bin& binAt(double x, double y) {
      const int ret = binIndexAt(x, y);
      if (ret == -1) throw RangeError("No bin found!!");
      return bin(ret);
    }

    /// Get the bin containing point (x, y) (const).
    const Bin& binAt(double x, double y) const {
      const int ret = binIndexAt(x, y);
      if (ret == -1) throw RangeError("No bin found!!");
      return bin(ret);
    }


    /// Return the total distribution (non-const)
    DBN& totalDbn() {
      return _dbn;
    }
    /// Return the total distribution (const)
    const DBN& totalDbn() const {
      return _dbn;
    }
    /// Set the total distribution: CAREFUL!
    void setTotalDbn(const DBN& dbn) {
      _dbn = dbn;
    }


    /// Return the bins vector (non-const)
    Bins& bins() {
      return _bins;
    }

    /// Return the bins vector (const)
    const Bins& bins() const {
      return _bins;
    }


    /// Equality operator (on binning only)
    /// @todo Change as discussed below if we expose the Axis classes for direct use
    // (DM: Doesn't this break the semantics of equality?  As it's used only
    // rarely, isn't there a real case for having a "binningsCompatible" or
    // similar method?)
    bool operator == (const Axis2D& other) const {
      if (numBins() != other.numBins()) return false;
      for (size_t i = 0; i < numBins(); i++)
        if (!(fuzzyEquals(bin(i).xMin(), other.bin(i).xMin()) &&
              fuzzyEquals(bin(i).xMax(), other.bin(i).xMax()) &&
              fuzzyEquals(bin(i).yMin(), other.bin(i).yMin()) &&
              fuzzyEquals(bin(i).yMax(), other.bin(i).yMax())))
          return false;
      return true;
    }

    /// Non-equality operator
    bool operator != (const Axis2D& other) const {
      return ! operator == (other);
    }


    /// Addition operator
    Axis2D<BIN2D, DBN>& operator += (const Axis2D<BIN2D, DBN>& toAdd) {
      if (*this != toAdd) {
        throw LogicError("YODA::Axis2D: Cannot add axes with different binnings.");
      }
      for (size_t i = 0; i < bins().size(); ++i) {
        bin(i) += toAdd.bin(i);
      }
      _dbn += toAdd._dbn;
      return *this;
    }

    /// Subtraction operator
    Axis2D<BIN2D, DBN>& operator -= (const Axis2D<BIN2D, DBN>& toSubtract) {
      if (*this != toSubtract) {
        throw LogicError("YODA::Axis2D: Cannot add axes with different binnings.");
      }
      for (size_t i = 0; i < bins().size(); ++i) {
        bin(i) -= toSubtract.bin(i);
      }
      _dbn -= toSubtract._dbn;
      return *this;
    }


  private:

    void _checkUnlocked(void) {
      // Ensure that axis is not locked
      if (_locked)
        throw LockError("Attempting to update a locked 2D axis");
    }


    /// Detect if there is a binning gap in the given bin index range
    /// @todo WRITE THIS!
    bool _gapInRange(size_t from, size_t to) {
      Bin& toBin = bin(to);
      Bin& fromBin = bin(from);
      return true;
    }


    void _updateAxis(Bins& bins) {
      // Deal with the case that there are no bins supplied (who called that?!)
      if (bins.size() == 0) {
        _binSearcherX = Utils::BinSearcher();
        _binSearcherY = Utils::BinSearcher();
        _nx = 0;
        _ny = 0;
        _xRange = std::make_pair(0, 0);
        _yRange = std::make_pair(0, 0);
      }

      // Sort the bins
      std::sort(bins.begin(), bins.end());

      // Create the edges
      std::vector<double> xedges, yedges, xwidths, ywidths;
      for (const Bin& bin : bins) {
        xedges.push_back(bin.xMin());
        xedges.push_back(bin.xMax());
        xwidths.push_back(bin.xWidth());
        yedges.push_back(bin.yMin());
        yedges.push_back(bin.yMax());
        ywidths.push_back(bin.yWidth());
      }

      // Sort the edges and widths
      std::sort(xedges.begin(), xedges.end());
      std::sort(yedges.begin(), yedges.end());
      std::sort(xwidths.begin(), xwidths.end());
      std::sort(ywidths.begin(), ywidths.end());

      // Obtain the median widths as a typical scale for uniqueness comparisons
      // const double medianxwidth = xwidths[ (xwidths.size()-1)/2 ];
      // const double medianywidth = ywidths[ (ywidths.size()-1)/2 ];
      const double minxwidth = xwidths[0];
      const double minywidth = ywidths[0];

      // Uniqueify the bin edges in the x- and y-cut vectors, with some numerical fuzziness
      xedges.resize(std::unique(xedges.begin(), xedges.end(), CmpFloats(1e-3, minxwidth)) - xedges.begin());
      yedges.resize(std::unique(yedges.begin(), yedges.end(), CmpFloats(1e-3, minywidth)) - yedges.begin());

      const size_t nx = xedges.size();
      const size_t ny = yedges.size();
      const size_t N = nx * ny;
      //std::cout << "Unique Axis2D edge list sizes: nx = " << nx << ", ny = " << ny << std::endl;
      assert(bins.size() <= (nx-1)*(ny-1) && "Input bins vector size must agree with computed number of unique bins");

      // Create a sea of indices, starting with an all-gaps configuration
      std::vector<ssize_t> indexes(N, -1);

      // Iterate through bins and find out which
      Utils::BinSearcher xSearcher(xedges);
      Utils::BinSearcher ySearcher(yedges);
      for (size_t i = 0; i < bins.size(); ++i) {
        Bin& bin = bins[i];

        // std::cout << "Bin #" << i << " edges: "
        //           << "[(" << bin.xMin() << "," << bin.xMax() << "), "
        //           << "(" << bin.yMin() << "," << bin.yMax() << ")] " << std::endl;

        const size_t xiMin= xSearcher.index(bin.xMin()) - 1;
        const size_t xiMax= xSearcher.index(bin.xMax()) - 1;
        const size_t yiMin = ySearcher.index(bin.yMin()) - 1;
        const size_t yiMax = ySearcher.index(bin.yMax()) - 1;

        // std::cout << "Sub-bin range indices: x = " << xiMin << ".." << xiMax << ", y = " << yiMin << ".." << yiMax << std::endl;

        // Loop over sub-bins in the edge list and assign indices / detect overlaps
        for (size_t xi = xiMin; xi < xiMax; xi++) {
          for (size_t yi = yiMin; yi < yiMax; yi++) {
            const size_t ii = _index(nx, xi, yi);
            if (indexes[ii] != -1) {
              std::stringstream ss;
              ss << "Bin edges overlap! Bin #" << i << " with edges "
                 << "[(" << bin.xMin() << "," << bin.xMax() << "), "
                 << "(" << bin.yMin() << "," << bin.yMax() << ")] "
                 << "overlaps bin #" << indexes[ii] << " in sub-bin #" << ii;
              throw RangeError(ss.str());
            }
            indexes[ii] = i;
          }
        }
      }

      // Job's a good'n - let's change our class.
      _nx = nx;
      _ny = ny;

      _xRange = std::make_pair(xedges.front(), xedges.back());
      _yRange = std::make_pair(yedges.front(), yedges.back());

      _indexes = indexes;
      _bins = bins;

      _binSearcherX = xSearcher;
      _binSearcherY = ySearcher;
    }


    /// Definition of global bin ID in terms of x and y bin IDs
    static size_t _index(size_t nx, size_t x, size_t y) {
      return y * nx + x;
    }

    /// @brief Get the outflow array index by x-index and y-index
    ///
    /// Indices are -1 = below range, 0 = in range, +1 = above range, e.g. (+1,
    /// -1) is in the "bottom right" position by being greater than the greatest
    /// x-edge and less than the lowest y-edge.
    static size_t _outflowIndex(int ix, int iy) {
      if (ix == 0 || iy == 0)
        throw UserError("The in-range (0,0) index pair is not a valid outflow specifier");
      ix += 1;
      iy += 1;
      if (ix > 2 || iy > 2)
        throw UserError("Outflow index out of range: valid indices are -1, 0, 1");
      size_t rtn = 3*ix + iy; // uncorrected for (0,0) index offset
      if (rtn > 4) rtn -= 1; // offset correction (note that raw rtn == 4 is not possible)
      return rtn;
    }


    /// @name Data structures
    /// @{

    /// Bins vector
    Bins _bins;

    /// Total distribution
    DBN _dbn;

    // Outflows
    Outflows _outflows;

    // Binsearcher, for searching bins
    Utils::BinSearcher _binSearcherX, _binSearcherY;

    EdgePair1D _xRange, _yRange;

    // Mapping from bin-searcher indices to bin indices (allowing gaps)
    std::vector<ssize_t> _indexes;

    // Numbers of edges on axes (necessary for bounds checking and indexing)
    size_t _nx, _ny;

    /// Whether modifying bin edges is permitted
    bool _locked;

    /// @}

  };

}

#endif
