// -*- C++ -*-
//
// This file is part of YODA -- Yet more Objects for Data Analysis
// Copyright (C) 2008-2017 The YODA collaboration (see AUTHORS for details)
//
#ifndef YODA_BINSEARCHER_H
#define YODA_BINSEARCHER_H

#include "YODA/Utils/fastlog.h"
#include "YODA/Utils/MathUtils.h"
#include <cstdlib>
#include <cmath>
#include <vector>
#include <limits>
#include <memory>

namespace YODA {
  namespace Utils {


    const size_t SEARCH_SIZE = 16;
    const size_t BISECT_LINEAR_THRESHOLD = 32;


    /// @brief Bin estimator
    ///
    /// Base class for guessing the right bin index for a given value. The
    /// better the guess, the less time spent looking.
    struct Estimator {

      /// Virtual destructor needed for inheritance
      virtual ~Estimator() {}

      /// Return offset bin index estimate, with 0 = underflow and Nbins+1 = overflow
      size_t estindex(double x) const {
        const int i = _est(x);
        if (i < 0) return 0;
        const size_t i2 = (size_t) i;
        if (i2 >= _N) return _N+1;
        return i2 + 1;
      }

      /// Return offset bin index estimate, with 0 = underflow and Nbins+1 = overflow
      size_t operator() (double x) const {
        return estindex(x);
      }

    protected:

      /// Make an int-valued estimate of bin index
      /// @note No range checking or underflow offset
      virtual int _est(double x) const = 0;

      /// Number of bins
      size_t _N;
    };


    /// @brief Linear bin estimator
    ///
    /// This class handles guessing a index bin with a hypothesis of uniformly
    /// spaced bins on a linear scale.
    struct LinEstimator : public Estimator {

      /// Constructor
      LinEstimator(size_t nbins, double xlow, double xhigh) {
        _N = nbins;
        _c = xlow;
        _m = (double) nbins / (xhigh - xlow);
      }

      /// Copy constructor
      LinEstimator(const LinEstimator& other) {
        _N = other._N;
        _c = other._c;
        _m = other._m;
      }

      /// Call operator returns estimated bin index (offset so 0 == underflow)
      int _est(double x) const {
        return (int) floor(_m * (x - _c));
      }

    protected:
      double _c, _m;
    };


    /// @brief Logarithmic bin estimator
    ///
    /// This class handles guessing a bin index with a hypothesis of uniformly
    /// spaced bins on a logarithmic scale.
    ///
    /// @todo Make a generalised version of this with a transform function
    class LogEstimator : public Estimator {
    public:

      /// Constructor
      LogEstimator(size_t nbins, double xlow, double xhigh) {
        _N = nbins;
        _c = log2(xlow);
        _m = nbins / (log2(xhigh) - _c);
      }

      /// Copy constructor
      LogEstimator(const LogEstimator& other) {
        _N = other._N;
        _c = other._c;
        _m = other._m;
      }

      /// Call operator returns estimated bin index (offset so 0 == underflow)
      int _est(double x) const {
        return (int) floor(_m * (fastlog2(x) - _c));
      }

    protected:
      double _c, _m;
    };


    /// @brief Bin searcher
    ///
    /// @author David Mallows
    /// @author Andy Buckley
    ///
    /// Handles low-level bin lookups using a hybrid algorithm that is
    /// considerably faster for regular (logarithmic or linear) and near-regular
    /// binnings. Comparable performance for irregular binnings.
    ///
    /// The reason this works is that linear search is faster than bisection
    /// search up to about 32-64 elements. So we make a guess, and we then do a
    /// linear search. If that fails, then we bisect on the remainder,
    /// terminating once bisection search has got the range down to about 32. So
    /// we actually pay for the fanciness of predicting the bin out of speeding
    /// up the bisection search by finishing it with a linear search. So in most
    /// cases, we get constant-time lookups regardless of the space.
    ///
    class BinSearcher {
    public:

      /// Default constructor
      /// @todo What's the point? Remove?
      BinSearcher() {
        _est = std::make_shared<LinEstimator>(0, 0, 1);
      }

      /// Copy constructor
      BinSearcher(const BinSearcher& bs) {
        _est = bs._est;
        _edges = bs._edges;
      }

      // /// Explicit constructor, specifying the edges and estimation strategy
      // BinSearcher(const std::vector<double>& edges, bool log) {
      //   _updateEdges(edges);
      //   // Internally use a log or linear estimator as requested
      //   if (log) {
      //     _est.reset(new LogEstimator(edges.size()-1, edges.front(), edges.back()));
      //   } else {
      //     _est.reset(new LinEstimator(edges.size()-1, edges.front(), edges.back()));
      //   }
      // }

      /// Fully automatic constructor: give bin edges and it does the rest!
      BinSearcher(const std::vector<double>& edges) {
        _updateEdges(edges);

        if (edges.empty()) {
          _est = std::make_shared<LinEstimator>(0, 0, 1);
        } else if (edges.front() <= 0.0) {
          _est = std::make_shared<LinEstimator>(edges.size()-1, edges.front(), edges.back());
        } else {
          LinEstimator linEst(edges.size()-1, edges.front(), edges.back());
          LogEstimator logEst(edges.size()-1, edges.front(), edges.back());

          // Calculate mean index estimate deviations from the correct answers (for bin edges)
          double logsum = 0, linsum = 0;
          for (size_t i = 0; i < edges.size(); i++) {
            logsum += logEst(edges[i]) - i;
            linsum += linEst(edges[i]) - i;
          }
          const double log_avg = logsum / edges.size();
          const double lin_avg = linsum / edges.size();

          // This also implicitly works for NaN returned from the log There is a
          // subtle bug here if the if statement is the other way around, as
          // (nan < linsum) -> false always.  But (nan > linsum) -> false also.
          if (log_avg < lin_avg) { //< Use log estimator if its avg performance is better than lin
            _est = std::make_shared<LogEstimator>(logEst);
          } else { // Else use linear estimation
            _est = std::make_shared<LinEstimator>(linEst);
          }
        }
      }


      /// Look up a bin index
      /// @note Returned indices are offset by one, so 0 = underflow and Nbins+1 = overflow
      size_t index(double x) const {
        // Get initial estimate
        size_t index = std::min(_est->estindex(x),_edges.size()-1);
        // Return now if this is the correct bin
        if (x >= _edges[index] && x < _edges[index+1]) return index;

        // Otherwise refine the estimate, if x is not exactly on a bin edge
        if (x > _edges[index]) {
          const ssize_t newindex = _linsearch_forward(index, x, SEARCH_SIZE);
          index = (newindex > 0) ? newindex : _bisect(x, index, _edges.size()-1);
        } else if (x < _edges[index]) {
          const ssize_t newindex = _linsearch_backward(index, x, SEARCH_SIZE);
          index = (newindex > 0) ? newindex : _bisect(x, 0, index+1);
        }

        assert(x >= _edges[index] && (x < _edges[index+1] || std::isinf(x)));
        return index;
      }

      /// Look up an in-range bin index
      /// @note This returns a *normal* index starting with zero for the first in-range bin
      ssize_t index_inrange(double x) const {
        const size_t i = index(x);
        if (i == 0 || i == _edges.size()-1) return -1;
        return i;
      }


      /// Public access to the list of bin edges, including infinities at either end
      const std::vector<double>& edges() const { return _edges; }

      /// Public access to a bin edge value
      double edge(size_t i) const { return edges().at(i); }

      /// How many bin edges in this searcher?
      size_t size() const { return _edges.size(); }


      /// Check if two BinSearcher objects have the same edges
      bool same_edges(const BinSearcher& other) const {
        if (size() != other.size()) return false;
        for (size_t i = 1; i < size()-1; i++) {
          /// @todo Be careful about using fuzzyEquals... should be an exact comparison?
          if (!fuzzyEquals(edge(i), other.edge(i))) return false;
        }
        return true;
      }


      /// Find edges which are shared between BinSearcher objects, within numeric tolerance
      /// @note The return vector is sorted and includes -inf and inf
      std::vector<double> shared_edges(const BinSearcher& other) const {
        std::vector<double> rtn;
        rtn.push_back(-std::numeric_limits<double>::infinity());

        // Primarily loop over the smaller axis, since shared_edges \in {smaller}
        const int ndiff = size() - other.size();
        const BinSearcher& larger = (ndiff > 0) ? *this : other;
        const BinSearcher& smaller = (ndiff > 0) ? other : *this;
        size_t jmin = 1; //< current index in inner axis, to avoid unnecessary recomparisons (since vectors are sorted)
        for (size_t i = 1; i < smaller.size()-1; ++i) {
          const double x = smaller.edge(i);
          for (size_t j = jmin; j < larger.size()-1; ++j) {
            if (fuzzyEquals(x, larger.edge(j))) {
              rtn.push_back(x);
              jmin = j+1;
              break;
            }
          }
        }

        rtn.push_back(std::numeric_limits<double>::infinity());
        std::sort(rtn.begin(), rtn.end());
        return rtn;
      }


      // /// Check if two BinSearcher objects have compatible edges, i.e. one is a subset of the other
      // bool subset_edges(const BinSearcher& other) const {
      //   if (_edges.size() != other._edges.size()) return false;

      //   const int ndiff = numBins() - other.numBins();
      //   if (ndiff == 0) return same_edges(other);

      //   const BinSearcher& larger = (ndiff > 0) ? *this : other;
      //   const BinSearcher& smaller = (ndiff < 0) ? other : *this;
      //   size_t jmin = 1; //< current index in larger axis, to avoid unnecessary recomparisons
      //   for (size_t i = 1; i < smaller.size()-1; ++i) {
      //     bool found_match = false;
      //     for (size_t j = jmin; j < larger.size()-1; ++j) {
      //       if (fuzzyEquals(smaller.edge(i), larger.edge(j))) {
      //         found_match = true;
      //         jmin = j+1;
      //         break;
      //       }
      //     }
      //     if (!found_match) return false;
      //   }
      // }


    protected:

      /// Set the edges array and related member variables
      void _updateEdges(const std::vector<double>& edges) {
        // Array of in-range edges, plus underflow and overflow sentinels
        _edges.clear();
        _edges.resize(edges.size() + 2);

        // Copy vector with -+inf at ends
        _edges[0] = -std::numeric_limits<double>::infinity();
        for (size_t i = 0; i < edges.size(); i++) _edges[i+1] = edges[i];
        _edges[_edges.size()-1] = std::numeric_limits<double>::infinity();
      }


      /// @brief Linear search in the forward direction
      ///
      /// Return bin index or -1 if not found within linear search range. Assumes that edges[istart] <= x
      ssize_t _linsearch_forward(size_t istart, double x, size_t nmax) const {
        assert(x >= _edges[istart]); // assumption that x >= start is wrong
        for (size_t i = 0; i < nmax; i++) {
          const size_t j = istart + i + 1; // index of _next_ edge
          if (j > _edges.size()-1) return -1;
          if (x < _edges[j]) {
            assert(x >= _edges[j-1] && (x < _edges[j] || std::isinf(x)));
            return j-1; // note one more iteration needed if x is on an edge
          }
        }
        return -1;
      }

      /// @brief Linear search in the backward direction
      ///
      /// Return bin index or -1 if not found within linear search range. Assumes that edges[istart] > x
      ssize_t _linsearch_backward(size_t istart, double x, size_t nmax) const {
        assert(x < _edges[istart]); // assumption that x < start is wrong
        for (size_t i = 0; i < nmax; i++) {
          const int j = istart - i - 1; // index of _next_ edge (working backwards)
          if (j < 0) return -1;
          if (x >= _edges[j]) {
            assert(x >= _edges[j] && (x < _edges[j+1] || std::isinf(x)));
            return (ssize_t) j; // note one more iteration needed if x is on an edge
          }
        }
        return -1;
      }


      /// Truncated bisection search, adapted from C++ std lib implementation
      size_t _bisect(double x, size_t imin, size_t imax) const {
        size_t len = imax - imin;
        while (len >= BISECT_LINEAR_THRESHOLD) {
          const size_t half = len >> 1;
          const size_t imid = imin + half;
          if (x >= _edges[imid]) {
            if (x < _edges[imid+1]) return imid; // Might as well return directly if we get lucky!
            imin = imid;
          } else {
            imax = imid;
          }
          len = imax - imin;
        }
        assert(x >= _edges[imin] && (x < _edges[imax] || std::isinf(x)));
        return _linsearch_forward(imin, x, BISECT_LINEAR_THRESHOLD);
      }


    protected:

      /// Estimator object to be used for making fast bin index guesses
      std::shared_ptr<Estimator> _est;

      /// List of bin edges, including +- inf at either end
      std::vector<double> _edges;

    };


  }
}

#endif
