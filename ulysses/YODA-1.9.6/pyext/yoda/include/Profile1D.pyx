cimport util
cdef class Profile1D(AnalysisObject):
    """
    1D profile histogram, used to measure mean values of a y variable, binned in x.

    Complete histogram binning is supported, including uniform/regular binning,
    variable-width binning, unbinned gaps in the covered range, and
    under/overflows. Rebinning by integer factors, or by explicit merging of
    contiguous bins is also supported.

    Rescaling of weights and/or the x axis is permitted in-place: the result is
    still a valid Histo1D. Binning-compatible 1D histograms may be divided,
    resulting in a Scatter2D since further fills would not be meaningful.

    Several sets of arguments are tried by the constructor in the following
    order.

    Profile1D(path="", title="").
      Construct a histogram with optional path and title but no bins.

    Profile1D(nbins, low, high, path="", title="")
      Construct a histogram with optional path and title, and nbins bins
      uniformly distributed between low and high.

    Profile1D(B, path="", title="").
      Construct a histogram with optional path and title, from an
      iterator of bins, B.
    """

    cdef inline c.Profile1D* p1ptr(self) except NULL:
        return <c.Profile1D*> self.ptr()


    def __init__(self, *args, **kwargs):
        util.try_loop([self.__init2, self.__init3, self.__init5], *args, **kwargs)

    def __init2(self, path="", title=""):
        path  = path.encode('utf-8')
        title = title.encode('utf-8')
        cutil.set_owned_ptr(
            self, new c.Profile1D(<string>path, <string>title))

    # TODO: Is Cython clever enough that we can make 3a and 3b versions and let it do the type inference?
    def __init3(self, bins_or_edges, path="", title=""):
        # TODO: Do this type-checking better
        cdef vector[double] edges
        try:
            path  = path.encode('utf-8')
            title = title.encode('utf-8')
            ## If float conversions work for all elements, it's a list of edges:
            edges = list(float(x) for x in bins_or_edges)
            cutil.set_owned_ptr(self, new c.Profile1D(edges, <string>path, <string>title))
        except:
            ## Assume it's a list of HistoBin1D
            bins = bins_or_edges
            self.__init2(path, title)
            self.addBins(bins)

    def __init5(self, size_t nbins, double lower, double upper, path="", title=""):
        path  = path.encode('utf-8')
        title = title.encode('utf-8')
        cutil.set_owned_ptr(
            self, new c.Profile1D(nbins, lower, upper, <string>path, <string>title))


    def __len__(self):
        "Number of bins"
        return self.p1ptr().bins().size()

    def __getitem__(self, i):
        "Direct access to bins"
        cdef size_t ii = cutil.pythonic_index(i, self.p1ptr().bins().size())
        return cutil.new_borrowed_cls(ProfileBin1D, & self.p1ptr().bin(ii), self)


    def __repr__(self):
        return "<%s '%s' %d bins, sumw=%0.2g>" % \
               (self.__class__.__name__, self.path(),
                len(self.bins()), self.sumW())


    def reset(self):
        """None -> None.
        Reset the histogram but leave the bin structure."""
        self.p1ptr().reset()

    def clone(self):
        """None -> Profile1D.
        Clone this Profile1D."""
        return cutil.new_owned_cls(Profile1D, self.p1ptr().newclone())


    def fill(self, x, y, weight=1.0, fraction=1.0):
        """(x,y,[w]) -> None.
        Fill with given x & y values and optional weight."""
        self.p1ptr().fill(x, y, weight, fraction)

    def fillBin(self, size_t ix, double y, double weight=1.0, double fraction=1.0):
        """(ix,y,[w]) -> None.
        Fill bin ix with y value and optional weight."""
        self.p1ptr().fillBin(ix, y, weight, fraction)


    def totalDbn(self):
        """() -> Dbn2D
        The Dbn2D representing the total distribution."""
        return cutil.new_borrowed_cls(
            Dbn2D, &self.p1ptr().totalDbn(), self)

    def underflow(self):
        """() -> Dbn2D
        The Dbn2D representing the underflow distribution."""
        return cutil.new_borrowed_cls(
            Dbn2D, &self.p1ptr().underflow(), self)

    def overflow(self):
        """() -> Dbn2D
        The Dbn2D representing the overflow distribution."""
        return cutil.new_borrowed_cls(
            Dbn2D, &self.p1ptr().overflow(), self)


    def numEntries(self, includeoverflows=True):
        """([bool]) -> float
        Number of times this histogram was filled, optionally excluding the overflows."""
        return self.p1ptr().numEntries(includeoverflows)

    def effNumEntries(self, includeoverflows=True):
        """([bool]) -> float
        Effective number of times this histogram was filled, computed from weights and optionally excluding the overflows."""
        return self.p1ptr().effNumEntries(includeoverflows)

    def sumW(self, includeoverflows=True):
        """([bool]) -> float
        Sum of weights filled into this histogram."""
        return self.p1ptr().sumW(includeoverflows)

    def sumW2(self, includeoverflows=True):
        """([bool]) -> float
        Sum of weights filled into this histogram."""
        return self.p1ptr().sumW2(includeoverflows)


    def xMean(self, includeoverflows=True):
        """([bool]) -> float
        Mean x of the histogram, optionally excluding the overflows."""
        return self.p1ptr().xMean(includeoverflows)

    def xVariance(self, includeoverflows=True):
        """([bool]) -> float
        Variance in x of the histogram, optionally excluding the overflows."""
        return self.p1ptr().xVariance(includeoverflows)

    def xStdDev(self, includeoverflows=True):
        """([bool]) -> float
        Standard deviation in x of the histogram, optionally excluding the overflows."""
        return self.p1ptr().xStdDev(includeoverflows)

    def xStdErr(self, includeoverflows=True):
        """([bool]) -> float
        Standard error on the mean x of the histogram, optionally excluding the overflows."""
        return self.p1ptr().xStdErr(includeoverflows)

    def xRMS(self, includeoverflows=True):
        """([bool]) -> float
        RMS in x of the histogram, optionally excluding the overflows."""
        return self.p1ptr().xRMS(includeoverflows)


    def scaleW(self, double w):
        """(float) -> None.
        Rescale the weights in this histogram by the factor w."""
        self.p1ptr().scaleW(w)

    def scaleY(self, double f):
        """(float) -> None.
        Scale the y-direction (profiled value) in this histogram by the factor f."""
        self.p1ptr().scaleY(f)


    def xMin(self):
        """Low x edge of the histo."""
        return self.p1ptr().xMin()

    def xMax(self):
        """High x edge of the histo."""
        return self.p1ptr().xMax()

    def numBins(self):
        """() -> int
        Number of bins (not including overflows)."""
        return self.p1ptr().numBins()

    def numBinsX(self):
        """() -> int
        Number of bins on the x-axis (not including overflows)."""
        return self.p1ptr().numBinsX()

    def bins(self):
        """Access the ordered bins list."""
        return list(self)

    def bin(self, i):
        """Get the i'th bin"""
        # cdef size_t ii = cutil.pythonic_index(i, self.numBins())
        return cutil.new_borrowed_cls(ProfileBin1D, & self.p1ptr().bin(i), self)

    def binIndexAt(self, x):
        """Get the bin index containing position x"""
        return self.p1ptr().binIndexAt(x)

    def binAt(self, x):
        """Get the bin containing position x"""
        # TODO: what's the problem with this direct mapping? Produces compile error re. no default constructor...
        # return cutil.new_borrowed_cls(ProfileBin1D, & self.p1ptr().binAt(x), self)
        # TODO: need out-of-range check to return None?
        return self.bin(self.binIndexAt(x))

    def addBin(self, low, high):
        """Add a bin."""
        self.p1ptr().addBin(low, high)
        return self

    def addBins(self, edges):
        """Add several bins."""
        # TODO: simplify / make consistent
        cdef vector[double] cedges
        for i in edges:
            cedges.push_back(i)
        self.p1ptr().addBins(cedges)
        return self


    def mergeBins(self, a, b):
        """mergeBins(ia, ib) -> None.
        Merge bins from indices ia through ib."""
        self.p1ptr().mergeBins(a, b)

    def rebinBy(self, n, begin=0, end=None):
        """(n) -> None.
        Merge every group of n bins together (between begin and end, if specified)."""
        if end is None:
            end = self.numBins()
        self.p1ptr().rebinBy(int(n), begin, end)

    def rebinTo(self, edges):
        """([edges]) -> None.
        Merge bins to produce the given new edges... which must be a subset of the current ones."""
        self.p1ptr().rebinTo(edges)

    def rebin(self, arg, **kwargs):
        """(n) -> None or ([edges]) -> None
        Merge bins, like rebinBy if an int argument is given; like rebinTo if an iterable is given."""
        if hasattr(arg, "__iter__"):
            self.rebinTo(arg, **kwargs)
        else:
            self.rebinBy(arg, **kwargs)


    def mkScatter(self, usefocus=False, usestddev=False, uflow_binwidth=-1, oflow_binwidth=-1):
        """None -> Scatter2D.
        Convert this Profile1D to a Scatter2D, with y representing
        mean bin y values and their standard errors (or std deviations if usestddev=True).
        The remaining optional parameters allow under- and overflow points to be created."""
        cdef c.Scatter2D s2 = c.mkScatter_Profile1D(deref(self.p1ptr()), usefocus, usestddev, uflow_binwidth, oflow_binwidth)
        return cutil.new_owned_cls(Scatter2D, s2.newclone())


    def divideBy(self, Profile1D h):
        cdef c.Scatter2D s = c.Profile1D_div_Profile1D(deref(self.p1ptr()), deref(h.p1ptr()))
        return cutil.new_owned_cls(Scatter2D, s.newclone())


    def __iadd__(Profile1D self, Profile1D other):
        c.Profile1D_iadd_Profile1D(self.p1ptr(), other.p1ptr())
        return self
    def __isub__(Profile1D self, Profile1D other):
        c.Profile1D_isub_Profile1D(self.p1ptr(), other.p1ptr())
        return self

    def __add__(Profile1D self, Profile1D other):
        h = Profile1D()
        cutil.set_owned_ptr(h, c.Profile1D_add_Profile1D(self.p1ptr(), other.p1ptr()))
        return h
    def __sub__(Profile1D self, Profile1D other):
        h = Profile1D()
        cutil.set_owned_ptr(h, c.Profile1D_sub_Profile1D(self.p1ptr(), other.p1ptr()))
        return h

    def __div__(Profile1D self, Profile1D other):
        return self.divideBy(other)

    def __truediv__(Profile1D self, Profile1D other):
        return self.divideBy(other)


    ## Functions for array-based plotting, chi2 calculations, etc.

    # def sumWs(self):
    #     """All sumWs of the histo."""
    #     return [b.sumW() for b in self.bins()]

    # TODO: xyVals,Errs properties should be in a common Drawable2D (?) type (hmm, need a consistent nD convention...)
    # TODO: x bin properties should be in a common Binned1D type
    # TODO: add "useoverflows" optional args, move most into C++

    def _mknp(self, xs):
        try:
            import numpy
            return numpy.array(xs)
        except ImportError:
            return xs


    ## Geometric properties in x

    def xEdges(self):
        """All x edges of the histo."""
        return self._mknp(self.p1ptr().xEdges())

    def xWidths(self):
        """All x widths of the histo."""
        return self._mknp(self.p1ptr().xWidths())

    def xMins(self):
        """All x low edges of the histo."""
        return self._mknp([b.xMin() for b in self.bins()])

    def xMaxs(self):
        """All x high edges of the histo."""
        return self._mknp([b.xMax() for b in self.bins()])

    def xMids(self):
        """All x bin midpoints of the histo."""
        return self._mknp([b.xMid() for b in self.bins()])


    ## Filling properties in x

    def xFoci(self):
        """All x bin foci of the histo."""
        return self._mknp([b.xFocus() for b in self.bins()])

    def xVals(self, foci=False):
        return self.xFoci() if foci else self.xMids()

    def xErrs(self, foci=False):
        if foci:
            return [(b.xFocus()-b.xMin(), b.xMax()-b.xFocus()) for b in self.bins()]
        else:
            return [(b.xMid()-b.xMin(), b.xMax()-b.xMid()) for b in self.bins()]

    def xMin(self):
        """Lowest x value."""
        return min(self.xMins())

    def xMax(self):
        """Highest x value."""
        return max(self.xMaxs())


    def sumWs(self):
        """All sumW values of the histo."""
        rtn = self._mknp([b.sumW() for b in self.bins()])
        return rtn

    def yMeans(self):
        """All y heights y means."""
        return self._mknp([b.mean() for b in self.bins()])

    def yVals(self):
        return self.yMeans()


    def yStdErrs(self):
        """All standard errors on the y means."""
        return self._mknp([b.stdErr() for b in self.bins()])

    def yStdDevs(self):
        """All standard deviations of the y distributions."""
        return self._mknp([b.stdDev() for b in self.bins()])

    def yErrs(self, sd=False):
        return self.yStdDevs() if sd else self.yStdErrs()


    def yMins(self, sd=False):
        ys = self.yVals()
        es = self.yErrs(sd)
        return self._mknp([y-e for (y,e) in zip(ys,es)])

    def yMaxs(self, sd=False):
        ys = self.yVals()
        es = self.yErrs(sd)
        return self._mknp([y+e for (y,e) in zip(ys,es)])

    def yMin(self, sd=False):
        """Lowest y value."""
        return min(self.yMins(sd))

    def yMax(self, sd=False):
        """Highest y value."""
        return max(self.yMaxs(sd))


## Convenience alias
P1D = Profile1D
