import util

# cdef class Binned(AnalysisObject):
#     pass
# cdef class Fillable(AnalysisObject):
#     pass

cdef class Histo1D(AnalysisObject):
    """
    1D histogram, with distinction between bin areas and heights.

    Complete histogram binning is supported, including uniform/regular binning,
    variable-width binning, unbinned gaps in the covered range, and
    under/overflows. Rebinning by integer factors, or by explicit merging of
    contiguous bins is also supported.

    Rescaling of weights and/or the x axis is permitted in-place: the result is
    still a valid Histo1D. Binning-compatible 1D histograms may be divided,
    resulting in a Scatter2D since further fills would not be meaningful.

    Several sets of arguments are tried by the constructor in the
    following order.

    Histo1D(path="", title="").
      Construct a histogram with optional path and title but no bins.

    Histo1D(nbins, low, high, path="", title="")
      Construct a histogram with optional path and title, and nbins bins
      uniformly distributed between low and high.

    Histo1D(B, path="", title="").
      Construct a histogram with optional path and title, from an
      iterator of bins, B.
    """

    cdef inline c.Histo1D* h1ptr(self) except NULL:
        return <c.Histo1D*> self.ptr()

    def __init__(self, *args, **kwargs):
        util.try_loop([self.__init2, self.__init5, self.__init3], *args, **kwargs)

    def __init2(self, path="", title=""):
        path  = path.encode('utf-8')
        title = title.encode('utf-8')
        cutil.set_owned_ptr(self, new c.Histo1D(<string>path, <string>title))

    # TODO: Is Cython clever enough that we can make 3a and 3b versions and let it do the type inference?
    def __init3(self, bins_or_edges, path="", title=""):
        # TODO: Do this type-checking better
        cdef vector[double] edges
        try:
            path  = path.encode('utf-8')
            title = title.encode('utf-8')
            ## If float conversions work for all elements, it's a list of edges:
            edges = list(float(x) for x in bins_or_edges)
            cutil.set_owned_ptr(self, new c.Histo1D(edges, <string>path, <string>title))
        except:
            ## Assume it's a list of HistoBin1D
            bins = bins_or_edges
            self.__init2(path, title)
            self.addBins(bins)

    def __init5(self, nbins, low, high, path="", title=""):
        path  = path.encode('utf-8')
        title = title.encode('utf-8')
        cutil.set_owned_ptr(self, new c.Histo1D(nbins, low, high, <string>path, <string>title))


    def __len__(self):
        "Number of bins"
        return self.numBins()

    def __getitem__(self, i):
        "Direct access to bins"
        cdef size_t ii = cutil.pythonic_index(i, self.numBins())
        return cutil.new_borrowed_cls(HistoBin1D, & self.h1ptr().bin(ii), self)


    def __repr__(self):
        xmean = None
        if self.sumW() != 0:
            xmean = "%0.2e" % self.xMean()
        return "<%s '%s' %d bins, sumw=%0.2g, xmean=%s>" % \
               (self.__class__.__name__, self.path(),
                len(self.bins()), self.sumW(), xmean)


    def reset(self):
        """None -> None.
        Reset the histogram but leave the bin structure."""
        self.h1ptr().reset()


    def clone(self):
        """None -> Histo1D.
        Clone this Histo1D."""
        return cutil.new_owned_cls(Histo1D, self.h1ptr().newclone())


    def fill(self, x, weight=1.0, fraction=1.0):
        """(x,[w]) -> None.
        Fill with given x value and optional weight."""
        self.h1ptr().fill(x, weight, fraction)


    def fillBin(self, size_t ix, weight=1.0, fraction=1.0):
        """(ix,[w]) -> None.
        Fill bin ix and optional weight."""
        self.h1ptr().fillBin(ix, weight, fraction)


    def totalDbn(self):
        """None -> Dbn1D
        The Dbn1D representing the total distribution."""
        return cutil.new_borrowed_cls(Dbn1D, &self.h1ptr().totalDbn(), self)

    def underflow(self):
        """None -> Dbn1D
        The Dbn1D representing the underflow distribution."""
        return cutil.new_borrowed_cls(Dbn1D, &self.h1ptr().underflow(), self)

    def overflow(self):
        """None -> Dbn1D
        The Dbn1D representing the overflow distribution."""
        return cutil.new_borrowed_cls(Dbn1D, &self.h1ptr().overflow(), self)


    def integral(self, includeoverflows=True):
        """([bool]) -> float
        Histogram integral, optionally excluding the overflows."""
        return self.h1ptr().integral(includeoverflows)

    def integralRange(self, int ia, int ib):
        """(int, int) -> float
        Integral between bins ia..ib inclusive"""
        return self.h1ptr().integralRange(ia, ib)

    def integralTo(self, int ia, includeunderflow=True):
        """(int, [bool]) -> float
        Integral up to bin ia inclusive, optionally excluding the underflow"""
        return self.h1ptr().integralRange(ia, includeunderflow)


    def numEntries(self, includeoverflows=True):
        """([bool]) -> float
        Number of times this histogram was filled, optionally excluding the overflows."""
        return self.h1ptr().numEntries(includeoverflows)

    def effNumEntries(self, includeoverflows=True):
        """([bool]) -> float
        Effective number of times this histogram was filled, computed from weights, and optionally excluding the overflows."""
        return self.h1ptr().effNumEntries(includeoverflows)

    def sumW(self, includeoverflows=True):
        """([bool]) -> float
        Sum of weights filled into this histogram, optionally excluding the overflows."""
        return self.h1ptr().sumW(includeoverflows)

    def sumW2(self, includeoverflows=True):
        """([bool]) -> float
        Sum of weights filled into this histogram, optionally excluding the overflows."""
        return self.h1ptr().sumW2(includeoverflows)


    def xMean(self, includeoverflows=True):
        """([bool]) -> float
        Mean x of the histogram, optionally excluding the overflows."""
        return self.h1ptr().xMean(includeoverflows)

    def xVariance(self, includeoverflows=True):
        """([bool]) -> float
        Variance in x of the histogram, optionally excluding the overflows."""
        return self.h1ptr().xVariance(includeoverflows)

    def xStdDev(self, includeoverflows=True):
        """([bool]) -> float
        Standard deviation in x of the histogram, optionally excluding the overflows."""
        return self.h1ptr().xStdDev(includeoverflows)

    def xStdErr(self, includeoverflows=True):
        """([bool]) -> float
        Standard error on the mean x of the histogram, optionally excluding the overflows."""
        return self.h1ptr().xStdErr(includeoverflows)

    def xRMS(self, includeoverflows=True):
        """([bool]) -> float
        RMS in x of the histogram, optionally excluding the overflows."""
        return self.h1ptr().xRMS(includeoverflows)


    def scaleW(self, w):
        """ (float) -> None.
        Rescale the weights in this histogram by the factor w."""
        self.h1ptr().scaleW(w)

    def normalize(self, normto=1.0, includeoverflows=True):
        """ (float, bool) -> None.
        Normalize the histogram."""
        self.h1ptr().normalize(normto, includeoverflows)


    def xMin(self):
        """Low x edge of the histo."""
        return self.h1ptr().xMin()

    def xMax(self):
        """High x edge of the histo."""
        return self.h1ptr().xMax()

    def numBins(self):
        """() -> int
        Number of bins (not including overflows)."""
        return self.h1ptr().numBins()

    def numBinsX(self):
        """() -> int
        Number of x-axis bins (not including overflows)."""
        return self.h1ptr().numBinsX()

    def bins(self):
        """Access the ordered bins list."""
        return list(self)

    def bin(self, i):
        """Get the i'th bin (equivalent to bins[i]"""
        # cdef size_t ii = cutil.pythonic_index(i, self.h1ptr().numBins())
        return cutil.new_borrowed_cls(HistoBin1D, & self.h1ptr().bin(i), self)

    def binIndexAt(self, x):
        """Get the bin index containing position x"""
        return self.h1ptr().binIndexAt(x)

    def binAt(self, x):
        """Get the bin containing position x"""
        # TODO: what's the problem with this direct mapping? Produces compile error re. no default constructor...
        #return cutil.new_borrowed_cls(HistoBin1D, & self.h1ptr().binAt(x), self)
        # TODO: need out-of-range check to return None?
        return self.bin(self.binIndexAt(x))

    def addBin(self, low, high):
        """(low, high) -> None.
        Add a bin."""
        self.h1ptr().addBin(low, high)

    def addBins(self, edges_or_bins):
        """Add several bins."""
        # TODO: simplify / make consistent
        arg = list(edges_or_bins)
        util.try_loop([self.__addBins_edges,
                       self.__addBins_tuples,
                       self.__addBins_points,
                       self.__addBins_bins], arg)

    def __addBins_edges(self, edges):
        cdef vector[double] cedges
        for edge in edges:
            cedges.push_back(edge)
        if len(edges):
            self.h1ptr().addBins(cedges)

    def __addBins_bins(self, bins):
        self.__addBins_tuples([ b.xEdges for b in bins ])

    def __addBins_points(self, points):
        self.__addBins_tuples([ p.xWidth for p in points ])

    def __addBins_tuples(self, tuples):
        cdef double a, b
        for a, b in tuples:
            self.h1ptr().addBin(a, b)


    def mergeBins(self, ia, ib):
        """mergeBins(ia, ib) -> None.
        Merge bins from indices ia through ib."""
        self.h1ptr().mergeBins(ia, ib)

    def rebinBy(self, n, begin=0, end=None):
        """(n) -> None.
        Merge every group of n bins together (between begin and end, if specified)."""
        if end is None:
            end = self.numBins()
        self.h1ptr().rebinBy(int(n), begin, end)

    def rebinTo(self, edges):
        """([edges]) -> None.
        Merge bins to produce the given new edges... which must be a subset of the current ones."""
        self.h1ptr().rebinTo(edges)

    def rebin(self, arg, **kwargs):
        """(n) -> None or ([edges]) -> None
        Merge bins, like rebinBy if an int argument is given; like rebinTo if an iterable is given."""
        if hasattr(arg, "__iter__"):
            self.rebinTo(arg, **kwargs)
        else:
            self.rebinBy(arg, **kwargs)


    def mkScatter(self, usefocus=False, binwidthdiv=True, uflow_binwidth=-1, oflow_binwidth=-1):
        """None -> Scatter2D.
        Convert this Histo1D to a Scatter2D, with optional argument to control the x-positions
        of points within the bins, whether y represents sumW or density, and optional under-
        and overflow points."""
        cdef c.Scatter2D s2 = c.mkScatter_Histo1D(deref(self.h1ptr()), usefocus, binwidthdiv, uflow_binwidth, oflow_binwidth)
        return cutil.new_owned_cls(Scatter2D, s2.newclone())


    def toIntegral(self, efficiency=False, includeunderflow=True, includeoverflow=True):
        """None -> Scatter2D.

        Convert this Histo1D to a Scatter2D representing an integral (i.e. cumulative)
        histogram constructed from this differential one.

        The efficiency argument is used to construct an 'efficiency integral' histogram
        and the includeXXXflow bools determine whether under and overflows are included
        in computing the (efficiency) integral.
        """
        cdef c.Scatter2D s
        if not efficiency:
            s = c.Histo1D_toIntegral(deref(self.h1ptr()), includeunderflow)
        else:
            s = c.Histo1D_toIntegralEff(deref(self.h1ptr()), includeunderflow, includeoverflow)
        return cutil.new_owned_cls(Scatter2D, s.newclone())

    def divideBy(self, Histo1D h, efficiency=False):
        """Histo1D -> Scatter2D

        Divide this histogram by h, returning a Scatter2D. The optional 'efficiency'
        argument, if set True, will use a binomial efficiency treatment of the errors.
        """
        # if type(h) is not Histo1D:
        #     raise ValueError("Histograms must be of the same type to be divided")
        cdef c.Scatter2D s
        if not efficiency:
            s = c.Histo1D_div_Histo1D(deref(self.h1ptr()), deref(h.h1ptr()))
        else:
            s = c.Histo1D_eff_Histo1D(deref(self.h1ptr()), deref(h.h1ptr()))
        return cutil.new_owned_cls(Scatter2D, s.newclone())


    ## In-place special methods

    def __iadd__(Histo1D self, Histo1D other):
        c.Histo1D_iadd_Histo1D(self.h1ptr(), other.h1ptr())
        return self

    def __isub__(Histo1D self, Histo1D other):
        c.Histo1D_isub_Histo1D(self.h1ptr(), other.h1ptr())
        return self

    # def __imul__(Histo1D self, double x):
    #     c.Histo1D_imul_dbl(self.h1ptr(), x)
    #     return self

    # def __idiv__(Histo1D self, double x):
    #     c.Histo1D_idiv_dbl(self.h1ptr(), x)
    #     return self


    ## Unbound special methods

    # # TODO: only to bootstrap sum(), but doesn't work properly? Seems to treat *self* as the int...
    # def __radd__(Histo1D self, zero):
    #     #assert zero == 0
    #     print "FOO"
    #     return self.clone()

    def __add__(Histo1D self, Histo1D other):
        # print "BAR"
        h = Histo1D()
        cutil.set_owned_ptr(h, c.Histo1D_add_Histo1D(self.h1ptr(), other.h1ptr()))
        return h

    # TODO: Cython doesn't support type overloading for special functions?
    # def __add__(Histo1D self, int x):
    #     """
    #     Special operator support to allow use of sum(histos) which starts from 0.
    #     """
    #     assert(x == 0)
    #     return self

    # TODO: Cython doesn't support type overloading for special functions?
    # def __radd__(Histo1D self, int x):
    #     """
    #     Special operator support to allow use of sum(histos) which starts from 0.
    #     """
    #     assert(x == 0)
    #     return self

    def __sub__(Histo1D self, Histo1D other):
        h = Histo1D()
        cutil.set_owned_ptr(h, c.Histo1D_sub_Histo1D(self.h1ptr(), other.h1ptr()))
        return h

    # def __mul__(Histo1D self, double x):
    #     h = c.Histo1D_mul_dbl(self.h1ptr(), x)
    #     return h
    # def __rmul__(Histo1D self, double x):
    #     h = c.Histo1D_mul_dbl(self.h1ptr(), x)
    #     return h

    def __div__(Histo1D self, Histo1D other):
        return self.divideBy(other)

    def __truediv__(Histo1D self, Histo1D other):
        return self.divideBy(other)



    ## Functions for array-based plotting, chi2 calculations, etc.

    # def sumWs(self):
    #     """All sumWs of the histo."""
    #     return [b.sumW for b in self.bins]

    def _mknp(self, xs):
        try:
            import numpy
            return numpy.array(xs)
        except ImportError:
            return xs

    def xEdges(self):
        """All x edges of the histo."""
        return self._mknp(self.h1ptr().xEdges())

    def xWidths(self):
        """All x widths of the histo."""
        return self._mknp(self.h1ptr().xWidths())

    def xMins(self):
        """All x low edges of the histo."""
        return self._mknp([b.xMin() for b in self.bins()])

    def xMaxs(self):
        """All x high edges of the histo."""
        return self._mknp([b.xMax() for b in self.bins()])

    def xMids(self):
        """All x bin midpoints of the histo."""
        return self._mknp([b.xMid() for b in self.bins()])

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

    def heights(self):
        """All y heights of the histo."""
        return self._mknp([b.height() for b in self.bins()])

    def areas(self):
        """All areas of the histo."""
        return self._mknp([b.area() for b in self.bins()])

    def yVals(self, area=False):
        return self.areas() if area else self.heights()


    def heightErrs(self): #, asymm=False):
        """All height errors of the histo.

        TODO: asymm arg / heightErrsMinus/Plus?
        """
        return self._mknp([b.heightErr() for b in self.bins()])

    def areaErrs(self): #, asymm=False):
        """All area errors of the histo.

        TODO: asymm arg / areaErrsMinus/Plus?
        """
        # Use symmetrised errors by default, or return a list of (-,+) pairs if asymm is requested."""
        # if asymm:
        #    pass
        #else:
        return self._mknp([b.areaErr() for b in self.bins()])

    def relErrs(self): #, asymm=False):
        """All relative errors of the histo.

        TODO: asymm arg / areaErrsMinus/Plus?
        """
        return self._mknp([b.relErr() for b in self.bins()])

    def yErrs(self, area=False):
        return self.areaErrs() if area else self.heightErrs()


    def yMins(self, area=False):
        ys = self.yVals(area)
        es = self.yErrs(area)
        return self._mknp([y-e for (y,e) in zip(ys,es)])

    def yMaxs(self, area=False):
        ys = self.yVals(area)
        es = self.yErrs(area)
        return self._mknp([y+e for (y,e) in zip(ys,es)])

    def yMin(self):
        """Lowest x value."""
        return min(self.yMins())

    def yMax(self):
        """Highest y value."""
        return max(self.yMaxs())


## Convenience alias
H1D = Histo1D
