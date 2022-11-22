cimport util
cdef class Histo2D(AnalysisObject):
    """
    2D histogram.

    Complete histogramming is supported, including uniform/regular binning,
    variable-width bininng, unbinned gaps in the covered range, and outflows
    (under/overflows around all edges and corners).

    Rebinning by integer factors, or by explicit merging of contiguous bins is
    also supported, but in development.

    Rescaling of weights and/or the x axis is permitted in-place: the
    result is still a valid Histo2D. Binning-compatible 2D histograms
    may be divided, resulting in a Scatter3D since further fills would
    not be meaningful.

    Several sets of arguments are tried by the constructor in the
    following order.

    Histo2D(path="", title="").
      Construct a histogram with optional path and title but no bins.

    Histo2D(nxbins, xlow, xhigh, nybins, ylow, yhigh, path="", title="").
      Construct a histogram with nxbins on the x axis and nybins on the y
      axis, distributed linearly between the respective low--high limits.
    """

    cdef inline c.Histo2D* h2ptr(self) except NULL:
        return <c.Histo2D*> self.ptr()


    def __init__(self, *args, **kwargs):
        util.try_loop([self.__init2, self.__init4, self.__init8], *args, **kwargs)

    def __init2(Histo2D self, path="", title=""):
        path  = path.encode('utf-8')
        title = title.encode('utf-8')
        cutil.set_owned_ptr(self, new c.Histo2D(<string>path, <string>title))

    def __init4(Histo2D self, xedges,  yedges,  path="", title=""):
        path  = path.encode('utf-8')
        title = title.encode('utf-8')
        # TODO: Do some type-checking and allow iterables of HistoBin2D as well?
        cutil.set_owned_ptr(self, new c.Histo2D(xedges, yedges, <string>path, <string>title))

    def __init8(Histo2D self, nxbins, xlow, xhigh,  nybins, ylow, yhigh,  path="", title=""):
        path  = path.encode('utf-8')
        title = title.encode('utf-8')
        cutil.set_owned_ptr(self, new c.Histo2D(nxbins, xlow, xhigh,  nybins, ylow, yhigh,  <string>path, <string>title))


    def __len__(self):
        "Number of bins"
        return self.numBins()

    def __getitem__(self, py_ix):
        "Direct access to bins"
        cdef size_t i = cutil.pythonic_index(py_ix, self.numBins())
        return cutil.new_borrowed_cls(HistoBin2D, & self.h2ptr().bins().at(i), self)


    def __repr__(self):
        return "<%s '%s' %d bins, sumw=%.2g>" % (self.__class__.__name__, self.path(), len(self.bins()), self.sumW())


    def reset(self):
        """None -> None.
        Reset the histogram but leave the bin structure."""
        self.h2ptr().reset()

    def clone(self):
        """None -> Histo2D.
        Clone this Profile2D."""
        return cutil.new_owned_cls(Histo2D, self.h2ptr().newclone())


    def fill(self, double x, double y, weight=1.0, fraction=1.0):
        """(x,y,[w]) -> None.
        Fill with given x,y values and optional weight."""
        self.h2ptr().fill(x, y, weight, fraction)

    def fillBin(self, size_t i, weight=1.0, fraction=1.0):
        """(i,[w]) -> None.
        Fill bin i and optional weight."""
        self.h2ptr().fillBin(i, weight, fraction)


    #@property
    def totalDbn(self):
        """() -> Dbn2D
        The Dbn2D representing the total distribution."""
        return cutil.new_borrowed_cls(Dbn2D, &self.h2ptr().totalDbn(), self)

    # TODO: reinstate
    # def outflow(self, ix, iy):
    #     """(ix,iy) -> Dbn2D
    #     The Dbn2D representing the ix,iy outflow distribution."""
    #     return cutil.new_borrowed_cls(Dbn2D, &self.h2ptr().outflow(ix, iy), self)


    def integral(self, includeoverflows=True):
        """([bool]) -> float
        Histogram integral, optionally excluding the overflows."""
        return self.h2ptr().integral(includeoverflows)


    def numEntries(self, includeoverflows=True):
        """([bool]) -> float
        Number of times this histogram was filled, optionally excluding overflows."""
        return self.h2ptr().numEntries(includeoverflows)

    def effNumEntries(self, includeoverflows=True):
        """([bool]) -> float
        Effective number of times this histogram was filled, computed from weights and optionally excluding overflows."""
        return self.h2ptr().effNumEntries(includeoverflows)

    def sumW(self, includeoverflows=True):
        """([bool]) -> float
        Sum of weights filled into this histogram."""
        return self.h2ptr().sumW(includeoverflows)

    def sumW2(self, includeoverflows=True):
        """([bool]) -> float
        Sum of squared weights filled into this histogram."""
        return self.h2ptr().sumW2(includeoverflows)


    def xMean(self, includeoverflows=True):
        """([bool]) -> float
        Mean x of the histogram, optionally excluding the overflows."""
        return self.h2ptr().xMean(includeoverflows)

    def yMean(self, includeoverflows=True):
        """([bool]) -> float
        Mean y of the histogram, optionally excluding the overflows."""
        return self.h2ptr().yMean(includeoverflows)

    def xyMean(self, includeoverflows=True):
        """([bool]) -> (float,float)
        Mean (x,y) of the histogram, optionally excluding the overflows."""
        return util.XY(self.xMean(includeoverflows), self.yMean(includeoverflows))


    def xVariance(self, includeoverflows=True):
        """([bool]) -> float
        Variance in x of the histogram, optionally excluding the overflows."""
        return self.h2ptr().xVariance(includeoverflows)

    def yVariance(self, includeoverflows=True):
        """([bool]) -> float
        Variance in y of the histogram, optionally excluding the overflows."""
        return self.h2ptr().yVariance(includeoverflows)

    def xyVariance(self, includeoverflows=True):
        """([bool]) -> (float,float)
        Variances in (x,y) of the histogram, optionally excluding the overflows."""
        return util.XY(self.xVariance(includeoverflows), self.yVariance(includeoverflows))


    def xStdDev(self, includeoverflows=True):
        """([bool]) -> float
        Standard deviation in x of the histogram, optionally excluding the overflows."""
        return self.h2ptr().xStdDev(includeoverflows)

    def yStdDev(self, includeoverflows=True):
        """([bool]) -> float
        Standard deviation in y of the histogram, optionally excluding the overflows."""
        return self.h2ptr().yStdDev(includeoverflows)

    def xyStdDev(self, includeoverflows=True):
        """([bool]) -> (float,float)
        Standard deviations in (x,y) of the histogram, optionally excluding the overflows."""
        return util.XY(self.xStdDev(includeoverflows), self.yStdDev(includeoverflows))


    def xStdErr(self, includeoverflows=True):
        """([bool]) -> float
        Standard error on the mean x of the histogram, optionally excluding the overflows."""
        return self.h2ptr().xStdErr(includeoverflows)

    def yStdErr(self, includeoverflows=True):
        """([bool]) -> float
        Standard error on the mean y of the histogram, optionally excluding the overflows."""
        return self.h2ptr().yStdErr(includeoverflows)

    def xyStdErr(self, includeoverflows=True):
        """([bool]) -> (float,float)
        Standard errors on the mean (x,y) of the histogram, optionally excluding the overflows."""
        return util.XY(self.xStdErr(includeoverflows), self.yStdErr(includeoverflows))


    def xRMS(self, includeoverflows=True):
        """([bool]) -> float
        RMS in x of the histogram, optionally excluding the overflows."""
        return self.h2ptr().xRMS(includeoverflows)

    def yRMS(self, includeoverflows=True):
        """([bool]) -> float
        RMS in y of the histogram, optionally excluding the overflows."""
        return self.h2ptr().yRMS(includeoverflows)

    def xyRMS(self, includeoverflows=True):
        """([bool]) -> (float,float)
        RMS in (x,y) of the histogram, optionally excluding the overflows."""
        return util.XY(self.xRMS(includeoverflows), self.yRMS(includeoverflows))


    def scaleW(self, w):
        """(float) -> None.
        Rescale the weights in this histogram by the factor w."""
        self.h2ptr().scaleW(w)

    def normalize(self, double normto=1.0, bint includeoverflows=True):
        """(float, bool) -> None.
        Normalize the histogram."""
        self.h2ptr().normalize(normto, includeoverflows)


    def xMin(self):
        """Low x edge of the histo."""
        return self.h2ptr().xMin()

    def xMax(self):
        """High x edge of the histo."""
        return self.h2ptr().xMax()

    def yMin(self):
        """Low y edge of the histo."""
        return self.h2ptr().yMin()

    def yMax(self):
        """High y edge of the histo."""
        return self.h2ptr().yMax()


    def numBins(self):
        """() -> int
        Number of bins (not including overflows)."""
        return self.h2ptr().numBins()

    def numBinsX(self):
        """() -> int
        Number of bins (edges) along the x axis."""
        return self.h2ptr().numBinsX()

    def numBinsY(self):
        """() -> int
        Number of bins (edges) along the y axis."""
        return self.h2ptr().numBinsY()


    def bins(self):
        """Access the ordered bins list."""
        return [self.bin(i) for i in range( self.h2ptr().numBins())]

    def bin(self, i):
        """Get the i'th bin"""
        # cdef size_t ii = cutil.pythonic_index(i, self.h2ptr().numBins())
        return cutil.new_borrowed_cls(HistoBin2D, & self.h2ptr().bin(i), self)

    # TODO: it's more intuitive to have an index for each axis
    # def bin(self, i, j):
    #     """Get the (i,j)'th bin"""
    #     # cdef size_t ii = cutil.pythonic_index(i, self.h2ptr().numBins())
    #     # cdef size_t jj = cutil.pythonic_index(j, self.h2ptr().numBins())
    #     return cutil.new_borrowed_cls(HistoBin2D, & self.h2ptr().bin(i,j), self)

    def binIndexAt(self, x, y):
        """Get the bin index pair containing position (x,y)"""
        return self.h2ptr().binIndexAt(x, y)

    def binAt(self, x, y):
        """Get the bin containing position (x,y)"""
        # TODO: what's the problem with this direct mapping? Produces compile error re. no default constructor...
        #return cutil.new_borrowed_cls(HistoBin2D, & self.h2ptr().binAt(x,y), self)
        # TODO: need out-of-range check to return None?
        return self.bin(self.binIndexAt(x,y))


    def addBin(self, xlow, xhigh, ylow, yhigh):
        """Add a bin."""
        self.h2ptr().addBin(pair[double, double](xlow, xhigh),
                               pair[double, double](ylow, yhigh))
        return self

    def addBins(self, bins):
        """Add several bins."""
        cdef vector[c.HistoBin2D] cbins
        cdef c.HistoBin2D *cbin
        for elem in bins:
            if type(elem) is tuple:
                xlow, xhigh, ylow, yhigh = elem
                # Allocating on the stack requires a default constructor [1]
                # which HistoBin2D doesn't have.
                #
                # [1] https://github.com/cython/cython/wiki/WrappingCPlusPlus#declare-a-var-with-the-wrapped-c-class
                cbin = new c.HistoBin2D(xlow, xhigh, ylow, yhigh)
                # Assuming this doesn't throw exceptions, or we are going to
                # leak memory
                cbins.push_back(deref(cbin))
                del cbin
            elif type(elem) is HistoBin2D:
                cbin = (<HistoBin2D>elem).hb2ptr()
                cbins.push_back(deref(cbin))
            else:
                raise ValueError(
                    "Unsupported bin of type {}".format(type(elem))
                )
        self.h2ptr().addBins(cbins)

    # def mergeBins(self, size_t a, size_t b):
    #     self.h2ptr().mergeBins(a, b)

    # def rebin(self, int n):
    #     self.h2ptr().rebin(n)


    def mkScatter(self, usefocus=False, binareadiv=True):
        """None -> Scatter3D.
        Convert this Histo2D to a Scatter3D, with optional argument to control the positions
        of points within the bins, and whether z represents sumW or density."""
        cdef c.Scatter3D s3 = c.mkScatter_Histo2D(deref(self.h2ptr()), usefocus, binareadiv)
        return cutil.new_owned_cls(Scatter3D, s3.newclone())

    def divideBy(self, Histo2D h, efficiency=False):
        """Histo2D -> Scatter3D

        Divide this histogram by Histo2D h, returning a Scatter3D. The optional 'efficiency'
        argument, if set True, will use a binomial efficiency treatment of the errors.
        """
        # if type(h) is not Histo2D:
        #     raise ValueError("Histograms must be of the same type to be divided")
        # TODO: allow dividing profiles by histos, etc.? (But then what do the errors mean? Add in quad?)
        cdef c.Scatter3D s
        if not efficiency:
            s = c.Histo2D_div_Histo2D(deref(self.h2ptr()), deref(h.h2ptr()))
        else:
            s = c.Histo2D_eff_Histo2D(deref(self.h2ptr()), deref(h.h2ptr()))
        return cutil.new_owned_cls(Scatter3D, s.newclone())



    def __iadd__(Histo2D self, Histo2D other):
        c.Histo2D_iadd_Histo2D(self.h2ptr(), other.h2ptr())
        return self
    def __isub__(Histo2D self, Histo2D other):
        c.Histo2D_isub_Histo2D(self.h2ptr(), other.h2ptr())
        return self

    def __add__(Histo2D self, Histo2D other):
        h = Histo2D()
        cutil.set_owned_ptr(h, c.Histo2D_add_Histo2D(self.h2ptr(), other.h2ptr()))
        return h
    def __sub__(Histo2D self, Histo2D other):
        h = Histo2D()
        cutil.set_owned_ptr(h, c.Histo2D_sub_Histo2D(self.h2ptr(), other.h2ptr()))
        return h

    def __div__(Histo2D self, Histo2D other):
        return self.divideBy(other)

    def __truediv__(Histo2D self, Histo2D other):
        return self.divideBy(other)


    ## Functions for array-based plotting, chi2 calculations, etc.
    #
    # TODO: add "useoverflows" optional args, move most into C++

    def _mknp(self, xs):
        try:
            import numpy
            return numpy.array(xs)
        except ImportError:
            return xs


    ## Geometric properties in x

    def xEdges(self):
        """Unique x edges of the histo."""
        return self._mknp(self.h2ptr().xEdges())

    def xWidths(self):
        """All x widths of the histo."""
        return self._mknp(self.h2ptr().xWidths())

    def xMin(self):
        """Lowest x value."""
        return self.xEdges()[0]

    def xMax(self):
        """Highest x value."""
        return self.xEdges()[-1]

    def xMins(self, unique=True, asgrid=False):
        """Unique/all x low edges of the histo."""
        if unique:
            rtn = self.xEdges()[:-1]
        else:
            rtn = self._mknp([b.xMin() for b in self.bins()])
        if asgrid:
            rtn = rtn.reshape(self.numBinsX(), self.numBinsY())
        return rtn

    def xMaxs(self, unique=True, asgrid=False):
        """Unique/all x high edges of the histo."""
        if unique:
            rtn = self.xEdges()[1:]
        else:
            rtn = self._mknp([b.xMax() for b in self.bins()])
        if asgrid:
            rtn = rtn.reshape(self.numBinsX(), self.numBinsY())
        return rtn

    def xMids(self, unique=True, asgrid=False):
        """Unique/all x bin midpoints of the histo."""
        if unique:
            rtn = (self.xMins() + self.xMaxs())/ 2.
        else:
            rtn = self._mknp([b.xMid() for b in self.bins()])
        if asgrid:
            rtn = rtn.reshape(self.numBinsX(), self.numBinsY())
        return rtn


    ## Filling properties in x

    def xFoci(self, asgrid=False):
        """All x bin foci of the histo."""
        rtn = self._mknp([b.xFocus() for b in self.bins()])
        if asgrid:
            rtn = rtn.reshape(self.numBinsX(), self.numBinsY())
        return rtn

    def xVals(self, foci=False, asgrid=False):
        """All x values of the histo."""
        return self.xFoci(asgrid) if foci else self.xMids(False, asgrid)

    def xErrs(self, foci=False, asgrid=False):
        """All x errors of the histo."""
        # TODO: rewrite using xFoci/xMids/xMins/xMaxs
        if foci:
            rtn = [(b.xFocus()-b.xMin(), b.xMax()-b.xFocus()) for b in self.bins()]
        else:
            rtn = [(b.xMid()-b.xMin(), b.xMax()-b.xMid()) for b in self.bins()]
        if asgrid:
            rtn = rtn.reshape(self.numBinsX(), self.numBinsY())
        return rtn


    ## Geometric properties in y

    def yEdges(self):
        """Unique y edges of the histo."""
        return self._mknp(self.h2ptr().yEdges())

    def yWidths(self):
        """All y widths of the histo."""
        return self._mknp(self.h2ptr().yWidths())

    def yMin(self):
        """Lowest y value."""
        return self.yEdges()[0]

    def yMax(self):
        """Highest y value."""
        return self.yEdges()[-1]

    def yMins(self, unique=True, asgrid=False):
        """Unique/all y low edges of the histo."""
        if unique:
            rtn = self.yEdges()[:-1]
        else:
            rtn = self._mknp([b.yMin() for b in self.bins()])
        if asgrid:
            rtn = rtn.reshape(self.numBinsX(), self.numBinsY())
        return rtn

    def yMaxs(self, unique=True, asgrid=False):
        """Unique/all y high edges of the histo."""
        if unique:
            rtn = self.yEdges()[1:]
        else:
            rtn = self._mknp([b.yMax() for b in self.bins()])
        if asgrid:
            rtn = rtn.reshape(self.numBinsX(), self.numBinsY())
        return rtn

    def yMids(self, unique=True, asgrid=False):
        """Unique/all y bin midpoints of the histo."""
        if unique:
            rtn = (self.yMins() + self.yMaxs())/ 2.
        else:
            rtn = self._mknp([b.yMid() for b in self.bins()])
        if asgrid:
            rtn = rtn.reshape(self.numBinsX(), self.numBinsY())
        return rtn



    ## Filling properties in y

    def yFoci(self, asgrid=False):
        """All y bin foci of the histo."""
        rtn = self._mknp([b.yFocus() for b in self.bins()])
        if asgrid:
            rtn = rtn.reshape(self.numBinsX(), self.numBinsY())
        return rtn

    def yVals(self, foci=False, asgrid=False):
        """All y values of the histo."""
        return self.yFoci(asgrid) if foci else self.yMids(False, asgrid)

    def yErrs(self, foci=False, asgrid=False):
        """All y errors of the histo."""
        # TODO: rewrite using yFoci/yMids/yMins/yMaxs
        if foci:
            rtn = [(b.yFocus()-b.yMin(), b.yMax()-b.yFocus()) for b in self.bins()]
        else:
            rtn = [(b.yMid()-b.yMin(), b.yMax()-b.yMid()) for b in self.bins()]
        if asgrid:
            rtn = rtn.reshape(self.numBinsX(), self.numBinsY())
        return rtn


    ## Filling properties in general

    def sumWs(self, asgrid=False):
        """All sumW values of the histo."""
        rtn = self._mknp([b.sumW() for b in self.bins()])
        if asgrid:
            rtn = rtn.reshape(self.numBinsX(), self.numBinsY())
        return rtn

    def heights(self, asgrid=False):
        """All y heights of the histo."""
        rtn = self._mknp([b.height() for b in self.bins()])
        if asgrid:
            rtn = rtn.reshape(self.numBinsX(), self.numBinsY())
        return rtn

    def volumes(self, asgrid=False):
        """All volumes of the histo."""
        rtn = self._mknp([b.volume() for b in self.bins()])
        if asgrid:
            rtn = rtn.reshape(self.numBinsX(), self.numBinsY())
        return rtn

    def zVals(self, vol=False, asgrid=False):
        return self.volumes(asgrid) if vol else self.heights(asgrid)

    def heightErrs(self, asymm=False, asgrid=False):
        """All height errors of the histo.

        TODO: asymm arg / heightErrsMinus/Plus
        """
        rtn = self._mknp([b.heightErr() for b in self.bins()])
        if asgrid:
            rtn = rtn.reshape(self.numBinsX(), self.numBinsY())
        return rtn

    def volumeErrs(self, asymm=False, asgrid=False):
        """All volume errors of the histo.

        TODO: asymm arg / areaErrsMinus/Plus?
        """
        # Use symmetrised errors by default, or return a list of (-,+) pairs if asymm is requested."""
        # if asymm:
        #    pass
        #else:
        rtn = self._mknp([b.volumeErr() for b in self.bins()])
        if asgrid:
            rtn = rtn.reshape(self.numBinsX(), self.numBinsY())
        return rtn

    def relErrs(self, asymm=False, asgrid=False):
        """All relative errors of the histo.

        TODO: asymm arg / areaErrsMinus/Plus?
        """
        rtn = self._mknp([b.relErr() for b in self.bins()])
        if asgrid:
            rtn = rtn.reshape(self.numBinsX(), self.numBinsY())
        return rtn

    def zErrs(self, vol=False, asgrid=False):
        return self.volErrs(asgrid) if vol else self.heightErrs(asgrid)


    def zMins(self, area=False, asgrid=False):
        zs = self.zVals(area)
        es = self.zErrs(area)
        rtn = self._mknp([z-e for (z,e) in zip(zs,es)])
        if asgrid:
            rtn = rtn.reshape(self.numBinsX(), self.numBinsY())
        return rtn

    def zMaxs(self, area=False, asgrid=False):
        zs = self.zVals(area)
        es = self.zErrs(area)
        rtn = self._mknp([z+e for (z,e) in zip(zs,es)])
        if asgrid:
            rtn = rtn.reshape(self.numBinsX(), self.numBinsY())
        return rtn

    def zMin(self, area=False):
        """Lowest z value."""
        return min(self.zMins(area))

    def zMax(self, area=False):
        """Highest z value."""
        return max(self.zMaxs(area))


## Convenience alias
H2D = Histo2D
