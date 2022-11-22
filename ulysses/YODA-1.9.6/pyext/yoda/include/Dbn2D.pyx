cimport util
cdef class Dbn2D(util.Base):
    """
    A 2D distribution 'counter', used and exposed by 2D histograms and
    1D profiles and their bins.

    TODO: also provide normal scalar access to quantities like xRMS
    """

    cdef c.Dbn2D* d2ptr(self) except NULL:
        return <c.Dbn2D *> self.ptr()

    # TODO: remove!
    cdef c.Dbn2D* _Dbn2D(self) except NULL:
        return <c.Dbn2D *> self.ptr()

    def __dealloc__(self):
        cdef c.Dbn2D *p = self.d2ptr()
        if self._deallocate:
            del p


    def __init__(self):
        cutil.set_owned_ptr(self, new c.Dbn2D())

    def __repr__(self):
        mean = self.mean if self.sumW > 0 else None
        sd = self.stdDev if self.sumW > 0 else None
        return '<Dbn2D(mean=%s, stdDev=%s)>' % (mean, sd)


    def copy(self):
        return cutil.new_owned_cls(Dbn2D, new c.Dbn2D(deref(self.d2ptr())))

    def reset(self):
        """
        () -> None

        Reset the distribution counters to the unfilled state."""
        self.d2ptr().reset()


    def fill(self, x, y, weight=1.0, fraction=1.0):
        """
        (x, y, weight=1.0) -> None

        Fills the distribution with the given weight at given (x, y).

        """
        self.d2ptr().fill(x, y, weight, fraction)


    def scaleW(self, w):
        """
        (float) -> None

        Scale the weights by the given factor.
        """
        self.d2ptr().scaleW(w)

    def scaleX(self, x):
        """
        (float) -> None

        Scale the x dimension by the given factor.
        """
        self.d2ptr().scaleX(x)

    def scaleY(self, y):
        """
        (float) -> None

        Scale the y dimension by the given factor.
        """
        self.d2ptr().scaleY(y)

    def scaleXY(self, x, y):
        """
        (float, float) -> None

        Scale the x and y dimensions by the given factors.
        """
        self.d2ptr().scaleXY(x, y)


    # TODO: map direct properties from C++

    #@property
    def mean(self):
        """Weighted mean of x"""
        return util.XY(self.d2ptr().xMean(), self.d2ptr().yMean())

    #@property
    def variance(self):
        """Weighted variance of x"""
        return util.XY(self.d2ptr().xVariance(), self.d2ptr().yVariance())

    #@property
    def stdDev(self):
        """Weighted standard deviation of x"""
        return util.XY(self.d2ptr().xStdDev(), self.d2ptr().yStdDev())

    #@property
    def stdErr(self):
        """Weighted standard error on <x>"""
        return util.XY(self.d2ptr().xStdErr(), self.d2ptr().yStdErr())

    #@property
    def rms(self):
        """Weighted root mean squared (RMS) of x"""
        return util.XY(self.d2ptr().xRMS(), self.d2ptr().yRMS())


    #@property
    def numEntries(self):
        """The number of entries"""
        return self.d2ptr().numEntries()

    #@property
    def effNumEntries(self):
        """Effective number of entries (for weighted events)"""
        return self.d2ptr().effNumEntries()


    #@property
    def errW(self):
        """Error on sumW"""
        return self.d2ptr().errW()

    #@property
    def relErrW(self):
        """Relative error on sumW"""
        return self.d2ptr().relErrW()


    #@property
    def sumW(self):
        """sum(weights)"""
        return self.d2ptr().sumW()

    #@property
    def sumW2(self):
        """sum(weights * weights)"""
        return self.d2ptr().sumW2()

    #@property
    def sumWX(self):
        """sum(weights * xs)"""
        return self.d2ptr().sumWX()

    #@property
    def sumWY(self):
        """sum(weights * ys)"""
        return self.d2ptr().sumWY()

    #@property
    def sumWX2(self):
        """sum(weights * xs * xs)"""
        return self.d2ptr().sumWX2()

    #@property
    def sumWY2(self):
        """sum(weights * ys * ys)"""
        return self.d2ptr().sumWY2()

    #@property
    def sumWXY(self):
        """sum(weights xs * ys)"""
        return self.d2ptr().sumWXY()


    def __add__(Dbn2D self, Dbn2D other):
        return cutil.new_owned_cls(Dbn2D, new c.Dbn2D(deref(self.d2ptr()) + deref(other.d2ptr())))

    def __sub__(Dbn2D self, Dbn2D other):
        return cutil.new_owned_cls(Dbn2D, new c.Dbn2D(deref(self.d2ptr()) - deref(other.d2ptr())))
