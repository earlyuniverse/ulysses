cimport util
cdef class Dbn1D(util.Base):
    """
    A 1D distribution 'counter', used and exposed by 1D histograms and their bins.
    """

    cdef c.Dbn1D* d1ptr(self) except NULL:
        return <c.Dbn1D *> self.ptr()
    # TODO: remove!
    cdef c.Dbn1D *_Dbn1D(self) except NULL:
        return <c.Dbn1D *> self.ptr()

    def __dealloc__(self):
        cdef c.Dbn1D *p = self.d1ptr()
        if self._deallocate:
            del p


    def __init__(self):
        cutil.set_owned_ptr(self, new c.Dbn1D())

    def __repr__(self):
        mean = self.mean if self.sumW > 0 else None
        sd = self.stdDev if self.sumW > 0 else None
        return '<Dbn1D(mean=%s, stddev=%s)>' % (mean, sd)


    def copy(self):
        return cutil.set_owned_ptr(self, new c.Dbn1D(deref(self.d1ptr())))

    def reset(self):
        """
        () -> None

        Reset the distribution counters to the unfilled state.
        """
        self.d1ptr().reset()


    def fill(self, x, weight=1.0, fraction=1.0):
        """
        (float x, float weight=1.0) -> None

        Fills the distribution with the given weight at given x.
        """
        self.d1ptr().fill(x, weight, fraction)

    def scaleW(self, w):
        """
        (float) -> None

        Scale the weights by the given factor.
        """
        self.d1ptr().scaleW(w)

    def scaleX(self, x):
        """
        (float) -> None

        Scale the x dimension by the given factor.
        """
        self.d1ptr().scaleX(x)


    #@property
    def xMean(self):
        """Weighted mean of x"""
        return self.d1ptr().xMean()

    #@property
    def xVariance(self):
        """Weighted variance of x"""
        return self.d1ptr().xVariance()

    #@property
    def xStdDev(self):
        """Weighted standard deviation of x"""
        return self.d1ptr().xStdDev()

    #@property
    def xStdErr(self):
        """Weighted standard error on <x>"""
        return self.d1ptr().xStdErr()

    #@property
    def xRMS(self):
        """Weighted root mean squared (RMS) of x"""
        return self.d1ptr().xRMS()

    #@property
    def numEntries(self):
        """The number of entries"""
        return self.d1ptr().numEntries()

    #@property
    def effNumEntries(self):
        """Effective number of entries (for weighted events)"""
        return self.d1ptr().effNumEntries()


    #@property
    def errW(self):
        """Error on sumW"""
        return self.d1ptr().errW()

    #@property
    def relErrW(self):
        """Relative error on sumW"""
        return self.d1ptr().relErrW()


    #@property
    def sumW(self):
        """sum(weights)"""
        return self.d1ptr().sumW()

    #@property
    def sumW2(self):
        """sum(weights * weights)"""
        return self.d1ptr().sumW2()

    #@property
    def sumWX(self):
        """sum(weights * xs)"""
        return self.d1ptr().sumWX()

    #@property
    def sumWX2(self):
        """sum(weights * xs * xs)"""
        return self.d1ptr().sumWX2()


    def __add__(Dbn1D self, Dbn1D other):
        return cutil.new_owned_cls(Dbn1D, new c.Dbn1D(deref(self.d1ptr()) + deref(other.d1ptr())))

    def __sub__(Dbn1D self, Dbn1D other):
        return cutil.new_owned_cls(Dbn1D, new c.Dbn1D(deref(self.d1ptr()) - deref(other.d1ptr())))
