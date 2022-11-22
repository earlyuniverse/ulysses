cimport util
cdef class Dbn3D(util.Base):
    """
    A 3D distribution 'counter', used and exposed by 2D profiles and their bins.

    TODO: also provide normal scalar access to quantities like xRMS
    """

    cdef c.Dbn3D* d3ptr(self) except NULL:
        return <c.Dbn3D*> self.ptr()

    # TODO: remove
    cdef c.Dbn3D* _Dbn3D(self) except NULL:
        return <c.Dbn3D*> self.ptr()

    def __dealloc__(self):
        cdef c.Dbn3D *p = self.d3ptr()
        if self._deallocate:
            del p


    def __init__(self):
        cutil.set_owned_ptr(self, new c.Dbn3D())

    def __repr__(self):
        mean = self.mean if self.sumW > 0 else None
        sd = self.stdDev if self.sumW > 0 else None
        return 'Dbn3D(mean=%s, stddev=%s)' % (mean, sd)


    def copy(self):
        return cutil.new_owned_cls(Dbn3D, new c.Dbn3D(deref(self.d3ptr())))


    def fill(self, x, y, z, weight=1.0, fraction=1.0):
        """
        (x, y, z, weight=1.0) -> None

        Fills the distribution with the given weight at given (x, y).

        """
        self.d3ptr().fill(x, y, z, weight, fraction)

    def reset(self):
        """
        () -> None

        Reset the distribution counters to the unfilled state."""
        self.d3ptr().reset()


    def scaleW(self, w):
        """
        (float) -> None

        Scale the weights by the given factor.
        """
        self.d3ptr().scaleW(w)

    def scaleX(self, x):
        """
        (float) -> None

        Scale the x dimension by the given factor.
        """
        self.d3ptr().scaleX(x)

    def scaleY(self, y):
        """
        (float) -> None

        Scale the y dimension by the given factor.
        """
        self.d3ptr().scaleY(y)

    def scaleZ(self, z):
        """
        (float) -> None

        Scale the z dimension by the given factor.
        """
        self.d3ptr().scaleZ(z)

    def scaleXYZ(self, x, y, z):
        """
        (float, float, float) -> None

        Scale the x, y and z dimensions by the given factors.
        """
        self.d3ptr().scaleXYZ(x, y, z)


    # TODO: map direct properties from C++

    #@property
    def mean(self):
        """Weighted mean of x"""
        return util.XYZ(self.d3ptr().xMean(),
                        self.d3ptr().yMean(),
                        self.d3ptr().zMean())

    #@property
    def variance(self):
        """Weighted variance of x"""
        return util.XYZ(self.d3ptr().xVariance(),
                        self.d3ptr().yVariance(),
                        self.d3ptr().zVariance())

    #@property
    def stdDev(self):
        """Weighted standard deviation of x"""
        return util.XYZ(self.d3ptr().xStdDev(),
                        self.d3ptr().yStdDev(),
                        self.d3ptr().zStdDev())

    #@property
    def stdErr(self):
        """Weighted standard error on <x>"""
        return util.XYZ(self.d3ptr().xStdErr(),
                        self.d3ptr().yStdErr(),
                        self.d3ptr().zStdErr())

    #@property
    def rms(self):
        """Weighted root mean squared (RMS) of x"""
        return util.XYZ(self.d3ptr().xRMS(),
                        self.d3ptr().yRMS(),
                        self.d3ptr().zRMS())


    #@property
    def numEntries(self):
        """The number of entries"""
        return self.d3ptr().numEntries()

    #@property
    def effNumEntries(self):
        """Effective number of entries (for weighted events)"""
        return self.d3ptr().effNumEntries()


    #@property
    def errW(self):
        """Error on sumW"""
        return self.d3ptr().errW()

    #@property
    def relErrW(self):
        """Relative error on sumW"""
        return self.d3ptr().relErrW()


    #@property
    def sumW(self):
        """sum(weights)"""
        return self.d3ptr().sumW()

    #@property
    def sumW2(self):
        """sum(weights * weights)"""
        return self.d3ptr().sumW2()

    #@property
    def sumWX(self):
        """sum(weights * xs)"""
        return self.d3ptr().sumWX()

    #@property
    def sumWY(self):
        """sum(weights * ys)"""
        return self.d3ptr().sumWY()

    #@property
    def sumWZ(self):
        """sum(weights * zs)"""
        return self.d3ptr().sumWZ()

    #@property
    def sumWX2(self):
        """sum(weights * xs * xs)"""
        return self.d3ptr().sumWX2()

    #@property
    def sumWY2(self):
        """sum(weights * ys * ys)"""
        return self.d3ptr().sumWY2()

    #@property
    def sumWZ2(self):
        """sum(weights * zs * zs)"""
        return self.d3ptr().sumWZ2()

    #@property
    def sumWXY(self):
        """sum(weights * xs * ys)"""
        return self.d3ptr().sumWXY()

    #@property
    def sumWXZ(self):
        """sum(weights * xs * zs)"""
        return self.d3ptr().sumWXZ()

    #@property
    def sumWYZ(self):
        """sum(weights * ys * zs)"""
        return self.d3ptr().sumWYZ()


    def __add__(Dbn3D self, Dbn3D other):
        return cutil.new_owned_cls(Dbn3D, new c.Dbn3D(deref(self.d3ptr()) + deref(other.d3ptr())))

    def __sub__(Dbn3D self, Dbn3D other):
        return cutil.new_owned_cls(Dbn3D, new c.Dbn3D(deref(self.d3ptr()) - deref(other.d3ptr())))
