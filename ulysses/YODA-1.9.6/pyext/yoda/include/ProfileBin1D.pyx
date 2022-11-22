cdef class ProfileBin1D(Bin1D_Dbn2D):
    """
    A 1D profile bin, as stored inside Profile1D.

    Only one constructor:

    * ProfileBin1D(xlow, xhigh)
    """

    cdef inline c.ProfileBin1D* pb1ptr(self) except NULL:
        return <c.ProfileBin1D*> self.ptr()
    # TODO: remove
    cdef inline c.ProfileBin1D* _ProfileBin1D(self) except NULL:
        return <c.ProfileBin1D*> self.ptr()


    def __init__(self, double a, double b):
        cutil.set_owned_ptr(self, new c.ProfileBin1D(a, b))


    # def fill(self, x, y, weight=1.0, fraction=1.0):
    #     """
    #     (x, y, weight=1.0) -> None. Fill this bin with given values and weight.
    #     """
    #     self.pb1ptr().fill(x, y, weight, fraction)

    # def fillBin(self, y, weight=1.0, fraction=1.0):
    #     """
    #     (y, weight=1.0) -> None. Fill this bin with given y-value and weight.
    #     """
    #     self.pb1ptr().fillBin(y, weight, fraction)


    # def scaleY(self, ay):
    #     """
    #     float -> None

    #     Scale y values by ay.
    #     """
    #     self.pb1ptr().scaleY(ay)


    def mean(self):
        """The mean of the y-values that have filled the bin."""
        return self.pb1ptr().mean()

    def variance(self):
        """The variance of the y-values that have filled the bin."""
        return self.pb1ptr().variance()

    def stdDev(self):
        """The standard deviation of the y-values that have filled the bin."""
        return self.pb1ptr().stdDev()

    def stdErr(self):
        """The standard error of the y-values that have filled the bin."""
        return self.pb1ptr().stdErr()

    def rms(self):
        """The RMS of the y-values that have filled the bin."""
        return self.pb1ptr().rms()


    def sumWY(self):
        """sum(weights * ys)"""
        return self.pb1ptr().sumWY()

    def sumWY2(self):
        """sum(weights * ys * ys)"""
        return self.pb1ptr().sumWY2()


    def __iadd__(ProfileBin1D self, ProfileBin1D other):
        c.ProfileBin1D_iadd_ProfileBin1D(self.pb1ptr(), other.pb1ptr())
        return self

    def __isub__(ProfileBin1D self, ProfileBin1D other):
        c.ProfileBin1D_isub_ProfileBin1D(self.pb1ptr(), other.pb1ptr())
        return self

    def __add__(ProfileBin1D a, ProfileBin1D b):
        return cutil.new_owned_cls(ProfileBin1D,
                                  new c.ProfileBin1D(deref(a.pb1ptr()) + deref(b.pb1ptr())))

    def __sub__(ProfileBin1D a, ProfileBin1D b):
        return cutil.new_owned_cls(ProfileBin1D,
                                   new c.ProfileBin1D(deref(a.pb1ptr()) - deref(b.pb1ptr())))

    def __repr__(self):
        return 'ProfileBin1D(%g, %g; sumw=%g)' % (self.xEdges()[0], self.xEdges()[1], self.mean())
