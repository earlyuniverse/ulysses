#################################################
#############       WARNING      ################
#################################################
# This file has been automatically generated.
# Any changes you make here will get overridden.
# Instead, make your changes in Bin1D_DBN.pyx
#################################################
cdef class Bin1D_Dbn1D(Bin):
    """
    1D Bin based on an underlying Dbn1D.

    ProfileBin1D(xlow, xhigh)

    """

    cdef inline c.Bin1D_Dbn1D* b1ptr(self) except NULL:
        return <c.Bin1D_Dbn1D*> self.ptr()
    # TODO: remove
    cdef inline c.Bin1D_Dbn1D* _Bin1D(self) except NULL:
        return <c.Bin1D_Dbn1D*> self.ptr()


    def __init__(self, xlow, xhigh):
        cutil.set_owned_ptr(self, new c.Bin1D_Dbn1D(pair[double, double](xlow, xhigh)))

    def __repr__(self):
        return '<%s x=[%g, %g)>' % (self.__class__.__name__, self.xMin(), self.xMax())


    # def scaleX(self, ax):
    #     """
    #     float -> None
    #     Scale this bin's x arguments by ax.
    #     """
    #     self.b1ptr().scaleX(ax)

    # def scaleW(self, a):
    #     """
    #     float -> None
    #     Scale this bin's weights by a.
    #     """
    #     self.b1ptr().scaleW(a)


    def xEdges(self):
        """The lower and upper edges."""
        return (self.xMin(), self.xMax())

    def xMin(self):
        """The lower bin edge."""
        return self.b1ptr().xMin()

    def xMax(self):
        """The upper bin edge."""
        return self.b1ptr().xMax()

    def xMid(self):
        """The midpoint of the bin."""
        return self.b1ptr().xMid()

    def xWidth(self):
        """The width of the bin."""
        return self.b1ptr().xWidth()

    def xMean(self):
        """The mean of the x-values that have filled the bin."""
        return self.b1ptr().xMean()

    def xFocus(self):
        """
        The focus of the bin. If the bin has been filled, then this
        is the mean fill on this bin. If the bin has not been filled,
        then the focus is the midpoint of the bin.
        """
        return self.b1ptr().xFocus()

    def xVariance(self):
        """
        The variance of the x-values that have filled the bin.
        """
        return self.b1ptr().xVariance()

    def xStdDev(self):
        """
        The standard deviation of the x-values that have filled the bin.
        """
        return self.b1ptr().xStdDev()

    def xStdErr(self):
        """
        The standard error of the x-values that have filled the bin.
        """
        return self.b1ptr().xStdErr()

    def xRMS(self):
        """
        The root-mean-square of the x-values that have filled the bin.
        """
        return self.b1ptr().xRMS()



    def sumWX(self):
        """
        The sum of weights-times-x: sum(weights * x)
        """
        return self.b1ptr().sumWX()

    def sumWX2(self):
        """
        The sum of weights-times-x-squared: sum(weights * x * x)
        """
        return self.b1ptr().sumWX2()

    def merge(Bin1D_Dbn1D self, Bin1D_Dbn1D other):
        """
        merge(other) -> Bin1D_Dbn1D. Merge this bin with another of the
        same type. Only directly adjacent bins, i.e. those sharing a
        common edge, can be merged.
        """
        self.b1ptr().merge(deref(other.b1ptr()))
        return self



    def __add__(Bin1D_Dbn1D self, Bin1D_Dbn1D other):
        return cutil.new_owned_cls(Bin1D_Dbn1D,
                                   new c.Bin1D_Dbn1D(deref(self.b1ptr()) + deref(other.b1ptr())))

    def __sub__(Bin1D_Dbn1D self, Bin1D_Dbn1D other):
        return cutil.new_owned_cls(Bin1D_Dbn1D,
                                   new c.Bin1D_Dbn1D(deref(self.b1ptr()) - deref(other.b1ptr())))
