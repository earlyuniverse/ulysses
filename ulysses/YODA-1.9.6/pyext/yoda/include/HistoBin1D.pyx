cdef class HistoBin1D(Bin1D_Dbn1D):

    cdef inline c.HistoBin1D* hb1ptr(self) except NULL:
        return <c.HistoBin1D*> self.ptr()
    # TODO: remove
    cdef inline c.HistoBin1D* _HistoBin1D(self) except NULL:
        return <c.HistoBin1D*> self.ptr()


    def __init__(self, double a, double b):
        cutil.set_owned_ptr(self, new c.HistoBin1D(a, b))


    # def fill(self, value, double weight=1.0, fraction=1.0):
    #     """
    #     (value=None, weight=1.0)

    #     Fill this bin with the given value and given weight.
    #     """
    #     self.hb1ptr().fill(value, weight, fraction)

    # def fillBin(self, weight=1.0, fraction=1.0):
    #     """
    #     (weight=1.0) -> None. Fill this bin with given weight.
    #     """
    #     self.hb1ptr().fillBin(weight, fraction)


    #@property
    def area(self):
        """
        b.area <==> b.sumW

        The area of the bin is the sum of weights of the bin; it is
        independent of width.
        """
        return self.hb1ptr().area()

    #@property
    def height(self):
        """
        b.height <==> b.area / b.width

        The height of the bin is defined as the area divided by the
        width.
        """
        return self.hb1ptr().height()

    #@property
    def areaErr(self):
        """
        Error computed using binomial statistics on squared sum of bin
        weights, i.e. s.areaErr = sqrt(s.sumW2)
        """
        return self.hb1ptr().areaErr()

    #@property
    def heightErr(self):
        """
        Height error - scales the s.areaError by the reciprocal of the
        bin width.
        """
        return self.hb1ptr().heightErr()

    #@property
    def relErr(self):
        """
        Relative error - same for either area or height interpretations.
        """
        return self.hb1ptr().relErr()



    def __iadd__(HistoBin1D self, HistoBin1D other):
        c.HistoBin1D_iadd_HistoBin1D(self.hb1ptr(), other.hb1ptr())
        return self

    def __isub__(HistoBin1D self, HistoBin1D other):
        c.HistoBin1D_isub_HistoBin1D(self.hb1ptr(), other.hb1ptr())
        return self

    def __add__(HistoBin1D a, HistoBin1D b):
        return cutil.new_owned_cls(HistoBin1D,
                                   new c.HistoBin1D(deref(a.hb1ptr()) + deref(b.hb1ptr())))

    def __sub__(HistoBin1D a, HistoBin1D b):
        return cutil.new_owned_cls(HistoBin1D,
                                   new c.HistoBin1D(deref(a.hb1ptr()) - deref(b.hb1ptr())))

    def __repr__(self):
        return 'HistoBin1D(%g, %g; sumw=%g)' % (self.xEdges()[0], self.xEdges()[1], self.sumW())
