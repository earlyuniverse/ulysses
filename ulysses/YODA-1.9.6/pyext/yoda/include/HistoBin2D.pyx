# TODO: tidy once we have a working Histo2D
cdef class HistoBin2D(Bin2D_Dbn2D):

    cdef inline c.HistoBin2D* hb2ptr(self) except NULL:
        return <c.HistoBin2D*> self.ptr()
    # TODO: remove
    cdef inline c.HistoBin2D* _HistoBin2D(self) except NULL:
        return <c.HistoBin2D*> self.ptr()

    def __init__(self, xlow, xhigh, ylow, yhigh):
        cutil.set_owned_ptr(self, new c.HistoBin2D(xlow, xhigh, ylow, yhigh))


    # def fill(self, x, y, weight=1.0, fraction=1.0):
    #     self.hb2ptr().fill(x, y, weight, fraction)

    #@property
    def volume(self):
        return self.hb2ptr().volume()

    #@property
    def height(self):
        return self.hb2ptr().height()

    #@property
    def volumeErr(self):
        return self.hb2ptr().volumeErr()

    #@property
    def heightErr(self):
        return self.hb2ptr().heightErr()

    #@property
    def relErr(self):
        return self.hb2ptr().relErr()


    def __add__(HistoBin2D a, HistoBin2D b):
        return cutil.new_owned_cls(HistoBin2D, new c.HistoBin2D(deref(a.hb2ptr()) + deref(b.hb2ptr())))

    def __sub__(HistoBin2D a, HistoBin2D b):
        return cutil.new_owned_cls(HistoBin2D, new c.HistoBin2D(deref(a.hb2ptr()) - deref(b.hb2ptr())))

    def __repr__(self):
        return 'HistoBin2D(%g, %g; %g, %g; sumw=%g)' % (self.xEdges()[0], self.xEdges()[1],
                                                        self.yEdges()[0], self.yEdges()[1], self.sumW())
