cimport util
cdef class Dbn0D(util.Base):
    """
    A zero-dimensional 'counter', used and exposed by Counter.
    """

    cdef c.Dbn0D* d0ptr(self) except NULL:
        return <c.Dbn0D *> self.ptr()
    # TODO: remove!
    cdef c.Dbn0D *_Dbn0D(self) except NULL:
        return <c.Dbn0D *> self.ptr()

    def __dealloc__(self):
        cdef c.Dbn0D *p = self.d0ptr()
        if self._deallocate:
            del p


    def __init__(self):
        cutil.set_owned_ptr(self, new c.Dbn0D())

    def __repr__(self):
        return '<Dbn0D(val=%g, err=%g)>' % (self.val, self.err)


    def copy(self):
        return cutil.set_owned_ptr(self, new c.Dbn0D(deref(self.d0ptr())))

    def reset(self):
        """
        () -> None

        Reset the distribution counters to the unfilled state.
        """
        self.d0ptr().reset()


    def fill(self, weight=1.0, fraction=1.0):
        """
        (float weight=1.0) -> None

        Fills the distribution with the given weight at given x.
        """
        self.d0ptr().fill(weight, fraction)

    def scaleW(self, w):
        """
        (float) -> None

        Scale the weights by the given factor.
        """
        self.d0ptr().scaleW(w)


    #@property
    def numEntries(self):
        """The number of entries"""
        return self.d0ptr().numEntries()

    #@property
    def effNumEntries(self):
        """Effective number of entries (for weighted events)"""
        return self.d0ptr().effNumEntries()


    #@property
    def errW(self):
        """Error on sumW"""
        return self.d0ptr().errW()

    #@property
    def relErrW(self):
        """Relative error on sumW"""
        return self.d0ptr().relErrW()


    #@property
    def sumW(self):
        """sum(weights)"""
        return self.d0ptr().sumW()

    #@property
    def sumW2(self):
        """sum(weights * weights)"""
        return self.d0ptr().sumW2()


    def __add__(Dbn0D self, Dbn0D other):
        return cutil.new_owned_cls(Dbn0D, new c.Dbn0D(deref(self.d0ptr()) + deref(other.d0ptr())))

    def __sub__(Dbn0D self, Dbn0D other):
        return cutil.new_owned_cls(Dbn0D, new c.Dbn0D(deref(self.d0ptr()) - deref(other.d0ptr())))
