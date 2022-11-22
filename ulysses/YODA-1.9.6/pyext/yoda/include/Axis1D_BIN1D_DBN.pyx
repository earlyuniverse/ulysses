cimport util
# TODO (when there is absolutely nothing else to do) docstrings (but never will
# it be a user facing class... it's merely there for tests)
cdef class Axis1D_${BIN1D}_${DBN}(util.Base):

    cdef inline c.Axis1D[c.${BIN1D}, c.${DBN}]* a1dptr(self) except NULL:
        return <c.Axis1D[c.${BIN1D}, c.${DBN}]*> self.ptr()
    # TODO: remove
    cdef inline c.Axis1D[c.${BIN1D}, c.${DBN}]* _Axis1D(self) except NULL:
        return <c.Axis1D[c.${BIN1D}, c.${DBN}]*> self.ptr()

    def __dealloc__(self):
        cdef c.Axis1D[c.${BIN1D}, c.${DBN}]* p = self.a1ptr()
        if self._deallocate:
            del p


    def __init__(self):
        cutil.set_owned_ptr(self, new c.Axis1D[c.${BIN1D}, c.${DBN}]())

    def __repr__(self):
        return "<Axis1D with %d bins>" % self.numBins


    #@property
    def numBins(self):
        return self.a1ptr().bins().size()

    def __len__(self):
        return self.numBins

    # TODO: remove
    # def __getitem__(self, py_ix):
    #     cdef size_t i = cutil.pythonic_index(py_ix, self.a1ptr().bins().size())
    #     return cutil.new_borrowed_cls(${BIN1D}, & self.a1ptr().bins().at(i), self)


    def addBin(self, a, b):
        self.a1ptr().addBin(a, b)

    #@property
    def totalDbn(self):
        return cutil.new_borrowed_cls(${DBN}, &self.a1ptr().totalDbn(), self)

    #@property
    def underflow(self):
        return cutil.new_borrowed_cls(${DBN}, &self.a1ptr().underflow(), self)

    #@property
    def overflow(self):
        return cutil.new_borrowed_cls(${DBN}, &self.a1ptr().overflow(), self)

    def reset(self):
        self.a1ptr().reset()

    def eraseBin(self, i):
        self.a1ptr().eraseBin(i)

    def getBinIndex(self, x):
        return self.a1ptr().getBinIndex(x)

    def mergeBins(self, a, b):
        self.a1ptr().mergeBins(a, b)

    #def binAt(self, x):
    #    return self[self.a1ptr().getBinIndex(x)]
