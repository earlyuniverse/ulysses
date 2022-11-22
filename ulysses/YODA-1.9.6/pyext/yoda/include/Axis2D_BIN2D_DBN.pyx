cimport util
# TODO: docstrings
cdef class Axis2D_${BIN2D}_${DBN}(util.Base):

    cdef inline c.Axis2D[c.${BIN2D}, c.${DBN}]* a2ptr(self) except NULL:
        return <c.Axis2D[c.${BIN2D}, c.${DBN}]*> self.ptr()
    # TODO: remove
    cdef inline c.Axis2D[c.${BIN2D}, c.${DBN}]* _Axis2D(self) except NULL:
        return <c.Axis2D[c.${BIN2D}, c.${DBN}]*> self.ptr()

    def __dealloc__(self):
        cdef c.Axis2D[c.${BIN2D}, c.${DBN}]* p = self.a2ptr()
        if self._deallocate:
            del p


    def __init__(self, nx, xl, xu, ny, yl, yu):
        cutil.set_owned_ptr(self, new c.Axis2D[c.${BIN2D}, c.${DBN}](
            nx, pair[double, double](xl, xu),
            ny, pair[double, double](yl, yu)))


    #@property
    def numBins(self):
        return self._Axis1D().bins().size()

    def __len__(self):
        return self.numBins

    # TODO: remove
    # def __getitem__(self, py_ix):
    #     cdef size_t i = cutil.pythonic_index(py_ix, self.a2ptr().bins().size())
    #     return cutil.new_borrowed_cls(${BIN2D}, & self.a2ptr().bins().at(i), self)

    def __repr__(self):
        # TODO: improve
        return "<Axis2D with %d bins>" % self.numBins


    #@property
    def totalDbn(self):
        return cutil.new_owned_cls(
            ${DBN}, new c.${DBN}(self.a2ptr().totalDbn()))

    def addBin(self, a, b, c, d):
        self.a2ptr().addBin(a, b, c, d)

    #@property
    def outflow(self, ix, iy):
        return cutil.new_owned_cls(${DBN}, new c.${DBN}(self.a2ptr().outflow(ix, iy)))

    #@property
    def edges(self):
        return util.XY(
            util.EdgePair(self.a2ptr().xMin(), self.a2ptr().xMax()),
            util.EdgePair(self.a2ptr().yMin(), self.a2ptr().yMax())
        )

    def reset(self):
        self.a2ptr().reset()

    def binAt(self, x, y):
        cdef int ix = self.a2ptr().getBinIndex(x, y)
        if ix < 0:
            raise YodaExc_RangeError('No bin found!')
        return self[ix]
