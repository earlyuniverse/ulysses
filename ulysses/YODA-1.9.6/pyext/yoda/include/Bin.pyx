cimport util
cdef class Bin(util.Base):

    cdef inline c.Bin* bptr(self) except NULL:
        return <c.Bin*> self.ptr()

    def __dealloc__(self):
        cdef c.Bin* p = self.bptr()
        if self._deallocate:
            del p


    #@property
    def dim(self):
        """None -> int
        Dimension of the fill space (should match containing Histo/Profile)"""
        return self.bptr().dim()


    #@property
    def numEntries(self):
        """
        The number of entries that have filled the bin.
        """
        return self.bptr().numEntries()

    #@property
    def effNumEntries(self):
        """
        The effective number of entries in the bin.

        s.effNumEntries <==> (s.sumW ** 2) / s.sumW2
        """
        return self.bptr().effNumEntries()


    #@property
    def sumW(self):
        """
        The sum of weights: sum(weights).
        """
        return self.bptr().sumW()

    #@property
    def sumW2(self):
        """
        The sum of weights-squared: sum(weights * weights)
        """
        return self.bptr().sumW2()
