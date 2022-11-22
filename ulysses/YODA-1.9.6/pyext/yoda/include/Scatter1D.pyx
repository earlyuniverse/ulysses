cimport util
cdef class Scatter1D(AnalysisObject):
    """
    1D scatter plot, i.e. a collection of Point1D objects with positions and errors.

    Constructor calling idioms:

    Scatter1D(path="", title="")
      Create a new empty scatter, with optional path and title.

    Scatter1D(points, path="", title=""):
      Create a new empty scatter from an iterable of points, with optional path
      and title.

    TODO: more documentation!
    """

    cdef inline c.Scatter1D* s1ptr(self) except NULL:
        return <c.Scatter1D*> self.ptr()


    def __init__(self, *args, **kwargs):
        util.try_loop([self.__init_2, self.__init_3], *args, **kwargs)

    def __init_2(self, path="", title=""):
        path  = path.encode('utf-8')
        title = title.encode('utf-8')
        cutil.set_owned_ptr(self, new c.Scatter1D(<string>path, <string>title))

    def __init_3(self, points, path="", title=""):
        self.__init_2(path, title)
        self.addPoints(points)

    def clone(self):
        """() -> Scatter1D.
        Clone this Scatter1D."""
        return cutil.new_owned_cls(Scatter1D, self.s1ptr().newclone())

    def __repr__(self):
        return "<%s '%s' %d points>" % (self.__class__.__name__, self.path(), len(self.points()))


    def reset(self):
        "Reset the scatter, removing all points"
        self.s1ptr().reset()


    def numPoints(self):
        """() -> int
        Number of points in this scatter."""
        return self.s1ptr().numPoints()

    def __len__(self):
        return self.numPoints()


    def points(self):
        """Access the ordered list of points."""
        return [self.point(i) for i in range(self.numPoints())]

    def point(self, size_t i):
        """Access the i'th point."""
        return cutil.new_borrowed_cls(Point1D, &self.s1ptr().point(i), self)

    def __getitem__(self, py_ix):
        cdef size_t i = cutil.pythonic_index(py_ix, self.numPoints())
        return cutil.new_borrowed_cls(Point1D, &self.s1ptr().point(i), self)


    def addPoint(self, *args, **kwargs):
        """Add a new point.

        Provide either a single yoda.Point1D object, or the
        two args: x, xerrs=0.
        """
        try:
            self.__addPoint_point(*args, **kwargs)
        except TypeError:
            self.__addPoint_explicit(*args, **kwargs)

    def __addPoint_explicit(self, x, xerrs=0):
        self.__addPoint_point(Point1D(x, xerrs))

    def __addPoint_point(self, Point1D p):
        self.s1ptr().addPoint(p.p1ptr()[0])

    def addPoints(self, iterable):
        """Add several new points."""
        for row in iterable:
          try:
            self.addPoint(*row)
          except TypeError:
            self.addPoint(row)

    def rmPoint(self, idx):
        self.s1ptr().rmPoint(idx)

    def rmPoints(self, idxs):
        self.s1ptr().rmPoints(idxs)

    def combineWith(self, others):
        """Try to add points from other Scatter1Ds into this one."""
        cdef Scatter1D other
        try:
            # Can we type it as a Scatter1D?
            other = others
        except TypeError:
            # Could be an iterable...
            for other in others:
                self.s1ptr().combineWith(deref(other.s1ptr()))
        else:
            self.s1ptr().combineWith(deref(other.s1ptr()))


    def mkScatter(self):
        """None -> Scatter1D.
        Make a new Scatter1D. Exists to allow mkScatter calls on any AnalysisObject,
        even if it already is a scatter."""
        cdef c.Scatter1D s2 = c.mkScatter_Scatter1D(deref(self.s1ptr()))
        return cutil.new_owned_cls(Scatter1D, s2.newclone())


    def scaleX(self, a):
        """(float) -> None
        Scale the x values and errors of the points in this scatter by factor a."""
        self.s1ptr().scaleX(a)


    def scale(self, i, scale):
        """(int, float) -> None
        Scale values on axis i"""
        self.s1ptr().scale(i, scale)


    def transformX(self, f):
        """(fn) -> None
        Transform the x values and errors of the points in this scatter by function f."""
        import ctypes
        try:
            callback = ctypes.CFUNCTYPE(ctypes.c_double, ctypes.c_double)(f)
        except:
            raise RuntimeError("Callback is not of type (double) -> double")
        fptr = (<c.dbl_dbl_fptr*><size_t>ctypes.addressof(callback))[0]
        c.Scatter1D_transformX(deref(self.s1ptr()), fptr)
    
    def parseVariations(self):
        """None -> None
        Parse the YAML which contains the variations stored in the points of the Scatter.
        Only needs to be done once!"""
        return self.s1ptr().parseVariations()
    
    def updateTotalUncertainty(self):
        """None -> None
        Sum the error in the error map in quadrature and update the total ("") uncertainty.
        """
        return self.s1ptr().updateTotalUncertainty()
    
    def writeVariationsToAnnotations(self):
        """None -> None
        Parse the variations and add them to the annotations.
        """
        return self.s1ptr().writeVariationsToAnnotations()

    def variations(self):
        """None -> vector[string]
        Get the list of variations stored in the points of the Scatter"""
        cdef vector[string] vars = self.s1ptr().variations()
        return [var.decode('utf-8') for var in vars]

    def rmVariations(self):
        """None -> None
        Remove the variations stored in the points of the Scatter"""
        return self.s1ptr().rmVariations()
    
    def hasValidErrorBreakdown(self):
        """
        Check if the AO's error breakdown is not empty and has no bins with 0 uncertainty
        """
        counter = -1
        for p in self.points():
            counter += 1
            binErrs = p.errMap()
            if len(binErrs) < 2:
                return False
            binTotal = [0.,0.]
            for sys, err in binErrs.iteritems():
                binTotal[0] = (binTotal[0]**2 + err[0]**2)**0.5
                binTotal[1] = (binTotal[1]**2 + err[1]**2)**0.5
            if binTotal[0] == 0 and binTotal[1] == 0:
                return False
        return True

    # # TODO: remove?
    # def __add__(Scatter1D self, Scatter1D other):
    #     return cutil.new_owned_cls(Scatter1D, c.Scatter1D_add_Scatter1D(self.s1ptr(), other.s1ptr()))

    # # TODO: remove?
    # def __sub__(Scatter1D self, Scatter1D other):
    #     return cutil.new_owned_cls(Scatter1D, c.Scatter1D_sub_Scatter1D(self.s1ptr(), other.s1ptr()))


    def _mknp(self, xs):
        try:
            import numpy
            return numpy.array(xs)
        except ImportError:
            return xs

    def xVals(self):
        return self._mknp([p.x() for p in self.points()])

    def xMins(self):
        """All x low values."""
        return self._mknp([p.xMin() for p in self.points()])

    def xMaxs(self):
        """All x high values."""
        return self._mknp([p.xMax() for p in self.points()])

    # TODO: xErrs

    def xMin(self):
        """Lowest x value."""
        return min(self.xMins())

    def xMax(self):
        """Highest x value."""
        return max(self.xMaxs())


## Convenience alias
S1D = Scatter1D
