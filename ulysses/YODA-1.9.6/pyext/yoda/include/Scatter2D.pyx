cimport util
cdef class Scatter2D(AnalysisObject):
    """
    2D scatter plot, i.e. a collection of Point2D objects with positions and errors.

    Constructor calling idioms:

    Scatter2D(path="", title="")
      Create a new empty scatter, with optional path and title.

    Scatter2D(points, path="", title=""):
      Create a new empty scatter from an iterable of points, with optional path
      and title.

    TODO: more documentation!
    """

    cdef inline c.Scatter2D* s2ptr(self) except NULL:
        return <c.Scatter2D*> self.ptr()

    def __init__(self, *args, **kwargs):
        util.try_loop([self.__init_2, self.__init_3], *args, **kwargs)

    def __init_2(self, path="", title=""):
        path  = path.encode('utf-8')
        title = title.encode('utf-8')
        cutil.set_owned_ptr(self, new c.Scatter2D(<string>path, <string>title))

    def __init_3(self, points, path="", title=""):
        self.__init_2(path, title)
        self.addPoints(points)

    def clone(self):
        """() -> Scatter2D.
        Clone this Scatter2D."""
        return cutil.new_owned_cls(Scatter2D, self.s2ptr().newclone())

    def __repr__(self):
        return "<%s '%s' %d points>" % (self.__class__.__name__, self.path(), len(self.points()))


    def reset(self):
        "Reset the scatter, removing all points"
        self.s2ptr().reset()


    def numPoints(self):
        """() -> int
        Number of points in this scatter."""
        return self.s2ptr().numPoints()

    def __len__(self):
        return self.numPoints()


    def points(self):
        """Access the ordered list of points."""
        return [self.point(i) for i in range(self.numPoints())]

    def point(self, size_t i):
        """Access the i'th point."""
        return cutil.new_borrowed_cls(Point2D, &self.s2ptr().point(i), self)

    def __getitem__(self, py_ix):
        cdef size_t i = cutil.pythonic_index(py_ix, self.numPoints())
        return cutil.new_borrowed_cls(Point2D, &self.s2ptr().point(i), self)


    def addPoint(self, *args, **kwargs):
        """Add a new point.

        Provide either a single yoda.Point2D object, or the
        four args: x, y, xerrs=0, yerrs=0.
        """
        try:
            self.__addPoint_point(*args, **kwargs)
        except TypeError:
            self.__addPoint_explicit(*args, **kwargs)

    def __addPoint_explicit(self, x, y, xerrs=0, yerrs=0):
        self.__addPoint_point(Point2D(x, y, xerrs, yerrs))

    def __addPoint_point(self, Point2D p):
        self.s2ptr().addPoint(p.p2ptr()[0])

    def addPoints(self, iterable):
        """Add several new points."""
        for row in iterable:
          try:
            self.addPoint(*row)
          except TypeError:
            self.addPoint(row)

    def rmPoint(self, idx):
        self.s2ptr().rmPoint(idx)

    def rmPoints(self, idxs):
        self.s2ptr().rmPoints(idxs)

    def combineWith(self, others):
        """Try to add points from other Scatter2Ds into this one."""
        cdef Scatter2D other
        try:
            # Can we type it as a Scatter2D?
            other = others
        except TypeError:
            # Could be an iterable...
            for other in others:
                self.s2ptr().combineWith(deref(other.s2ptr()))
        else:
            self.s2ptr().combineWith(deref(other.s2ptr()))


    def mkScatter(self):
        """None -> Scatter2D.
        Make a new Scatter2D. Exists to allow mkScatter calls on any AnalysisObject,
        even if it already is a scatter."""
        cdef c.Scatter2D s2 = c.mkScatter_Scatter2D(deref(self.s2ptr()))
        return cutil.new_owned_cls(Scatter2D, s2.newclone())


    def scaleX(self, a):
        """(float) -> None
        Scale the x values and errors of the points in this scatter by factor a."""
        self.s2ptr().scaleX(a)

    def scaleY(self, a):
        """(float) -> None
        Scale the y values and errors of the points in this scatter by factor a."""
        self.s2ptr().scaleY(a)

    def scaleXY(self, ax=1.0, ay=1.0):
        """(float=1, float=1) -> None
        Scale the values and errors of the points in this scatter by factors ax, ay."""
        self.s2ptr().scaleXY(ax, ay)

    def scale(self, i, scale):
        """(int, float) -> None
        Scale values on axis i"""
        self.s2ptr().scale(i, scale)


    def transformX(self, f):
        """(fn) -> None
        Transform the x values and errors of the points in this scatter by function f."""
        import ctypes
        try:
            callback = ctypes.CFUNCTYPE(ctypes.c_double, ctypes.c_double)(f)
        except:
            raise RuntimeError("Callback is not of type (double) -> double")
        fptr = (<c.dbl_dbl_fptr*><size_t>ctypes.addressof(callback))[0]
        c.Scatter2D_transformX(deref(self.s2ptr()), fptr)

    def transformY(self, f):
        """(fn) -> None
        Transform the y values and errors of the points in this scatter by function f."""
        import ctypes
        try:
            callback = ctypes.CFUNCTYPE(ctypes.c_double, ctypes.c_double)(f)
        except:
            raise RuntimeError("Callback is not of type (double) -> double")
        fptr = (<c.dbl_dbl_fptr*><size_t>ctypes.addressof(callback))[0]
        c.Scatter2D_transformY(deref(self.s2ptr()), fptr)

    def parseVariations(self):
        """None -> None
        Parse the YAML which contains the variations stored in the points of the Scatter.
        Only needs to be done once!"""
        return self.s2ptr().parseVariations()

    def updateTotalUncertainty(self):
        """None -> None
        Sum the error in the error map in quadrature and update the total ("") uncertainty.
        """
        return self.s2ptr().updateTotalUncertainty()

    def writeVariationsToAnnotations(self):
        """None -> None
        Parse the variations and add them to the annotations.
        """
        return self.s2ptr().writeVariationsToAnnotations()

    def variations(self):
        """None -> vector[string]
        Get the list of variations stored in the points of the Scatter"""
        vars = self.s2ptr().variations()
        return [var.decode('utf-8') for var in vars]

    def rmVariations(self):
        """None -> None
        Remove the variations stored in the points of the Scatter"""
        return self.s2ptr().rmVariations()
    
    def _mknp(self, xs):
        try:
            import numpy
            return numpy.array(xs)
        except ImportError:
            return xs


    def covarianceMatrix(self, ignoreOffDiagonalTerms=False):
        """bool -> vector[vector[float]]
        Construct the covariance matrix"""
        return self._mknp(self.s2ptr().covarianceMatrix(ignoreOffDiagonalTerms))

    # # TODO: remove?
    # def __add__(Scatter2D self, Scatter2D other):
    #     return cutil.new_owned_cls(Scatter2D, c.Scatter2D_add_Scatter2D(self.s2ptr(), other.s2ptr()))

    # # TODO: remove?
    # def __sub__(Scatter2D self, Scatter2D other):
    #     return cutil.new_owned_cls(Scatter2D, c.Scatter2D_sub_Scatter2D(self.s2ptr(), other.s2ptr()))

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


    def correlationMatrix(self):
        """
        `covMatrix` numpy matrix
         Convert a covariance matrix to a correlation matrix (ie normalise entry in i,j by uncertainty of bin i * uncertainty in bin j)
        """
        covMatrix = self.covarianceMatrix()
        nbins = len(covMatrix)
        corr = [[0 for i in range(nbins)] for j in range (nbins)]
        for i in range(nbins):
            sigma_i = covMatrix[i][i]
            for j in range(nbins):
                sigma_j = covMatrix[j][j]
                corr[i][j] = covMatrix[i][j] / (sigma_i * sigma_j)**0.5
        return self._mknp(corr)


    def xVals(self):
        return self._mknp([p.x() for p in self.points()])

    def xMins(self):
        """All x low values."""
        # TODO: add extra dimensionality for multiple errors?
        return self._mknp([p.xMin() for p in self.points()])

    def xMaxs(self):
        """All x high values."""
        # TODO: add extra dimensionality for multiple errors?
        return self._mknp([p.xMax() for p in self.points()])

    def xErrs(self):
        """All x error pairs"""
        # TODO: add extra dimensionality for multiple errors?
        return self._mknp([p.xErrs() for p in self.points()])

    def xErrAvgs(self):
        """All x average errors"""
        # TODO: add extra dimensionality for multiple errors?
        return self._mknp([p.xErrAvg() for p in self.points()])

    def xMin(self):
        """Lowest x value."""
        # TODO: add extra dimensionality for multiple errors?
        return min(self.xMins())

    def xMax(self):
        """Highest x value."""
        # TODO: add extra dimensionality for multiple errors?
        return max(self.xMaxs())


    def yVals(self):
        return self._mknp([p.y() for p in self.points()])

    def yMins(self):
        """All y low values."""
        # TODO: add extra dimensionality for multiple errors?
        return self._mknp([p.yMin() for p in self.points()])

    def yMaxs(self):
        """All y high values."""
        # TODO: add extra dimensionality for multiple errors?
        return self._mknp([p.yMax() for p in self.points()])

    def yErrs(self):
        """All y error pairs"""
        # TODO: add extra dimensionality for multiple errors?
        return self._mknp([p.yErrs() for p in self.points()])

    def yErrAvgs(self):
        """All y average errors"""
        # TODO: add extra dimensionality for multiple errors?
        return self._mknp([p.yErrAvg() for p in self.points()])

    def yMin(self):
        """Lowest x value."""
        # TODO: add extra dimensionality for multiple errors?
        return min(self.yMins())

    def yMax(self):
        """Highest y value."""
        # TODO: add extra dimensionality for multiple errors?
        return max(self.yMaxs())


## Convenience alias
S2D = Scatter2D
