cimport util
cdef class Scatter3D(AnalysisObject):
    """
    3D scatter plot, i.e. a collection of Point3D objects with positions and errors.

    Constructor calling idioms:

    Scatter3D(path="", title="")
      Create a new empty scatter, with optional path and title.

    Scatter3D(points, path="", title=""):
      Create a new empty scatter from an iterable of points, with optional path
      and title.

    TODO: more documentation!
    """

    cdef inline c.Scatter3D* s3ptr(self) except NULL:
        return <c.Scatter3D*> self.ptr()

    def __init__(self, *args, **kwargs):
        util.try_loop([self.__init_2, self.__init_3], *args, **kwargs)

    def __init_2(self, path="", title=""):
        path  = path.encode('utf-8')
        title = title.encode('utf-8')
        cutil.set_owned_ptr(self, new c.Scatter3D(<string>path, <string>title))

    def __init_3(self, points, char* path="", char* title=""):
        self.__init_2(path, title)
        self.addPoints(points)

    def clone(self):
        """() -> Scatter3D.
        Clone this Scatter3D."""
        return cutil.new_owned_cls(Scatter3D, self.s3ptr().newclone())

    def __repr__(self):
        return "<%s '%s' %d points>" % (self.__class__.__name__, self.path(), len(self.points()))


    def reset(self):
        "Reset the scatter, removing all points"
        self.s3ptr().reset()


    def numPoints(self):
        """() -> int
        Number of points in this scatter."""
        return self.s3ptr().numPoints()

    def __len__(self):
        return self.numPoints()


    def points(self):
        """Access the ordered list of points."""
        return [self.point(i) for i in range(self.numPoints())]

    def point(self, size_t i):
        """Access the i'th point."""
        return cutil.new_borrowed_cls(Point3D, &self.s3ptr().point(i), self)

    def __getitem__(self, py_ix):
        cdef size_t i = cutil.pythonic_index(py_ix, self.numPoints())
        return cutil.new_borrowed_cls(Point3D, &self.s3ptr().point(i), self)


    def addPoint(self, *args, **kwargs):
        """Add a new point.

        Provide either a single yoda.Point3D object, or the
        3-6 args: x, y, z, xerrs=0, yerrs=0, zerrs=0.
        """
        try:
            self.__addPoint_point(*args, **kwargs)
        except TypeError:
            self.__addPoint_explicit(*args, **kwargs)

    def __addPoint_explicit(self, x, y, z, xerrs=0, yerrs=0, zerrs=0):
        self.__addPoint_point(Point3D(x, y, z, xerrs, yerrs, zerrs))

    def __addPoint_point(self, Point3D p):
        self.s3ptr().addPoint(p.p3ptr()[0])

    def addPoints(self, iterable):
        """Add several new points."""
        for row in iterable:
          try:
            self.addPoint(*row)
          except TypeError:
            self.addPoint(row)

    def rmPoint(self, idx):
        self.s3ptr().rmPoint(idx)

    def rmPoints(self, idxs):
        self.s3ptr().rmPoints(idxs)

    def combineWith(self, others):
        """Try to add points from other Scatter3Ds into this one."""
        cdef Scatter3D other
        try:
            # Can we type it as a Scatter3D?
            other = others
        except TypeError:
            # Could be an iterable...
            for other in others:
                self.s3ptr().combineWith(deref(other.s3ptr()))
        else:
            self.s3ptr().combineWith(deref(other.s3ptr()))


    def mkScatter(self):
        """None -> Scatter3D.
        Make a new Scatter3D. Exists to allow mkScatter calls on any AnalysisObject,
        even if it already is a scatter."""
        cdef c.Scatter3D s3 = c.mkScatter_Scatter3D(deref(self.s3ptr()))
        return cutil.new_owned_cls(Scatter3D, s3.newclone())


    def scaleX(self, a):
        """(float) -> None
        Scale the x values and errors of the points in this scatter by factor a."""
        self.s3ptr().scaleX(a)

    def scaleY(self, a):
        """(float) -> None
        Scale the y values and errors of the points in this scatter by factor a."""
        self.s3ptr().scaleY(a)

    def scaleZ(self, a):
        """(float) -> None
        Scale the z values and errors of the points in this scatter by factor a."""
        self.s3ptr().scaleZ(a)

    def scaleXYZ(self, ax=1, ay=1, az=1):
        """(float=1, float=1, float=1) -> None
        Scale the values and errors of the points in this scatter by factors ax, ay, az."""
        self.s3ptr().scaleXYZ(ax, ay, az)

    def scale(self, i, scale):
        """(int, float) -> None
        Scale values on axis i"""
        self.s3ptr().scale(i, scale)


    def transformX(self, f):
        """(fn) -> None
        Transform the x values and errors of the points in this scatter by function f."""
        import ctypes
        try:
            callback = ctypes.CFUNCTYPE(ctypes.c_double, ctypes.c_double)(f)
        except:
            raise RuntimeError("Callback is not of type (double) -> double")
        fptr = (<c.dbl_dbl_fptr*><size_t>ctypes.addressof(callback))[0]
        c.Scatter3D_transformX(deref(self.s3ptr()), fptr)

    def transformY(self, f):
        """(fn) -> None
        Transform the y values and errors of the points in this scatter by function f."""
        import ctypes
        try:
            callback = ctypes.CFUNCTYPE(ctypes.c_double, ctypes.c_double)(f)
        except:
            raise RuntimeError("Callback is not of type (double) -> double")
        fptr = (<c.dbl_dbl_fptr*><size_t>ctypes.addressof(callback))[0]
        c.Scatter3D_transformY(deref(self.s3ptr()), fptr)

    def transformZ(self, f):
        """(fn) -> None
        Transform the z values and errors of the points in this scatter by function f."""
        import ctypes
        try:
            callback = ctypes.CFUNCTYPE(ctypes.c_double, ctypes.c_double)(f)
        except:
            raise RuntimeError("Callback is not of type (double) -> double")
        fptr = (<c.dbl_dbl_fptr*><size_t>ctypes.addressof(callback))[0]
        c.Scatter3D_transformZ(deref(self.s3ptr()), fptr)


    # # TODO: remove?
    # def __add__(Scatter3D self, Scatter3D other):
    #     return cutil.new_owned_cls(Scatter3D, c.Scatter3D_add_Scatter3D(self.s3ptr(), other.s3ptr()))

    # # TODO: remove?
    # def __sub__(Scatter3D self, Scatter3D other):
    #     return cutil.new_owned_cls(Scatter3D, c.Scatter3D_sub_Scatter3D(self.s3ptr(), other.s3ptr()))
    
    def updateTotalUncertainty(self):
        """None -> None
        Sum the error in the error map in quadrature and update the total ("") uncertainty.
        """
        return self.s3ptr().updateTotalUncertainty()
    
    def writeVariationsToAnnotations(self):
        """None -> None
        Parse the variations and add them to the annotations.
        """
        return self.s3ptr().writeVariationsToAnnotations()

    def variations(self):
        """None -> vector[string]
        Get the list of variations stored in the points of the Scatter"""
        vars = self.s3ptr().variations()
        return [var.decode('utf-8') for var in vars]
    
    def rmVariations(self):
        """None -> None
        Remove the variations stored in the points of the Scatter"""
        return self.s3ptr().rmVariations()
    
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

    def covarianceMatrix(self, *ignoreOffDiagonalTerms):
        """bool -> vector[vector[float]]
        Construct the covariance matrix"""
        return self._mknp(self.s2ptr().covarianceMatrix(ignoreOffDiagonalTerms))


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

    def xErrs(self):
        """All x errors"""
        return self._mknp([p.xErrs() for p in self.points()])

    def xMin(self):
        """Lowest x value."""
        return min(self.xMins())

    def xMax(self):
        """Highest x value."""
        return max(self.xMaxs())


    def yVals(self):
        return self._mknp([p.y() for p in self.points()])

    def yMins(self):
        """All x low values."""
        return self._mknp([p.yMin() for p in self.points()])

    def yMaxs(self):
        """All x high values."""
        return self._mknp([p.yMax() for p in self.points()])

    def yErrs(self):
        """All y errors"""
        return self._mknp([p.yErrs() for p in self.points()])

    def yMin(self):
        """Lowest x value."""
        return min(self.yMins())

    def yMax(self):
        """Highest y value."""
        return max(self.yMaxs())


    def zVals(self):
        return self._mknp([p.z() for p in self.points()])

    def zMins(self):
        """All z low values."""
        return self._mknp([p.zMin() for p in self.points()])

    def zMaxs(self):
        """All z high values."""
        return self._mknp([p.zMax() for p in self.points()])

    def zErrs(self):
        """All z errors"""
        return self._mknp([p.zErrs() for p in self.points()])

    def zMin(self):
        """Lowest z value."""
        return min(self.zMins())

    def zMax(self):
        """Highest z value."""
        return max(self.zMaxs())


## Convenience alias
S3D = Scatter3D
