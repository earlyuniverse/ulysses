cimport util
cdef class Point2D(Point):
    """
    A 2D point with errors, used by the Scatter2D class.
    """

    cdef c.Point2D* p2ptr(self) except NULL:
        return <c.Point2D*> self.ptr()


    def __init__(self, x=0, y=0, xerrs=0, yerrs=0, source=""):
        if source==None: source=""
        cutil.set_owned_ptr(self, new c.Point2D())
        self.setX(x)
        self.setY(y)
        self.setXErrs(xerrs)
        self.setYErrs(yerrs, source)

    def copy(self):
        return cutil.new_owned_cls(Point2D, new c.Point2D(deref(self.p2ptr())))

    # TODO: add clone() as mapping to (not yet existing) C++ newclone()?

    # LC: Can't get this to work :/
    # use addPoint method instead to
    # ensure parentage is correctly set
    #def setParent(self, Scatter2D scatter):
    #    """Set the parent scatter of this point"""
    #    #return self.p2ptr().setParent(<Scatter*> scatter.s2ptr()[0])
    #    self.p2ptr().setParent(deref(scatter.s2ptr()[0]))
    #    #return None

    def x(self):
        """The x value"""
        return self.p2ptr().x()
    def setX(self, x):
        """Set the x value"""
        self.p2ptr().setX(x)

    def xErrs(self):
        """The x errors"""
        return util.read_error_pair(self.p2ptr().xErrs())

    def setXErrs(self, val):
        """Set the x errors"""
        self.p2ptr().setXErrs(util.read_symmetric(val))

    def xMin(self):
        """The minimum x position, i.e. lowest error"""
        return self.p2ptr().xMin()
    def xMax(self):
        """The maximum x position, i.e. highest error"""
        return self.p2ptr().xMax()

    def xErrAvg(self):
        return self.p2ptr().xErrAvg()



    def y(self):
        """The y value"""
        return self.p2ptr().y()
    def setY(self, y):
        """Set the y value"""
        self.p2ptr().setY(y)

    def yErrs(self):
        """The (total) y errors"""
        return util.read_error_pair(self.p2ptr().yErrs())

    def yErrsFromSource(self, source):
        """The y errors from a particular named source of uncertainty in the breakdown"""
        if isinstance(source, str):
           source = source.encode('utf-8')
        return util.read_error_pair(self.p2ptr().yErrs(source))

    def setYErrs(self, *es):
        """
        (float,) -> None
        ([float, float]) -> None
        (float, float) -> None
        (float, string) -> None
        ([float, float], string) -> None
        (float, float, string) -> None

        Set asymmetric errors on y-axis with an optional string argument to
        specify which named source of uncertainty in the error breakdown should
        be set. By default, if no source is provided, the total uncertainty is set.

        TODO: simplify, this is too much for the Python wrapper
        """
        source = None
        es = list(es)
        if type(es[-1]) is str:
            source = es[-1]
            es = es[:-1]
        errs = es
        if source is None:
            source = ""
        if isinstance(source, str):
            source = source.encode('utf-8')
        if len(errs) == 1:
            if not hasattr(errs[0], "__iter__"):
                self.setErr(2,errs[0], source)
                return
            errs = errs[0]
        # assert len(errs) == 2:
        self.pptr().setErrs(2, tuple(errs), source)

    def setYErrs(self, val, source):
        """(float, string) -> None
        Set symmetric errors on y-axis with an optional string argument to
        specify which named source of uncertainty in the error breakdown should
        be set. By default, if no source is provided, the total uncertainty is set"""
        if source is None:
            source = ""
        self.p2ptr().setYErrs(util.read_symmetric(val), source)

    def yMin(self):
        """The minimum y position, i.e. lowest error"""
        return self.p2ptr().yMin()
    def yMax(self):
        """The maximum y position, i.e. highest error"""
        return self.p2ptr().yMax()

    def yErrAvg(self):
        return self.p2ptr().yErrAvg()



    # property xy:
    #     """x and y coordinates as a tuple"""
    #     def __get__(self):
    #         return util.XY(self.x, self.y)
    #     def __set__(self, val):
    #         self.x, self.y = val




    def scaleX(self, a):
        """(float) -> None
        Scale the x values and errors by factor a."""
        self.p2ptr().scaleX(a)

    def scaleY(self, a):
        """(float) -> None
        Scale the y values and errors by factor a."""
        self.p2ptr().scaleY(a)

    def scaleXY(self, x=1.0, y=1.0):
        """
        (float=1, float=1) -> None
        Scale the point coordinates by the given factors.
        """
        self.p2ptr().scaleXY(x, y)

    # TODO: remove
    def scale(self, x=1.0, y=1.0):
        """
        (float=1, float=1) -> None
        DEPRECATED! Use scaleXY
        Scale the point coordinates by the given factors.
        """
        self.p2ptr().scaleXY(x, y)


    def __repr__(self):
        return '<Point2D(x=%g, y=%g)>' % (self.x(), self.y())

    def __richcmp__(Point2D self, Point2D other, int op):
        if op == 0:
            return deref(self.p2ptr()) < deref(other.p2ptr())
        elif op == 1:
            return deref(self.p2ptr()) <= deref(other.p2ptr())
        elif op == 2:
            return deref(self.p2ptr()) == deref(other.p2ptr())
        elif op == 3:
            return deref(self.p2ptr()) != deref(other.p2ptr())
        elif op == 4:
            return deref(self.p2ptr()) > deref(other.p2ptr())
        elif op == 5:
            return deref(self.p2ptr()) >= deref(other.p2ptr())
