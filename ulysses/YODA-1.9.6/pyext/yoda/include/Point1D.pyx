cimport util
cdef class Point1D(Point):
    """
    A 1D point with errors, used by the Scatter1D class.
    """

    cdef c.Point1D* p1ptr(self) except NULL:
        return <c.Point1D*> self.ptr()


    def __init__(self, x=0, xerrs=0, source=""):
        if source==None: source=""
        cutil.set_owned_ptr(self, new c.Point1D())
        self.setX(x)
        self.setXErrs(xerrs, source)

    def copy(self):
        return cutil.new_owned_cls(Point1D, new c.Point1D(deref(self.p1ptr())))

    # TODO: add clone() as mapping to (not yet existing) C++ newclone()?


    def setXErrs(self, *es):
        """
        (float,) -> None
        ([float, float]) -> None
        (float, float) -> None
        (float, string) -> None
        ([float, float], string) -> None
        (float, float, string) -> None

        Set asymmetric errors on x-axis with an optional string argument to
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
                self.setErr(1, errs[0], source) #< duplicate of last line?
                return
            errs = errs[0]
        # assert len(errs) == 2:
        self.pptr().setErrs(1, tuple(errs), source)

    def setXErrs(self, val, source):
        """(float, string) -> None

        Set symmetric errors on x-axis with an optional string argument to
        specify which named source of uncertainty in the error breakdown should
        be set. By default, if no source is provided, the total uncertainty is set"""
        if source is None:
            source = ""
        self.p1ptr().setXErrs(util.read_symmetric(val), source)


    # property x:
    #     """x coordinate"""
    #     def __get__(self):
    #         return self.p1ptr().x()
    #     def __set__(self, x):
    #         self.p1ptr().setX(x)

    # property xErrs:
    #     """The x errors"""
    #     def __get__(self):
    #         return util.read_error_pair(self.p1ptr().xErrs())
    #     def __set__(self, val):
    #         self.p1ptr().setXErrs(util.read_symmetric(val))

    def x(self):
        """The x value"""
        return self.p1ptr().x()
    def setX(self, x):
        """Set the x value"""
        self.p1ptr().setX(x)

    def xErrs(self):
        """The x errors"""
        return util.read_error_pair(self.p1ptr().xErrs())

    def xErrsFromSource(self, source):
        """The y errors"""
        if isinstance(source, str):
           source = source.encode('utf-8')
        return util.read_error_pair(self.p1ptr().xErrs(source))

    def setXErrs(self, *es):
        """(int, float) -> None
           (int, [float, float]) -> None
           (int, float, float) -> None
        Set asymmetric errors on axis i"""
        source = None
        es = list(es)
        if type(es[-1]) is str:
            source = es[-1]
            es = es[:-1]
        else:
            pass
        errs = es
        if source is None:
            source = ""
        if len(errs) == 1:
            if not hasattr(errs[0], "__iter__"):
                self.setErr(1,errs[0], source)
                return
            errs = errs[0]
        # assert len(errs) == 2:
        if isinstance(source, str):
           source = source.encode('utf-8')
        self.pptr().setErrs(1, tuple(errs), source)

    def setYErrs(self, val, source):
        if source is None:
            source = ""
        self.p1ptr().setXErrs(util.read_symmetric(val))

    #@property
    def xMin(self):
        """The minimum x position, i.e. lowest error"""
        return self.p1ptr().xMin()
    #@property
    def xMax(self):
        """The maximum x position, i.e. highest error"""
        return self.p1ptr().xMax()

    def xErrAvg(self):
        return self.p1ptr().xErrAvg()


    def scaleX(self, a):
        """(float) -> None
        Scale the x values and errors by factor a."""
        self.p1ptr().scaleX(a)


    def __repr__(self):
        return '<Point1D(x=%g)>' % self.x()

    def __richcmp__(Point1D self, Point1D other, int op):
        if op == 0:
            return deref(self.p1ptr()) < deref(other.p1ptr())
        elif op == 1:
            return deref(self.p1ptr()) <= deref(other.p1ptr())
        elif op == 2:
            return deref(self.p1ptr()) == deref(other.p1ptr())
        elif op == 3:
            return deref(self.p1ptr()) != deref(other.p1ptr())
        elif op == 4:
            return deref(self.p1ptr()) > deref(other.p1ptr())
        elif op == 5:
            return deref(self.p1ptr()) >= deref(other.p1ptr())
