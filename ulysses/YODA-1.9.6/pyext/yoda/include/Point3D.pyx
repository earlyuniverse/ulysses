cimport util
cdef class Point3D(Point):
    """
    A 3D point with errors, used by the Scatter3D class.
    """

    cdef c.Point3D* p3ptr(self) except NULL:
        return <c.Point3D*> self.ptr()


    def __init__(self, x=0, y=0, z=0, xerrs=0, yerrs=0, zerrs=0, source=""):
        if source==None: source=""
        cutil.set_owned_ptr(self, new c.Point3D())
        # TODO: need shortcuts
        self.setX(x)
        self.setY(y)
        self.setZ(z)
        self.setXErrs(xerrs)
        self.setYErrs(yerrs)
        self.setZErrs(zerrs, source)

    def copy(self):
        return cutil.new_owned_cls(Point3D, new c.Point3D(deref(self.p3ptr())))



    def x(self):
        """The x value"""
        return self.p3ptr().x()
    def setX(self, x):
        """Set the x value"""
        self.p3ptr().setX(x)

    def xErrs(self):
        """The x errors"""
        return util.read_error_pair(self.p3ptr().xErrs())

    def setXErrs(self, val):
        """Set the x errors"""
        self.p3ptr().setXErrs(util.read_symmetric(val))

    def xMin(self):
        """The minimum x position, i.e. lowest error"""
        return self.p3ptr().xMin()
    def xMax(self):
        """The maximum x position, i.e. highest error"""
        return self.p3ptr().xMax()

    def xErrAvg(self):
        return self.p3ptr().xErrAvg()


    def y(self):
        """The y value"""
        return self.p3ptr().y()
    def setY(self, y):
        """Set the y value"""
        self.p3ptr().setY(y)

    def yErrs(self):
        """The y errors"""
        return util.read_error_pair(self.p3ptr().yErrs())

    def setYErrs(self, val):
        """Set the y errors"""
        self.p3ptr().setYErrs(util.read_symmetric(val))

    def yMin(self):
        """The minimum y position, i.e. lowest error"""
        return self.p3ptr().yMin()
    def yMax(self):
        """The maximum y position, i.e. highest error"""
        return self.p3ptr().yMax()

    def yErrAvg(self):
        return self.p3ptr().yErrAvg()


    def z(self):
        """The z value"""
        return self.p3ptr().z()
    def setZ(self, z):
        """Set the z value"""
        self.p3ptr().setZ(z)

    def zErrs(self):
        """The z errors"""
        return util.read_error_pair(self.p3ptr().zErrs())

    def zErrsFromSource(self, source):
        """The z errors"""
        if isinstance(source, str):
           source = source.encode('utf-8')
        return util.read_error_pair(self.p3ptr().zErrs(source))

    def setZErrs(self, *es):
        """
        (float,) -> None
        ([float, float]) -> None
        (float, float) -> None
        (float, string) -> None
        ([float, float], string) -> None
        (float, float, string) -> None

        Set asymmetric errors on z-axis with an optional string argument to
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
                self.setErr(3, errs[0], source)
                return
            errs = errs[0]
        # assert len(errs) == 2:
        self.pptr().setErrs(3, tuple(errs), source)

    def setZErrs(self, val, source):
        """(float, string) -> None
        Set symmetric errors on z-axis with an optional string argument to
        specify which named source of uncertainty in the error breakdown should
        be set. By default, if no source is provided, the total uncertainty is set"""
        if source is None:
            source = ""
        self.p3ptr().setZErrs(util.read_symmetric(val), source)

    def zMin(self):
        """The minimum z position, i.e. lowest error"""
        return self.p3ptr().zMin()
    def zMax(self):
        """The maximum z position, i.e. highest error"""
        return self.p3ptr().zMax()

    def zErrAvg(self):
        return self.p3ptr().zErrAvg()


    # property xyz:
    #     def __get__(self):
    #         return util.XYZ(self.x, self.y, self.z)
    #     def __set__(self, val):
    #         self.x, self.y, self.z = val


    # # TODO: How does this fit into the multi-error API? Still useful, but just reports first errs... how to get _all_ +- err pairs?
    # # LC: This is fine because preserntly only the highest dimension supports multi-errors
    # property xErrs:
    #     def __get__(self):
    #         return util.read_error_pair(self.p3ptr().xErrs())
    #     def __set__(self, val):
    #         self.p3ptr().setXErrs(util.read_symmetric(val))

    # # TODO: How does this fit into the multi-error API? Still useful, but just reports first errs... how to get _all_ +- err pairs?
    # # LC: This is fine because preserntly only the highest dimension supports multi-errors
    # property yErrs:
    #     def __get__(self):
    #         return util.read_error_pair(self.p3ptr().yErrs())
    #     def __set__(self, val):
    #         self.p3ptr().setYErrs(util.read_symmetric(val))

    # # TODO: How does this fit into the multi-error API? Still useful, but just reports first errs... how to get _all_ +- err pairs?
    # # LC: I think it's Ok to leave this like this, for most users the nominal is what they want anyway,
    # # and for those who want multi-errs, they can set using a method eg setErrs(dim,(ed,eu),source) and access using errs(dim,(ed,eu),source)
    # property zErrs:
    #     def __get__(self):
    #         return util.read_error_pair(self.p3ptr().zErrs())
    #     def __set__(self, val):
    #         self.p3ptr().setZErrs(util.read_symmetric(val))



    def scaleX(self, ax):
        """
        (float) -> None
        Scale the x point coordinates by the given factor.
        """
        self.p3ptr().scaleX(ax)

    def scaleY(self, ay):
        """
        (float) -> None
        Scale the y point coordinates by the given factor.
        """
        self.p3ptr().scaleY(ay)

    def scaleZ(self, az):
        """
        (float) -> None
        Scale the z point coordinates by the given factor.
        """
        self.p3ptr().scaleZ(az)

    def scaleXYZ(self, ax=1.0, ay=1.0, az=1.0):
        """
        (float=1.0, float=1.0, float=1.0) -> None
        Scale the point coordinates by the given factors.
        """
        self.p3ptr().scaleXYZ(ax, ay, az)

    # TODO: remove
    def scaleXYZ(self, ax=1.0, ay=1.0, az=1.0):
        """
        (double=1.0, double=1.0, double=1.0) -> None
        DEPRECATED: USE scaleXYZ
        Scale the point coordinates by the given factors.
        """
        self.scaleXYZ(ax, ay, az)


    # TODO: transformX,Y,Z


    def __repr__(self):
        return '<Point3D(x=%g, y=%g, z=%g)>' % (self.x(), self.y(), self.z())

    def __richcmp__(Point3D self, Point3D other, int op):
        if op == 0:
            return deref(self.p3ptr()) < deref(other.p3ptr())
        elif op == 1:
            return deref(self.p3ptr()) <= deref(other.p3ptr())
        elif op == 2:
            return deref(self.p3ptr()) == deref(other.p3ptr())
        elif op == 3:
            return deref(self.p3ptr()) != deref(other.p3ptr())
        elif op == 4:
            return deref(self.p3ptr()) > deref(other.p3ptr())
        elif op == 5:
            return deref(self.p3ptr()) >= deref(other.p3ptr())
