cimport util
cdef class Point(util.Base):
    """
    A generic point with errors, used by the Scatter classes.
    """

    cdef c.Point* pptr(self) except NULL:
        return <c.Point*> self.ptr()

    def __dealloc__(self):
        cdef c.Point *p = self.pptr()
        if self._deallocate:
            del p

    # def __init__(self):
    #     cutil.set_owned_ptr(self, new c.Point())

    # def copy(self):
    #     return cutil.new_owned_cls(Point, new c.Point(deref(self.pptr())))

    # TODO: add clone() as mapping to (not yet existing) C++ newclone()?


    #@property
    def dim(self):
        """None -> int

        Space dimension of the point (should match containing Scatter)"""
        return self.pptr().dim()


    def val(self, i):
        """int -> float
        Value on axis i"""
        return self.pptr().val(i)

    def setVal(self, i, val):
        """(int, float) -> None

        Value on axis i"""
        self.pptr().setVal(i, val)


    def errs(self, i, source=""):
        """int, string -> float

        Errors on axis i. Optional string argument allows to
        access a particular named uncertainty source from the
        breakdown.
        """
        if source is None: source = ""
        if isinstance(source, str):
           source = source.encode('utf-8')
        return util.read_error_pair(self.pptr().errs(i,source))

    def setErr(self, i, e, source=""):
        """(int, float, string) -> None

        Set symmetric errors on axis i, with optional string argument
        to specify which named uncertainty source from the breakdown to
        change. By default, the total uncertainty is changed.
        """
        if source is None: source = ""
        if isinstance(source, str):
           source = source.encode('utf-8')
        self.pptr().setErr(i, e, source)

    def setErrs(self, i, *es):
        """
        (int, float) -> None
        (int, [float, float]) -> None
        (int, float, float) -> None
        (int, float, string) -> None
        (int, [float, float], string) -> None
        (int, float, float, string) -> None

        Set asymmetric errors on axis i.

        Optional string argument to specify which named uncertainty
        source from the breakdown to change. By default, the total
        uncertainty is changed.
        """
        source=None
        es=list(es)
        if type(es[-1]) is str:
            source=es[-1]
            es=es[:-1]
        else:
            try:
                source = es[-1].decode('utf-8')
                es=es[:-1]
            except:
                pass
        errs = es
        if source is None: source=""
        if len(errs) == 1:
            if not hasattr(errs[0], "__iter__"):
                self.setErr(i,errs[0], source)
                return
            errs=errs[0]
        # assert len(errs) == 2:
        if isinstance(source, str):
             source = source.encode('utf-8')
        self.pptr().setErrs(i, tuple(errs), source)


    def errMinus(self, i, source=""):
        """int, string -> float

        Minus error on axis i.

        Optional string argument to specify which named uncertainty
        source from the breakdown to change. By default, the total
        uncertainty is changed.
        """
        if source is None: source = ""
        if isinstance(source, str):
           source = source.encode('utf-8')
        return self.pptr().errMinus(i ,source)

    def setErrMinus(self, i, e, source=""):
        """(int, float, string) -> None

        Set minus error on axis i.

        Optional string argument to specify which named uncertainty
        source from the breakdown to change. By default, the total
        uncertainty is changed.
        """
        if source is None: source = ""
        if isinstance(source, str):
           source = source.encode('utf-8')
        self.pptr().setErrMinus(i, e, source)


    def errPlus(self, i, source=""):
        """int, string -> float

        Plus error on axis i.

        Optional string argument to specify which named uncertainty
        source from the breakdown to access. By default, the total
        uncertainty is accessed.
        """
        if source is None: source = ""
        if isinstance(source, str):
           source = source.encode('utf-8')
        return self.pptr().errPlus(i, source)

    def setErrPlus(self, i, e, source=""):
        """(int, float, string) -> None

        Set plus error on axis i.

        Optional string argument to specify which named uncertainty
        source from the breakdown to change. By default, the total
        uncertainty is changed.
        """
        if source is None: source = ""
        if isinstance(source, str):
           source = source.encode('utf-8')
        self.pptr().setErrPlus(i, e, source)


    def errAvg(self, i, source=""):
        """int -> float

        Average error on axis i.

        Optional string argument to specify which named uncertainty
        source from the breakdown to access. By default, the total
        uncertainty is accessed.
        """
        if source is None: source = ""
        if isinstance(source, str):
           source = source.encode('utf-8')
        return self.pptr().errAvg(i, source)


    def set(self, i, val, *es, source=""):
        """
        (int, float, float) -> None
        (int, float, [float, float]) -> None
        (int, float, float, float) -> None
        (int, float, float, string) -> None
        (int, float, [float, float],string) -> None
        (int, float, float, float,string) -> None

        Set value and errors on axis i.

        Optional string argument to specify which named uncertainty
        source from the breakdown to change. By default, the total
        uncertainty is changed.
        """
        errs = es
        if source is None: source = ""
        if len(es) == 1:
            if hasattr(es[0], "__iter__"):
                errs = [es[0], es[0]]
            else:
                errs = es[0]
        # assert len(errs) == 2:
        if isinstance(source, str):
           source = source.encode('utf-8')
        self.pptr().set(i, val, errs, source)

    def errMap(self):
        """None -> {string: [float,float]}

        Error map of this point
        """
        err_map = self.pptr().errMap()
        err_map = dict((k.decode('utf-8'), v) for k, v in err_map.items())
        return err_map


    def rmVariations(self):
        """None -> None

        Remove variations in the error map of this point
        """
        return self.pptr().rmVariations()


    def scale(self, i, scale):
        """(int, float) -> None

        Scale values on axis i
        """
        self.pptr().scale(i, scale)


    # def __repr__(self):
    #     return '<Point(x=%g)>' % self.x
