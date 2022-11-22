cimport util

cdef class AnalysisObject(util.Base):
    """
    AnalysisObject is the base class of the main user-facing objects, such as
    the Histo, Profile and Scatter classes.
    """

    # Pointer upcasting mechanism
    cdef inline c.AnalysisObject* aoptr(self) except NULL:
        return <c.AnalysisObject*> self.ptr()

    # # Pointer upcasting mechanism
    # # DEPRECATED
    # cdef inline c.AnalysisObject* _AnalysisObject(self) except NULL:
    #     return <c.AnalysisObject*> self.ptr()

    # Deallocator (only needed as a base class)
    def __dealloc__(self):
        p = self.aoptr()
        if self._deallocate:
            del p


    #@property
    def type(self):
        "String identifier for this type"
        return self.aoptr().type().decode('utf-8')

    #@property
    def dim(self):
        "Fill dimension or plot dimension of this object, for fillables and scatters respectively"
        return self.aoptr().dim()

    #@property
    def annotations(self):
        """() -> list[str]
        A list of all annotation/metadata keys."""
        return [ a.decode('utf-8') for a in self.aoptr().annotations() ]

    #@property
    def annotationsDict(self):
        """() -> dict[str->str]
        A dict of all annotations/metadata entries."""
        # TODO: add a map equivalent to C++?
        return dict((k.lower(), self.annotation(k)) for k in self.annotations())

    def annotation(self, k, default=None):
        """Get annotation k from this object (falling back to default if not set).

        The annotation string will be automatically converted to Python
        native types as far as possible -- more complex types are possible
        via the ast and yaml modules."""
        try:
            rtn = self.aoptr().annotation(<string>k.encode('utf-8')).decode('utf-8')
            try:
                import yaml
                rtn = yaml.full_load(rtn)
            except:
                rtn = util._autotype(rtn, True)
        except:
            rtn = default
        return rtn

    def setAnnotation(self, k, v):
        """Set annotation k on this object."""
        self.aoptr().setAnnotation(<string>k.encode('utf-8'),
                                   <string>util._autostr(v).encode('utf-8'))

    def hasAnnotation(self, k):
        """Check if this object has annotation k."""
        return self.aoptr().hasAnnotation(<string>k.encode('utf-8'))

    def rmAnnotation(self, k):
        """Remove annotation k from this object."""
        self.aoptr().rmAnnotation(<string>k.encode('utf-8'))

    def clearAnnotations(self):
        """Clear the annotations dictionary."""
        self.aoptr().clearAnnotations()


    def dump(self):
        """A human readable representation of this object."""
        try:
            from cStringIO import StringIO
        except ImportError:
            from io import StringIO
        f = StringIO()
        writeFLAT([self], f)
        f.seek(0)
        return f.read().strip()


    #@property
    def name(self):
        """
        Return the histogram name, i.e. the last part of the path (which may be empty).
        """
        return self.aoptr().name().decode('utf-8')


    def path(self):
        """
        Used for persistence and as a unique identifier. Must begin with
        a '/' if not the empty string.
        """
        return self.aoptr().path().decode('utf-8')

    def setPath(self, path):
        """
        Used for persistence and as a unique identifier. Must begin with
        a '/' if not the empty string.
        """
        self.aoptr().setPath(<string>path.encode('utf-8'))

    # property path:
    #     """
    #     Used for persistence and as a unique identifier. Must begin with
    #     a '/' if not the empty string.
    #     """
    #     def __get__(self):
    #         return self.aoptr().path().decode('utf-8')

    #     def __set__(self, path):
    #         self.aoptr().setPath(<string>path.encode('utf-8'))


    def title(self):
        """
        Histogram title
        """
        return self.aoptr().title().decode('utf-8')

    def setTitle(self, title):
        """
        Set the histogram title (optional)
        """
        self.aoptr().setTitle(<string>title.encode('utf-8'))

    # property title:
    #     """
    #     Convenient access to the histogram title (optional).
    #     """
    #     def __get__(self):
    #         return self.aoptr().title().decode('utf-8')

    #     def __set__(self, title):
    #         self.aoptr().setTitle(<string>title.encode('utf-8'))


    def __repr__(self):
        return "<%s '%s'>" % (self.__class__.__name__, self.path)


## Convenience alias
AO = AnalysisObject
