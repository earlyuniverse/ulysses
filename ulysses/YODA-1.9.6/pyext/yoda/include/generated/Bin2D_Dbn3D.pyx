#################################################
#############       WARNING      ################
#################################################
# This file has been automatically generated.
# Any changes you make here will get overridden.
# Instead, make your changes in Bin2D_DBN.pyx
#################################################
cimport util
# TODO: docstrings
cdef class Bin2D_Dbn3D(Bin):
    """2D Bin class templated on a Dbn3D"""

    cdef inline c.Bin2D_Dbn3D* b2ptr(self) except NULL:
        return <c.Bin2D_Dbn3D*> self.ptr()
    # TODO: remove
    cdef inline c.Bin2D_Dbn3D* _Bin2D(self) except NULL:
        return <c.Bin2D_Dbn3D*> self.ptr()


    def __init__(self, xlow, xhigh, ylow, yhigh):
        cutil.set_owned_ptr(self, new c.Bin2D_Dbn3D( pair[double, double](xlow, xhigh),
                                                      pair[double, double](ylow, yhigh) ))

    def __repr__(self):
        return '<%s x=[%g, %g), y=[%g, %g)>' % (self.__class__.__name__, self.xMin(), self.xMax(), self.yMin(), self.yMax())


    # def scaleXY(self, x=1.0, y=1.0):
    #     self.b2ptr().scaleXY(x, y)

    # def scaleW(self, w):
    #     self.b2ptr().scaleW(w)


    def xEdges(self):
        """
        The lower and upper x edges.
        """
        cdef pair[double, double] x = self.b2ptr().xEdges()
        return (x.first, x.second)

    def yEdges(self):
        """
        The lower and upper y edges.
        """
        cdef pair[double, double] y = self.b2ptr().yEdges()
        return (y.first, y.second)

    def xyEdges(self):
        """
        The lower and upper x,y edge pairs.
        """
        return util.XY(self.xEdges, self.yEdges)


    def xMin(self):
        """Low edge in x."""
        return self.b2ptr().xMin()

    def yMin(self):
        """Low edge in y."""
        return self.b2ptr().yMin()

    def xyMin(self):
        """Low edges in x,y."""
        return util.XY(self.xMin(), self.yMin())


    def xMax(self):
        """High edge in x."""
        return self.b2ptr().xMax()

    def yMax(self):
        """High edge in y."""
        return self.b2ptr().yMax()

    def xyMax(self):
        """High edges in x,y."""
        return util.XY(self.xMax, self.yMax)


    def xMid(self):
        """Geometric centre of the bin in x"""
        return self.b2ptr().xMid()

    def yMid(self):
        """Geometric centre of the bin in y"""
        return self.b2ptr().yMid()

    def xyMid(self):
        """Geometric centre of the bin"""
        return util.XY(self.xMid, self.yMid)


    def xWidth(self):
        """Width of the bin in x"""
        return self.b2ptr().xWidth()

    def yWidth(self):
        """Width of the bin in y"""
        return self.b2ptr().yWidth()

    def xyWidths(self):
        """The widths of this bin in the x- and y-dimensions."""
        return util.XY(self.xWidth, self.yWidth)


    def area(self):
        """The area of this bin in the x-y plane."""
        return self.b2ptr().area()


    def xFocus(self):
        """Focus of the bin in x"""
        return self.b2ptr().xFocus()

    def yFocus(self):
        """Focus of the bin in y"""
        return self.b2ptr().yFocus()

    def xyFocus(self):
        """The focus of the bin in the x- and y-dimensions"""
        return util.XY(self.xFocus, self.yFocus)


    def xMean(self):
        return self.b2ptr().xMean()

    def yMean(self):
        return self.b2ptr().xMean()

    def xyMean(self):
        return util.XY(self.xMean, self.yMean)


    def xVariance(self):
        return self.b2ptr().xVariance()

    def yVariance(self):
        return self.b2ptr().xVariance()

    def xyVariance(self):
        return util.XY(self.xVariance, self.yVariance)


    def xStdDev(self):
        return self.b2ptr().xStdDev()

    def yStdDev(self):
        return self.b2ptr().yStdDev()

    def xyStdDev(self):
        return util.XY(self.xStdDev, self.yStdDev)


    def xStdErr(self):
        return self.b2ptr().xStdErr()

    def yStdErr(self):
        return self.b2ptr().yStdErr()

    def xyStdErr(self):
        return util.XY(self.xStdErr, self.yStdErr)


    def xRMS(self):
        return self.b2ptr().xRMS()

    def yRMS(self):
        return self.b2ptr().yRMS()

    def xyRMS(self):
        return util.XY(self.xRMS, self.yRMS)


    # Raw statistics #
    ##################

    def sumWX(self):
        return self.b2ptr().sumWX()

    def sumWY(self):
        return self.b2ptr().sumWY()

    def sumWXY(self):
        return self.b2ptr().sumWXY()

    def sumWX2(self):
        return self.b2ptr().sumWX2()

    def sumWY2(self):
        return self.b2ptr().sumWY2()


    #def merge(Bin2D_Dbn3D self, Bin2D_Dbn3D other):
    #    self.b2ptr().merge(deref(other.b2ptr()))
    #    return self

    def adjacentTo(Bin2D_Dbn3D self, Bin2D_Dbn3D other):
        return self.b2ptr().adjacentTo(deref(other.b2ptr()))


    def __add__(Bin2D_Dbn3D self, Bin2D_Dbn3D other):
        return cutil.new_owned_cls(
            Bin2D_Dbn3D,
            new c.Bin2D_Dbn3D(deref(self.b2ptr()) + deref(other.b2ptr())))

    def __sub__(Bin2D_Dbn3D self, Bin2D_Dbn3D other):
        return cutil.new_owned_cls(
            Bin2D_Dbn3D,
            new c.Bin2D_Dbn3D(deref(self.b2ptr()) - deref(other.b2ptr())))
