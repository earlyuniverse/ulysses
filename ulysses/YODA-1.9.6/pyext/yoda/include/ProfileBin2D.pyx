#TODO improve this once we have a working Profile2D
cdef class ProfileBin2D(Bin2D_Dbn3D):

    cdef inline c.ProfileBin2D* pb2ptr(self) except NULL:
        return <c.ProfileBin2D*> self.ptr()
    # TODO: remove
    cdef inline c.ProfileBin2D* _ProfileBin2D(self) except NULL:
        return <c.ProfileBin2D*> self.ptr()

    def __init__(self, xlow, xhigh, ylow, yhigh):
        cutil.set_owned_ptr(self, new c.ProfileBin2D(xlow, xhigh, ylow, yhigh))

    # def fill(self, x, y, z, weight=1.0, fraction=1.0):
    #     self.pb2ptr().fill(x, y, z, weight, fraction)

    # def fill_bin(self, z, weight=1.0, fraction=1.0):
    #     self.pb2ptr().fillBin(z, weight, fraction)

    #@property
    def mean(self):
        return self.pb2ptr().mean()

    #@property
    def stdDev(self):
        return self.pb2ptr().stdDev()

    #@property
    def variance(self):
        return self.pb2ptr().variance()

    #@property
    def stdErr(self):
        return self.pb2ptr().stdErr()

    #@property
    def rms(self):
        return self.pb2ptr().rms()

    #@property
    def sumWZ(self):
        return self.pb2ptr().sumWZ()

    #@property
    def sumWZ2(self):
        return self.pb2ptr().sumWZ2()

    def __add__(ProfileBin2D a, ProfileBin2D b):
        return cutil.new_owned_cls(ProfileBin2D, new c.ProfileBin2D(deref(a.pb2ptr()) + deref(b.pb2ptr())))

    def __sub__(ProfileBin2D a, ProfileBin2D b):
        return cutil.new_owned_cls(ProfileBin2D, new c.ProfileBin2D(deref(a.pb2ptr()) - deref(b.pb2ptr())))

    def __repr__(self):
        # TODO: Urk
        return 'ProfileBin2D(%g, %g, %g, %g)' % (self.edges.x + self.edges.y)
