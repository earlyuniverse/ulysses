def mkScatter(ao, usefocus=False, p_usestddev=False, h_binsizediv=True, uflow_binwidth=-1., oflow_binwidth=-1.):
    """AnalysisObject -> Scatter{1,2,3}D
    Convert an AnalysisObject to a Scatter, using the logic of the bound mkScatter methods.

    All args other than the AO itself should be supplied as keywords rather than
    positional, to avoid trouble since several (prefixed with a target-type
    letter) only apply to specific types of supplied AO.
    """
    s = None
    if ao.type() == "Histo1D":
        s = ao.mkScatter(usefocus, h_binsizediv, uflow_binwidth, oflow_binwidth)
    elif ao.type() == "Histo2D":
        s = ao.mkScatter(usefocus, h_binsizediv)

    elif ao.type() == "Profile1D":
        s = ao.mkScatter(usefocus, p_usestddev, uflow_binwidth, oflow_binwidth)
    elif ao.type() == "Profile2D":
        s = ao.mkScatter(usefocus, p_usestddev)

    else: # Counter and Scatters
        s = ao.mkScatter()

    return s


def divide(ao1, ao2):
    """(AnalysisObject, AnalysisObject) -> Scatter{1,2,3}D
    Divide one AnalysisObject by another, producing a Scatter of appropriate dimension by using the logic of the bound divideBy methods."""
    return ao1.divideBy(ao2)


# def divide(ao1, ao2, efficiency=False):
#     """(AO, AO, [bool]) -> Scatter{1,2,3}D
#     Convert two AnalysisObjects (either two HistoNDs or two ProfileNDs) to a Scatter(N+1)D.

#     The efficiency boolean is used to enable correlated division for Histos, using binomial
#     statistics as appropriate when the numerator contains a subset of the fills in the denominator.

#     TODO: Implement efficiency division mapping.
#     TODO: Implement 2D/2D division mapping.
#     TODO: Implement scatter division.
#     """
#     cdef c.Scatter2D s2
#     if type(ao1) is not type(ao2):
#         raise ValueError("Histograms must be of the same type to be divided")
#     if type(ao1) is Histo1D:
#         s2 = c.Histo1D_div_Histo1D(deref(ao1._Histo1D()), deref(ao2._Histo1D()))
#         return cutil.new_owned_cls(Scatter2D, s2.newclone())
#     elif type(ao1) is Profile1D:
#         s2 = c.Profile1D_div_Profile1D(deref(ao1._Profile1D()), deref(ao2._Profile1D()))
#         return cutil.new_owned_cls(Scatter2D, s2.newclone())
#     raise ValueError("TODO: Only division of Histo1D and Profile1D supported so far... please contact the developers!")


def linspace(nbins, xmin, xmax):
    """(int, float, float) -> list[float]
    Make a list of n+1 bin edges linearly spaced between xmin and xmax, with the first and
    last edges on those boundaries."""
    return c.linspace(nbins, xmin, xmax)


def logspace(nbins, xmin, xmax):
    """(int, float, float) -> list[float]
    Make a list of n+1 bin edges linearly spaced on the interval log(xmin..xmax), with
    the first and last edges on those boundaries."""
    return c.logspace(nbins, xmin, xmax)


def pdfspace(nbins, xmin, xmax, fn, nsample=10000):
    """(int, float, float, [int]) -> list[float]
    Make a list of n+1 bin edges spaced with density proportional to fn(x) between
    xmin and xmax, with the first and last edges on those boundaries.

    The density is manually integrated by the Trapezium Rule, using nsample linspace points.

    Note: manually implemented in Python here rather than mapping the C++ version, since that
    requires some awkward Cython shim work:
    https://stackoverflow.com/questions/39044063/pass-a-closure-from-cython-to-c
    https://github.com/cython/cython/tree/master/Demos/callback
    """
    dx = (xmax-xmin)/float(nsample)
    xs = linspace(xmin, xmax, nsample+1)
    ys = [max(fn(x), 0) for x in xs]
    areas = [(ys[i] + ys[i+1])*dx/2. for i in range(nsample)]
    #areas = (ys[:-1] + ys[1:])*dx/2
    da = sum(areas)/nbins
    asum = 0
    xedges = [xmin]
    for i in range(nsample-1):
        asum += areas[i]
        if asum > da:
            asum = 0
            xedges.append(xs[i+1])
    xedges.append(xmax)
    assert(len(xedges) == nbins+1)
    return xedges


def index_between(x, binedges):
    """(float, list[float]) -> int
    Return the index of the bin which would contain x, or -1 if there is no enclosing
    bin in the given set of n+1 bin edges."""
    return c.index_between(x, binedges)


def mean(sample):
    """(list[float]) -> float
    Return the unweighted mean of the entries in the provided sample list."""
    return c.mean(sample)


def covariance(sample1, sample2):
    """(list[float], list[float]) -> float
    Return the unweighted covariance of the two provided sample lists."""
    return c.covariance(sample1, sample2)


def correlation(sample1, sample2):
    """(float, list[float]) -> int
    Return the unweighted correlation of the two provided sample lists."""
    return c.correlation(sample1, sample2)
