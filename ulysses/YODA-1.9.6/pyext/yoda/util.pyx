from collections import namedtuple
from operator import itemgetter


def as_bool(x):
    if type(x) is bool:
        return x
    s = str(x)
    if s.lower() in ("true", "yes", "on", "1", "1.0"):
        return True
    if s.lower() in ("false", "no", "off", "0", "0.0"):
        return False
    raise Exception("'{}' cannot be parsed as a boolean flag".format(s))


def _autotype(var, autobool=False):
    """Automatically convert strings to numerical types if possible."""
    if type(var) is not str:
        return var
    ## Convert via Python ast parser
    try:
        import ast
        var = ast.literal_eval(var)
    except:
        # TODO: print a warning?
        pass
    ## Try friendly string conversions to bool
    if autobool and type(var) is str:
        try:
            var = as_bool(var)
        except:
            pass
    ## Finally return
    return var


# def _autonp(var):
#     """Automatically return lists as numpy arrays if numpy is imported"""
#     if "numpy" in dir():
#         return numpy.array(var)
#     elif "np" in dir():
#         return np.array(var)
#     else:
#         return var


def _autostr(var, precision=8):
    """Automatically format numerical types as the right sort of string."""
    if type(var) is float:
        return ("% ." + str(precision) + "e") % var
    elif not isinstance(var, (list,tuple)):
        return str(var)
    else:
        return ",".join(_autostr(subval) for subval in var)




cdef class Base:
    pass


def try_loop(fs, *args, char *_msg='Invalid arguments', **kwargs):
    for f in fs:
        try:
            f(*args, **kwargs)
            return
        except (TypeError, AttributeError):
            pass
    raise TypeError(_msg)


XY = namedtuple('XY', ('x', 'y'))
XYZ = namedtuple('XYZ', ('x', 'y', 'z'))
EdgePair = namedtuple('EdgePair', ('low', 'high'))
ErrorPair = namedtuple('ErrorPair', ('minus', 'plus'))


## Utils for handling error conversions to/from std::pair
from libcpp.pair cimport pair

def read_edge_pair(pair[double, double] es):
    return EdgePair(es.first, es.second)

def read_error_pair(pair[double, double] es):
    return ErrorPair(es.first, es.second)

def read_symmetric(val):
    try:
        a, b = val
    except TypeError:
        a = b = val
    return pair[double, double](a, b)
