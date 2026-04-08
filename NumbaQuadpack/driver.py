import ctypes as ct
import importlib.util
from numba import njit, types
import numpy as np

quadpack_sig = types.double(types.double,
                            types.CPointer(types.double))

# Locate the compiled extension via Python's import machinery.
# This works on all platforms regardless of the platform-specific suffix
# (.cpython-310-x86_64-linux-gnu.so, .cp310-win_amd64.pyd, etc.)
_spec = importlib.util.find_spec('NumbaQuadpack.libcquadpack')
if _spec is None or _spec.origin is None:
    raise ImportError(
        "NumbaQuadpack.libcquadpack extension not found. "
        "Please reinstall ulysses: pip install ulysses"
    )

libquadpack = ct.CDLL(_spec.origin)

dqags_ = libquadpack.dqags
dqags_.argtypes = [ct.c_void_p, ct.c_double, ct.c_double, ct.c_double, ct.c_double, \
                 ct.c_void_p, ct.c_void_p, ct.c_void_p, ct.c_void_p,]
dqags_.restype = ct.c_double

@njit
def dqags(funcptr, a, b, data = np.array([0.0], np.float64), epsabs = 1.49e-08, epsrel = 1.49e-08):
    abserr = np.array(0.0,np.float64)
    neval = np.array(0,np.int32)
    ier = np.array(0,np.int32)
    
    sol = dqags_(funcptr, a, b, epsabs, epsrel, \
                 abserr.ctypes.data, neval.ctypes.data, \
                 ier.ctypes.data, data.ctypes.data)
    
    success = True
    if ier != 0:
        success = False
        
    return sol, abserr.item(), success
