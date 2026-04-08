import ctypes as ct
import glob
import os
from numba import njit, types
import numpy as np

quadpack_sig = types.double(types.double,
                            types.CPointer(types.double))

# Locate the compiled extension via Python's import machinery.
# This works on all platforms regardless of the platform-specific suffix
# (.cpython-310-x86_64-linux-gnu.so, .cp310-win_amd64.pyd, etc.)
# Search every entry on sys.path for the compiled extension.
# This finds it in site-packages even when the source directory is
# earlier on sys.path (e.g. running a notebook from the repo root).
import sys
_candidates = []
for _p in sys.path:
    _candidates = [f for f in glob.glob(os.path.join(_p, 'NumbaQuadpack', 'libcquadpack.*'))
                   if not f.endswith('.c')]
    if _candidates:
        break

if not _candidates:
    raise ImportError(
        "NumbaQuadpack.libcquadpack extension not found. "
        "Please reinstall: pip install ulysses"
    )

libquadpack = ct.CDLL(_candidates[0])

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
