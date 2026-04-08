import ctypes as ct
import glob
import os
from numba import njit, types
import numpy as np

quadpack_sig = types.double(types.double,
                            types.CPointer(types.double))

# Look for the compiled extension in the same directory as this file.
# This works for all install modes:
#   - regular install: driver.py and the .so are both in site-packages/NumbaQuadpack/
#   - editable install: setuptools copies the .so into the source NumbaQuadpack/
#   - build_ext --inplace: .so lands next to driver.py in the source tree
_pkg_dir = os.path.dirname(os.path.abspath(__file__))
_candidates = [f for f in glob.glob(os.path.join(_pkg_dir, 'libcquadpack.*'))
               if not f.endswith('.c')]

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
