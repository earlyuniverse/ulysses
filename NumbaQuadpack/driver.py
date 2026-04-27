import ctypes as ct
import glob
import os
import sys
from numba import njit, types
import numpy as np

quadpack_sig = types.double(types.double,
                            types.CPointer(types.double))

# Search for the compiled extension. We first check the directory containing
# this file (covers editable installs and build_ext --inplace), then fall back
# to scanning all sys.path entries for a NumbaQuadpack/ sub-directory.
# The fallback handles the case where Python loads driver.py from the source
# tree (e.g. repo root in cwd) but the extension was installed into a venv's
# site-packages — as happens with 'uv pip install .' or running from the repo
# root after a regular 'pip install .'.
def _find_libcquadpack():
    import sysconfig
    native_suffix = sysconfig.get_config_var('EXT_SUFFIX') or ''

    search_dirs = [os.path.dirname(os.path.abspath(__file__))]
    search_dirs += [os.path.join(p, 'NumbaQuadpack') for p in sys.path if p]

    # Collect ALL candidates across every search directory before choosing.
    # This avoids returning a foreign-platform .so from the source tree when
    # the correctly compiled native extension lives further down in sys.path
    # (e.g. site-packages), which happens when the repo root is on sys.path
    # via an editable install .pth file or Jupyter's working-directory injection.
    all_hits = []
    seen = set()
    for d in search_dirs:
        d = os.path.abspath(d)
        if d in seen:
            continue
        seen.add(d)
        hits = [f for f in glob.glob(os.path.join(d, 'libcquadpack.*'))
                if not f.endswith('.c') and not f.endswith('.h')]
        all_hits.extend(hits)

    if not all_hits:
        return None
    # Prefer the file whose suffix exactly matches this Python interpreter.
    native = [f for f in all_hits if f.endswith(native_suffix)]
    return native[0] if native else all_hits[0]

_ext_path = _find_libcquadpack()
if _ext_path is None:
    raise ImportError(
        "NumbaQuadpack.libcquadpack extension not found. "
        "Please reinstall: pip install ulysses"
    )

libquadpack = ct.CDLL(_ext_path)

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
