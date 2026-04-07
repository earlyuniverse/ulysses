import ctypes as ct
import numba as nb
from numba import njit, types
import numpy as np
import os
import platform 

quadpack_sig = types.double(types.double,
                            types.CPointer(types.double))

rootdir = os.path.dirname(os.path.realpath(__file__))+'/'

# Determine the correct library filename based on the OS
system = platform.system()
if system == "Windows":
    lib_name = "cquadpack.dll"
elif system == "Darwin": # macOS
    lib_name = "libcquadpack.dylib"
else: # Linux and others
    lib_name = "libcquadpack.so"

lib_path = os.path.join(rootdir, lib_name)

# Check if file exists to give a better error message
if not os.path.exists(lib_path):
    raise ImportError(f"Could not find {lib_name} at {lib_path}. "
                      "Please ensure the library was compiled correctly.")

libquadpack = ct.CDLL(lib_path)

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
