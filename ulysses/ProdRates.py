import pandas as pd
import os
import warnings
from functools import lru_cache
import numpy as np
from scipy.integrate import quad
from mpmath import polylog
from scipy import interpolate
from scipy.special import expit
from scipy.special import kn
from scipy.special import zeta
from numba import njit
import numba as nb
import math
from NumbaQuadpack import quadpack_sig, dqags
from concurrent.futures import ThreadPoolExecutor

# Get the directory of the current Python file
current_dir = os.path.dirname(os.path.abspath(__file__))

# Define the path to the data/ directory
data_dir = os.path.join(current_dir, 'data')

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx, array[idx]

"""
Defining some constants
"""

g1  = 0.35
g2  = 0.65
ht  = 0.993
lam = 0.129
gw  = 2
Tew = 160 #GeV
v   = 174 #GeV

def vev(T):
    return np.sqrt((1 - T**2/Tew**2)*np.heaviside(Tew-T, 0.5)*v**2)

def mH(T):
    return np.sqrt(2*lam*vev(T)**2)

def mW(T):
    return vev(T)*g2/np.sqrt(2)

def mZ(T):
    return vev(T)*np.sqrt(g1**2 + g2**2)/np.sqrt(2)

def m_ell(T):
    return (T/4) * np.sqrt(g1**2 + 3*g2**2)

def m_phi(T):
    return (T/4) * np.sqrt(g1**2 + 3*g2**2 + 4*ht**2 + 8*lam)

#(m_ell/M)**2
def xell(z): 
    return (1/(16*z**2)) * (g1**2 + 3*g2**2) 
    
#(m_phi/M)**2
def xphi(z):
    return (1/(16*z**2)) * (g1**2 + 3*g2**2 + 4*ht**2 + 8*lam)

def lambda_kin(x, y, z):
    return (x**2 + y**2 + z**2 - 2*x*y - 2*x*z - 2*y*z)

def f_F(x):
    return expit(-x) #= 1/(1 + np.exp(x)) avoiding overflows

def f_F_2(x):
    return 1/(1 + np.exp(x))

def f_B(x):
    return 1/(np.exp(x) - 1)

def y0_onshell(y, z):
    if z != 0:
        return z * np.sqrt(1 + y**2 /z**2)
    else:
        return y
    
@njit
def _log_abs_ratio_nb(num, den):
    eps = 1e-14
    return math.log(math.fabs(num) + eps) - math.log(math.fabs(den) + eps)

@njit
def _f_F_nb(x):
    return 1.0 / (1.0 + math.exp(x))

@njit
def _f_B_nb(x):
    if x > 40.0:
        return math.exp(-x)
    if x < 1e-10:
        return 1e10
    return 1.0 / (math.exp(x) - 1.0)


"""
HELICITY PROJECTORS
"""


def yp(y, z):
    y0 = y0_onshell(y, z)
    return y0 + y

def ym(y, z):
    y0 = y0_onshell(y, z)
    return y0 - y
    

"""
Naive Contribution -- Symmetric Phase (s) -- Neutrino Decay part (N)
"""
#Neglecting the lepton mass in the 1to2 contribution does not lead to differences
#in the final result, but makes the plot smoother
#to include the lepton mass just replace xl with xell(z) (or xl = xell(z))
xl = 0 

def Ndecay_condition(z):
    return np.heaviside(1-np.sqrt(xl) - np.sqrt(xphi(z)), 0.5)

def xp_sN(y, z):
    y0 = y0_onshell(y, z) #On-shell condition
    return (1/2) * (y0 * (1 + xl - xphi(z)) + y * np.sqrt(lambda_kin(1,xl ,xphi(z))))

def xm_sN(y, z):
    y0 = y0_onshell(y, z) #On-shell condition
    return (1/2) * (y0 * (1 + xl  - xphi(z)) - y * np.sqrt(lambda_kin(1,xl ,xphi(z))))

#Integrals from PhysRevD.87.085009
def In_num(n, y, z):
    y0 = y0_onshell(y, z) #On-shell condition
    def integrand(x):
        return x**n * (1 - f_F(x) + f_B(y0-x))
    return quad(integrand, xm_sN(y, z), xp_sN(y, z))[0]

def In(n, y, z):
    y0 = y0_onshell(y, z)
    ym = xm_sN(y,z)
    yp = xp_sN(y,z)
    if n == 0:
        l1 = np.log(1 - np.exp(ym-y0))
        l2 = np.log(1 - np.exp(yp-y0))
        l3 = np.logaddexp(0, yp) - np.logaddexp(0, ym)
        return l1 - l2 + l3
    if n == 1:
        l1 = yp * (np.logaddexp(0, yp) - np.log(1 - np.exp(yp - y0)))
        l2 = -ym * (np.logaddexp(0, ym) - np.log(1 - np.exp(ym - y0)))
        l3 = polylog(2, np.exp(ym - y0)) - polylog(2, np.exp(yp - y0))
        l4 = polylog(2, -np.exp(yp)) - polylog(2, -np.exp(ym))
        return l1 + l2 + l3 + l4

#Sigma_0 / T
def Sigma_sN0(y, z):
    if Ndecay_condition(z) != 0:    
        return (gw/(16*np.pi*y)) * In(1, y, z)
    else: 
        return 0

#\sum_i (Sigma_i * k_i /(|k| T))
def Sigma_sNvec(y, z):
    y0 = y0_onshell(y, z) #On-shell condition
    if Ndecay_condition(z) != 0:
        return (gw/(16*np.pi * y**2)) * (y0 * In(1, y, z) - (z**2 / 2) * (1 + xl - xphi(z)) * In(0, y, z))
    else:
        return 0

#\Sigma_- / T
def Sigma_sNm(y, z):
    return Sigma_sN0(y, z) - Sigma_sNvec(y, z)

#\Sigma_+ / T
def Sigma_sNp(y, z):
    return Sigma_sN0(y, z) + Sigma_sNvec(y, z)

#\gamma_- / T
def gamma_sNm(y, z):
    y0 = y0_onshell(y, z)
    return Sigma_sNp(y, z) * ym(y, z) / y0

#\gamma_- / T
def gamma_sNp(y, z):
    y0 = y0_onshell(y, z)
    return Sigma_sNm(y, z) * yp(y, z) / y0

def ImPi_sN(y, z):
    f1 = Sigma_sNp(y, z) * ym(y, z)
    f2 = Sigma_sNm(y, z) * yp(y, z)
    return f1 + f2
    

"""
Naive Contribution -- Symmetric Phase (s) -- Higgs Decay part (H)
"""
def Hdecay_condition(z):
    return np.heaviside(np.sqrt(xphi(z)) - np.sqrt(xl) - 1, 0.5)

def xp_sH(y, z):
    y0 = y0_onshell(y, z) #On-shell condition
    return (1/2) * (y0 * (xphi(z) - xl - 1) + y * np.sqrt(lambda_kin(1,xl,xphi(z))))

def xm_sH(y, z):
    y0 = y0_onshell(y, z) #On-shell condition
    return (1/2) * (y0 * (xphi(z) - xl - 1) - y * np.sqrt(lambda_kin(1,xl,xphi(z))))


def Jn_num(n, y, z):
    y0 = y0_onshell(y, z) #On-shell condition
    def integrand(x):
        return x**n * (f_F(x) + f_B(y0+x))
    return quad(integrand, xm_sH(y, z), xp_sH(y, z))[0]

def Jn(n, y, z):
    y0 = y0_onshell(y, z)
    ym = xm_sH(y,z)
    yp = xp_sH(y,z)
    if n == 0:

            # Fermionic part: log(1 + e^{-x})
            Fm = np.log1p(np.exp(-ym))
            Fp = np.log1p(np.exp(-yp))

            # Bosonic part: log(1 - e^{-(x+y0)})
            Bm = np.log(-np.expm1(-(ym + y0)))
            Bp = np.log(-np.expm1(-(yp + y0)))

            fermion = Fm - Fp
            boson   = Bp - Bm

            return fermion + boson
    if n == 1:
        l1 = -yp * (np.logaddexp(0, -yp) - np.log(1 - np.exp(-yp - y0)))
        l2 = ym * (np.logaddexp(0, -ym) - np.log(1 - np.exp(-ym - y0)))
        l3 = polylog(2, np.exp(-ym - y0)) - polylog(2, np.exp(-yp - y0))
        l4 = polylog(2, -np.exp(-yp)) - polylog(2, -np.exp(-ym))
        return l1 + l2 + l3 + l4
    
#Sigma_0 / T
def Sigma_sH0(y, z):
    if Hdecay_condition(z) != 0:
        return (gw/(16*np.pi*y)) * Jn(1, y, z)
    else:
        return 0

#\sum_i (Sigma_i * k_i /(|k| T))
def Sigma_sHvec(y, z):
    y0 = y0_onshell(y, z) #On-shell condition
    if Hdecay_condition(z) != 0:
        return (gw/(16*np.pi * y**2)) * (y0 * Jn(1, y, z) - (z**2 / 2) * (xphi(z) - xl - 1) * Jn(0, y, z))
    else:
        return 0
    
#Sigma_- / T
def Sigma_sHm(y, z):
    return Sigma_sH0(y, z) - Sigma_sHvec(y, z)

#Sigma_+ / T
def Sigma_sHp(y, z):
    return Sigma_sH0(y, z) + Sigma_sHvec(y, z)

#\gamma_- / T
def gamma_sHm(y, z):
    y0 = y0_onshell(y, z)
    return Sigma_sHp(y, z) * ym(y, z) / y0

#\gamma_- / T
def gamma_sHp(y, z):
    y0 = y0_onshell(y, z)
    return Sigma_sHm(y, z) * yp(y, z) / y0

#\Im\Pi/T^2
def ImPi_sH(y, z):
    f1 = Sigma_sHp(y, z) * ym(y, z)
    f2 = Sigma_sHm(y, z) * yp(y, z)
    return f1 + f2


"""
Naive contribution -- Broken Phase (b) -- N to Higgs Decay Contribution (H) 
"""

def xL(M, T):
    return 0#m_ell(T)/M

def xH(M, T):
    return (mH(T)/M)**2

def Ndecay_condition_H(M, T):
    return np.heaviside(1-np.sqrt(xL(M, T)) - np.sqrt(xH(M, T)), 0.5)

def xp_bH(y, M, T):
    y0 = y0_onshell(y, M/T) #On-shell condition
    return (1/2) * (y0 * (1 + xL(M, T) - xH(M, T)) + y * np.sqrt(lambda_kin(1,xL(M, T) ,xH(M, T))))

def xm_bH(y, M, T):
    y0 = y0_onshell(y, M/T) #On-shell condition
    return (1/2) * (y0 * (1 + xL(M, T)  - xH(M, T)) - y * np.sqrt(lambda_kin(1,xL(M, T) ,xH(M, T))))

#Integrals from PhysRevD.87.085009
def In_bH_num(n, y, M, T):
    y0 = y0_onshell(y, M/T) #On-shell condition
    def integrand(x):
        return x**n * (1 - f_F(x) + f_B(y0-x))
    return quad(integrand, xm_bH(y, M, T), xp_bH(y, M, T))[0]

def In_bH(n, y, M, T):
    z = M/T
    y0 = y0_onshell(y, z)
    ym = xm_bH(y,M, T)
    yp = xp_bH(y,M, T)
    if n == 0:
        l1 = np.log(1 - np.exp(ym-y0))
        l2 = np.log(1 - np.exp(yp-y0))
        l3 = np.logaddexp(0, yp) - np.logaddexp(0, ym)
        return l1 - l2 + l3
    if n == 1:
        l1 = yp * (np.logaddexp(0, yp) - np.log(1 - np.exp(yp - y0)))
        l2 = -ym * (np.logaddexp(0, ym) - np.log(1 - np.exp(ym - y0)))
        l3 = polylog(2, np.exp(ym - y0)) - polylog(2, np.exp(yp - y0))
        l4 = polylog(2, -np.exp(yp)) - polylog(2, -np.exp(ym))
        return l1 + l2 + l3 + l4

#Sigma_0 / T
def Sigma_bH0(y, M, T):
    if Ndecay_condition_H(M, T) != 0:    
        return (1/(16*np.pi*y)) * In_bH(1, y, M, T)
    else: 
        return 0

#\sum_i (Sigma_i * k_i /(|k| T))
def Sigma_bHvec(y, M, T):
    y0 = y0_onshell(y, M/T) #On-shell condition
    if Ndecay_condition_H(M, T) != 0:
        return (1/(16*np.pi * y**2)) * (y0 * In_bH(1, y, M, T) - ((M**2 /(2*T**2)) * (1 + xL(M, T) - xH(M, T)) * In_bH(0, y, M, T)))
    else:
        return 0

#\Sigma_- / T
def Sigma_bHm(y, M, T):
    return (1/2) * (Sigma_bH0(y, M, T) - Sigma_bHvec(y, M, T))

#\Sigma_+ / T
def Sigma_bHp(y, M, T):
    return (1/2) * (Sigma_bH0(y, M, T) + Sigma_bHvec(y, M, T))

#\gamma_- / T
def gamma_bHm(y, M, T):
    y0 = y0_onshell(y, M/T)
    return Sigma_bHp(y, M, T) * ym(y, M/T) / y0

#\gamma_- / T
def gamma_bHp(y, M, T):
    y0 = y0_onshell(y, M/T)
    return Sigma_bHm(y, M, T) * yp(y, M/T) / y0

def ImPi_bH(y, M, T):
    f1 = Sigma_bHp(y, M, T) * ym(y, M/T)
    f2 = Sigma_bHm(y, M, T) * yp(y, M/T)
    return f1 + f2


"""
Naive contribution -- Broken Phase (b) -- N to W Decay Contribution (W) 
"""

def xW(M, T):
    return (mW(T)/M)**2

def Ndecay_condition_W(M, T):
    return np.heaviside(1-np.sqrt(xL(M, T)) - np.sqrt(xW(M, T)), 0.5)

def xp_bW(y, M, T):
    y0 = y0_onshell(y, M/T) #On-shell condition
    return (1/2) * (y0 * (1 + xL(M, T) - xW(M, T)) + y * np.sqrt(lambda_kin(1,xL(M, T) ,xW(M, T))))

def xm_bW(y, M, T):
    y0 = y0_onshell(y, M/T) #On-shell condition
    return (1/2) * (y0 * (1 + xL(M, T)  - xW(M, T)) - y * np.sqrt(lambda_kin(1,xL(M, T) ,xW(M, T))))

#Integrals from PhysRevD.87.085009
def In_bW_num(n, y, M, T):
    y0 = y0_onshell(y, M/T) #On-shell condition
    def integrand(x):
        return x**n * (1 - f_F(x) + f_B(y0-x))
    return quad(integrand, xm_bW(y, M, T), xp_bW(y, M, T))[0]

def In_bW(n, y, M, T):
    z = M/T
    y0 = y0_onshell(y, z)
    ym = xm_bW(y, M, T)
    yp = xp_bW(y, M, T)
    if n == 0:
        l1 = np.log(1 - np.exp(ym-y0))
        l2 = np.log(1 - np.exp(yp-y0))
        l3 = np.logaddexp(0, yp) - np.logaddexp(0, ym)
        return l1 - l2 + l3
    if n == 1:
        l1 = yp * (np.logaddexp(0, yp) - np.log(1 - np.exp(yp - y0)))
        l2 = -ym * (np.logaddexp(0, ym) - np.log(1 - np.exp(ym - y0)))
        l3 = polylog(2, np.exp(ym - y0)) - polylog(2, np.exp(yp - y0))
        l4 = polylog(2, -np.exp(yp)) - polylog(2, -np.exp(ym))
        return l1 + l2 + l3 + l4

#Sigma_0 / T
def Sigma_bW0(y, M, T):
    if Ndecay_condition_W(M, T) != 0:    
        return (1/(16*np.pi*y)) * In_bW(1, y, M, T)
    else: 
        return 0

#\sum_i (Sigma_i * k_i /(|k| T))
def Sigma_bWvec(y, M, T):
    y0 = y0_onshell(y, M/T) #On-shell condition
    if Ndecay_condition_W(M, T) != 0:
        return (1/(16*np.pi * y**2)) * (y0 * In_bW(1, y, M, T) - ((M**2 /(2*T**2)) * (1 + xL(M, T) - xW(M, T)) * In_bW(0, y, M, T)))
    else:
        return 0

#\Sigma_- / T
def Sigma_bWm(y, M, T):
    return (1/2) * (Sigma_bW0(y, M, T) - Sigma_bWvec(y, M, T))

#\Sigma_+ / T
def Sigma_bWp(y, M, T):
    return (1/2) * (Sigma_bW0(y, M, T) + Sigma_bWvec(y, M, T))

#\gamma_- / T
def gamma_bWm(y, M, T):
    y0 = y0_onshell(y, M/T)
    return Sigma_bWp(y, M, T) * ym(y, M/T) / y0

#\gamma_- / T
def gamma_bWp(y, M, T):
    y0 = y0_onshell(y, M/T)
    return Sigma_bWm(y, M, T) * yp(y, M/T) / y0

def ImPi_bW(y, M, T):
    f1 = Sigma_bWp(y, M, T) * ym(y, M/T)
    f2 = Sigma_bWm(y, M, T) * yp(y, M/T)
    return f1 + f2

"""
Naive contribution -- Broken Phase (b) -- N to Z Decay Contribution (Z) 
"""

def xZ(M, T):
    return (mZ(T)/M)**2

def Ndecay_condition_Z(M, T):
    return np.heaviside(1-np.sqrt(xL(M, T)) - np.sqrt(xZ(M, T)), 0.5)

def xp_bZ(y, M, T):
    y0 = y0_onshell(y, M/T) #On-shell condition
    return (1/2) * (y0 * (1 + xL(M, T) - xZ(M, T)) + y * np.sqrt(lambda_kin(1,xL(M, T) ,xZ(M, T))))

def xm_bZ(y, M, T):
    y0 = y0_onshell(y, M/T) #On-shell condition
    return (1/2) * (y0 * (1 + xL(M, T)  - xZ(M, T)) - y * np.sqrt(lambda_kin(1,xL(M, T) ,xZ(M, T))))

#Integrals from PhysRevD.87.085009
def In_bZ_num(n, y, M, T):
    y0 = y0_onshell(y, M/T) #On-shell condition
    def integrand(x):
        return x**n * (1 - f_F(x) + f_B(y0-x))
    return quad(integrand, xm_bZ(y, M, T), xp_bZ(y, M, T))[0]

def In_bZ(n, y, M, T):
    z = M/T
    y0 = y0_onshell(y, z)
    ym = xm_bZ(y,M, T)
    yp = xp_bZ(y,M, T)
    if n == 0:
        l1 = np.log(1 - np.exp(ym-y0))
        l2 = np.log(1 - np.exp(yp-y0))
        l3 = np.logaddexp(0, yp) - np.logaddexp(0, ym)
        return l1 - l2 + l3
    if n == 1:
        l1 = yp * (np.logaddexp(0, yp) - np.log(1 - np.exp(yp - y0)))
        l2 = -ym * (np.logaddexp(0, ym) - np.log(1 - np.exp(ym - y0)))
        l3 = polylog(2, np.exp(ym - y0)) - polylog(2, np.exp(yp - y0))
        l4 = polylog(2, -np.exp(yp)) - polylog(2, -np.exp(ym))
        return l1 + l2 + l3 + l4

#Sigma_0 / T
def Sigma_bZ0(y, M, T):
    if Ndecay_condition_Z(M, T) != 0:    
        return (1/(16*np.pi*y)) * In_bZ(1, y, M, T)
    else: 
        return 0

#\sum_i (Sigma_i * k_i /(|k| T))
def Sigma_bZvec(y, M, T):
    y0 = y0_onshell(y, M/T) #On-shell condition
    if Ndecay_condition_Z(M, T) != 0:
        return (1/(16*np.pi * y**2)) * (y0 * In_bZ(1, y, M, T) - ((M**2 /(2*T**2)) * (1 + xL(M, T) - xZ(M, T)) * In_bZ(0, y, M, T)))
    else:
        return 0

#\Sigma_- / T
def Sigma_bZm(y, M, T):
    return (1/2) * (Sigma_bZ0(y, M, T) - Sigma_bZvec(y, M, T))

#\Sigma_+ / T
def Sigma_bZp(y, M, T):
    return (1/2) * (Sigma_bZ0(y, M, T) + Sigma_bZvec(y, M, T))

#\gamma_- / T
def gamma_bZm(y, M, T):
    y0 = y0_onshell(y, M/T)
    return Sigma_bZp(y, M, T) * ym(y, M/T) / y0

#\gamma_- / T
def gamma_bZp(y, M, T):
    y0 = y0_onshell(y, M/T)
    return Sigma_bZm(y, M, T) * yp(y, M/T) / y0

def ImPi_bZ(y, M, T):
    f1 = Sigma_bZp(y, M, T) * ym(y, M/T)
    f2 = Sigma_bZm(y, M, T) * yp(y, M/T)
    return f1 + f2

"""
Broken Phase (b) -- Full Decay Contribution (D)
"""

def Sigma_bDm(y, M, T):
    return (Sigma_bHm(y, M, T) + Sigma_bZm(y, M, T) + 2 * (Sigma_bWm(y, M, T) ))

def Sigma_bDp(y, M, T):
    return (Sigma_bHp(y, M, T) + Sigma_bZp(y, M, T)+ 2 * (Sigma_bWp(y, M, T)  ))

def gamma_bDm(y, M, T):
    return (gamma_bHm(y, M, T) + gamma_bZm(y, M, T) + 2 * (gamma_bWm(y, M, T) ))

def gamma_bDp(y, M, T):
    return (gamma_bHp(y, M, T) + gamma_bZp(y, M, T) + 2 * (gamma_bWp(y, M, T) ))

"""
Broken Phase -- Indirect Contribution
"""

def V(m, y, T):
    A = m**2/(y*T)
    B = np.sqrt(m**2/T**2 + A - y)
    def Integrand1(x):
        return f_F(x) * (x * T + A * np.log(np.abs((2*A - x)/(2*A + x))))
    
    def Integrand2(x):
        return f_B(x) * (T * np.sqrt(x**2 - (m**2/T**2)) 
                         + A * np.log(np.abs((2*A - y 
                                              - np.sqrt(x**2 - (m**2/T**2)))/((2*A - y 
                                              + np.sqrt(x**2 - (m**2/T**2)))))))
    return (-T/(4*np.pi**2)) * (quad(Integrand1, 0, A)[0] + quad(Integrand1, A, 100)[0]
                                + quad(Integrand2, m/T, B)[0] + quad(Integrand2, B, 100)[0])

def SigmaH_m_full(y, M, T):
    y0 = y0_onshell(y, M/T)
    return (2 * g2**2 * V(mW(T), y, T) + (g1**2 + g2**2) * V(mZ(T), y, T))/(2*y0*T)

def SigmaA_m(y, M, T):
    return g2**2 * T /(2*np.pi)

#Valid for \pi*T >> mW, mZ, as in our case
def SigmaH_m(y, M, T):
    y0 = y0_onshell(y, M/T)
    return -m_ell(T)**2 /(2*y0*T)


#See 2103.16545, 1703.06085 and 1605.07720
def gamma_indirect_m_(y, M, T):
    y0 = y0_onshell(y, M/T)
    return (vev(T)**2/(T)) * (ym(y, M/T)/y0) * SigmaA_m(y, M, T) /(SigmaA_m(y, M, T)**2 + (T*ym(y, M/T) - SigmaH_m(y, M, T))**2)

def gamma_indirect_p(y, M, T):
    y0 = y0_onshell(y, M/T)
    bnu = -m_ell(T)**2 /(2*y0*T)
    gamma_nu = g2**2 * T /(2*np.pi)
    return (vev(T)**2/(T)) * (yp(y, M/T)/y0) * (gamma_nu/2) /((T*yp(y, M/T) - bnu)**2 + (gamma_nu/2)**2)

def gamma_indirect_m(y, M, T):
    y0 = y0_onshell(y, M/T)
    bnu = -m_ell(T)**2 /(2*y0*T)
    gamma_nu = g2**2 * T /(2*np.pi)
    return (vev(T)**2/(T)) * (ym(y, M/T)/y0) * (gamma_nu/2) /((T*ym(y, M/T) - bnu)**2 + (gamma_nu/2)**2)


@nb.cfunc(quadpack_sig)
def _A_integrand(w, data_):
    data = nb.carray(data_, (6,))
    y0 = data[0]; y = data[1]; D1 = data[2]; D2 = data[3]; z1 = data[4]; z2 = data[5]
    o1 = math.sqrt(w*w + z1*z1)
    o2 = math.sqrt(w*w + z2*z2)
    L1p = (_log_abs_ratio_nb(D1 - 2*o1*y0 - 2*y*w, D1 - 2*o1*y0 + 2*y*w)
         + _log_abs_ratio_nb(D1 + 2*o1*y0 - 2*y*w, D1 + 2*o1*y0 + 2*y*w))
    L1m = (_log_abs_ratio_nb(D1 - 2*o1*y0 - 2*y*w, D1 - 2*o1*y0 + 2*y*w)
         - _log_abs_ratio_nb(D1 + 2*o1*y0 - 2*y*w, D1 + 2*o1*y0 + 2*y*w))
    L2p = (_log_abs_ratio_nb(D2 + 2*o2*y0 + 2*y*w, D2 + 2*o2*y0 - 2*y*w)
         + _log_abs_ratio_nb(D2 - 2*o2*y0 + 2*y*w, D2 - 2*o2*y0 - 2*y*w))
    L2m = (_log_abs_ratio_nb(D2 + 2*o2*y0 + 2*y*w, D2 + 2*o2*y0 - 2*y*w)
         - _log_abs_ratio_nb(D2 - 2*o2*y0 + 2*y*w, D2 - 2*o2*y0 - 2*y*w))
    term_B = (-(D2 + 2.0*y*y) / (2.0*y) * w/o2 * L2p
              - y0*w/y * L2m
              + 4.0*w*w/o2) * _f_B_nb(o2)
    term_F = (D1 / (2.0*y) * w/o1 * L1p
              - y0*w/y * L1m
              + 4.0*w*w/o1) * _f_F_nb(o1)
    return (term_B + term_F) / (8.0 * math.pi * math.pi)

@nb.cfunc(quadpack_sig)
def _B_integrand(w, data_):
    data = nb.carray(data_, (6,))
    y0 = data[0]; y = data[1]; D1 = data[2]; D2 = data[3]; z1 = data[4]; z2 = data[5]
    zN2 = y0*y0 - y*y
    o1 = math.sqrt(w*w + z1*z1)
    o2 = math.sqrt(w*w + z2*z2)
    L1p = (_log_abs_ratio_nb(D1 - 2*o1*y0 - 2*y*w, D1 - 2*o1*y0 + 2*y*w)
         + _log_abs_ratio_nb(D1 + 2*o1*y0 - 2*y*w, D1 + 2*o1*y0 + 2*y*w))
    L1m = (_log_abs_ratio_nb(D1 - 2*o1*y0 - 2*y*w, D1 - 2*o1*y0 + 2*y*w)
         - _log_abs_ratio_nb(D1 + 2*o1*y0 - 2*y*w, D1 + 2*o1*y0 + 2*y*w))
    L2p = (_log_abs_ratio_nb(D2 + 2*o2*y0 + 2*y*w, D2 + 2*o2*y0 - 2*y*w)
         + _log_abs_ratio_nb(D2 - 2*o2*y0 + 2*y*w, D2 - 2*o2*y0 - 2*y*w))
    L2m = (_log_abs_ratio_nb(D2 + 2*o2*y0 + 2*y*w, D2 + 2*o2*y0 - 2*y*w)
         - _log_abs_ratio_nb(D2 - 2*o2*y0 + 2*y*w, D2 - 2*o2*y0 - 2*y*w))
    term_B = (zN2/y * w * L2m
             + D2/2.0 * y0/y * w/o2 * L2p
             - 4.0*y0*w*w/o2) * _f_B_nb(o2)
    term_F = (zN2/y * w * L1m
             - D1/2.0 * y0/y * w/o1 * L1p
             - 4.0*y0*w*w/o1) * _f_F_nb(o1)
    return (term_B + term_F) / (8.0 * math.pi * math.pi)

# .address gives the raw integer pointer — dqags passes it to the underlying
# C function as void*, which Numba handles correctly for integer arguments.
_A_addr = _A_integrand.address
_B_addr = _B_integrand.address

@njit
def _dqags_A(y0, y, D1, D2, z1, z2):
    data = np.empty(6)
    data[0] = y0; data[1] = y; data[2] = D1
    data[3] = D2; data[4] = z1; data[5] = z2
    result, _, _ = dqags(_A_addr, 1e-6, 80.0, data=data)
    return result

@njit
def _dqags_B(y0, y, D1, D2, z1, z2):
    data = np.empty(6)
    data[0] = y0; data[1] = y; data[2] = D1
    data[3] = D2; data[4] = z1; data[5] = z2
    result, _, _ = dqags(_B_addr, 1e-6, 80.0, data=data)
    return result


def A(y, z, z1, z2):
    y0 = y0_onshell(y, z)
    delta = z2**2 - z1**2
    D1 = y0**2 - y**2 - delta
    D2 = y0**2 - y**2 + delta
    return _dqags_A(y0, y, D1, D2, float(z1), float(z2)) / y**2

def B(y, z, z1, z2):
    y0 = y0_onshell(y, z)
    delta = z2**2 - z1**2
    D1 = y0**2 - y**2 - delta
    D2 = y0**2 - y**2 + delta
    return _dqags_B(y0, y, D1, D2, float(z1), float(z2)) / y**2

def anu_approx(y, M, T):
    yell = np.sqrt(g1**2 + 3*g2**2)/4
    y0   = y0_onshell(y, M/T)
    return (yell**2/y**2) * (1 - y0/(2*y) * np.log((y0+y)/(y0-y))) 

def bnu_approx(y, M, T):
    yell = np.sqrt(g1**2 + 3*g2**2)/4
    y0   = y0_onshell(y, M/T)
    return (yell**2/y**2) * (-y0 + (y0**2-y**2)/(2*y) * np.log((y0+y)/(y0-y))) 

_quad_pool = ThreadPoolExecutor(max_workers=4)

@lru_cache(maxsize=32768)
def _anu_full_cached(y, M, T):
    z = M / T
    zZ, zW = mZ(T) / T, mW(T) / T
    fZ = _quad_pool.submit(A, y, z, 0.0, zZ)
    fW = _quad_pool.submit(A, y, z, 0.0, zW)
    return (g1**2 + g2**2) / 4 * fZ.result() + g2**2 / 4 * fW.result()

@lru_cache(maxsize=32768)
def _bnu_full_cached(y, M, T):
    z = M / T
    zZ, zW = mZ(T) / T, mW(T) / T
    fZ = _quad_pool.submit(B, y, z, 0.0, zZ)
    fW = _quad_pool.submit(B, y, z, 0.0, zW)
    return (g1**2 + g2**2) / 4 * fZ.result() + g2**2 / 4 * fW.result()


def anu_full(y, M, T):
    y = float(y)
    M = float(M)
    T = float(T)
    return _anu_full_cached(y, M, T)


def bnu_full(y, M, T):
    y = float(y)
    M = float(M)
    T = float(T)
    return _bnu_full_cached(y, M, T)

# ---------- Thermal helpers (dimensionless, normalized by T) ----------

def _safe_pos(x, eps=1e-14):
    return x if x > eps else eps

def _log_ratio_1plus4y2_over_m2(y, m2_num, m2_den, eps=1e-14):
    # computes log[(1 + 4 y^2/m2_num)/(1 + 4 y^2/m2_den)]
    a = 1.0 + 4.0 * y**2 / _safe_pos(m2_num, eps)
    b = 1.0 + 4.0 * y**2 / _safe_pos(m2_den, eps)
    return np.log(a / b)

def _mix_angles_and_masses_tilde_bar(zW2, zZ2, g1, g2, nS=1.0, nG=3.0):
    # Debye masses divided by T^2
    e1 = (nS/6.0 + 5.0*nG/9.0) * g1**2
    e2 = (2.0/3.0 + nS/6.0 + nG/3.0) * g2**2

    # vacuum weak mixing
    s2t = 2.0 * g1 * g2 / (g1**2 + g2**2)
    c2t = (g2**2 - g1**2) / (g1**2 + g2**2)

    # tilde sector
    disc_t = np.sqrt(s2t**2 * zZ2**2 + (c2t * zZ2 + e2 - e1)**2)
    s2tt = s2t * zZ2 / _safe_pos(disc_t)
    c2tt = (c2t * zZ2 + e2 - e1) / _safe_pos(disc_t)

    # cos^2(theta-tilde), sin^2(theta-tilde) from double-angle data only
    c_d_t = 0.5 * (1.0 + c2t * c2tt + s2t * s2tt)
    s_d_t = 1.0 - c_d_t

    zt_plus2  = 0.5 * (zZ2 + e1 + e2 + disc_t)
    zt_minus2 = 0.5 * (zZ2 + e1 + e2 - disc_t)

    ztW2 = zW2 + e2
    ztZ2 = zt_plus2
    ztQ2 = zt_minus2

    # bar sector
    disc_b = np.sqrt(s2t**2 * zZ2**2 + (c2t * zZ2 + 0.5*(e2 - e1))**2)
    s2tb = s2t * zZ2 / _safe_pos(disc_b)
    c2tb = (c2t * zZ2 + 0.5*(e2 - e1)) / _safe_pos(disc_b)

    c_d_b = 0.5 * (1.0 + c2t * c2tb + s2t * s2tb)
    s_d_b = 1.0 - c_d_b

    zb_plus2  = 0.5 * (zZ2 + 0.5*(e1 + e2) + disc_b)
    zb_minus2 = 0.5 * (zZ2 + 0.5*(e1 + e2) - disc_b)

    zbW2 = zW2 + 0.5*e2
    zbZ2 = zb_plus2
    zbQ2 = zb_minus2

    return {
        "c_d_t": c_d_t, "s_d_t": s_d_t, "ztW2": ztW2, "ztZ2": ztZ2, "ztQ2": ztQ2,
        "c_d_b": c_d_b, "s_d_b": s_d_b, "zbW2": zbW2, "zbZ2": zbZ2, "zbQ2": zbQ2,
    }

# ---------- HTL pieces: hat means divided by T ----------

def Gamma_u_HTL_hat(y, T, nS=1.0, nG=3.0):
    zW2 = (mW(T)/T)**2
    zZ2 = (mZ(T)/T)**2
    th = _mix_angles_and_masses_tilde_bar(zW2, zZ2, g1, g2, nS=nS, nG=nG)

    su2 = 2.0 * g2**2 * _log_ratio_1plus4y2_over_m2(y, zW2, th["ztW2"])
    neutral = (g1**2 + g2**2) * (
        th["c_d_t"] * _log_ratio_1plus4y2_over_m2(y, zZ2, th["ztZ2"]) +
        th["s_d_t"] * _log_ratio_1plus4y2_over_m2(y, zZ2, th["ztQ2"])
    )
    return (su2 + neutral) / (16.0 * np.pi)


def Gamma_up2kK_HTL_hat(y, T, nS=1.0, nG=3.0):
    zW2 = (mW(T)/T)**2
    zZ2 = (mZ(T)/T)**2
    th = _mix_angles_and_masses_tilde_bar(zW2, zZ2, g1, g2, nS=nS, nG=nG)

    su2 = 2.0 * g2**2 * _log_ratio_1plus4y2_over_m2(y, zW2, th["zbW2"])
    neutral = (g1**2 + g2**2) * (
        th["c_d_b"] * _log_ratio_1plus4y2_over_m2(y, zZ2, th["zbZ2"]) +
        th["s_d_b"] * _log_ratio_1plus4y2_over_m2(y, zZ2, th["zbQ2"])
    )
    return (su2 + neutral) / (8.0 * np.pi)

# ---------- Born pieces (dimensionless chemical potentials mu/T) ----------

def l1b(x):
    # ln(1 - e^{-x})
    return np.log1p(-np.exp(-_safe_pos(x)))

def l1f(x):
    # ln(1 + e^{-x})
    return np.log1p(np.exp(-x))

def l2b(x):
    # Li_2(+e^{-x})
    return float(polylog(2, np.exp(-x)))

def l2f(x):
    # Li_2(-e^{-x})
    return float(polylog(2, -np.exp(-x)))


def Gamma_tilde_u_Born_hat(y, zM, mu1_hat=0.0, mu2_hat=0.0):
    # \tilde\Gamma_u / T
    a = zM**2 / (4.0*y)
    return (zM**2 / (32.0*np.pi*y**2)) * (
        l1f(a + mu1_hat) - l1b(y + a - mu2_hat)
    )


def Gamma_tilde_up2kK_Born_hat(y, zM, mu1_hat=0.0, mu2_hat=0.0):
    # (\tilde\Gamma_u + 2k\tilde\Gamma_K) / T
    a = zM**2 / (4.0*y)
    return (1.0 / (8.0*np.pi*y)) * (
        l2b(y + a - mu2_hat) - l2f(a + mu1_hat)
    )


def Gamma_u_Born_hat(y, T, mu_a_hat=0.0, mu_Q_hat=0.0):
    zZ = mZ(T)/T
    zW = mW(T)/T
    return ((g1**2 + g2**2) * Gamma_tilde_u_Born_hat(y, zZ, mu_a_hat, 0.0) +
            2.0*g2**2 * Gamma_tilde_u_Born_hat(y, zW, mu_a_hat - mu_Q_hat, mu_Q_hat))


def Gamma_up2kK_Born_hat(y, T, mu_a_hat=0.0, mu_Q_hat=0.0):
    zZ = mZ(T)/T
    zW = mW(T)/T
    return ((g1**2 + g2**2) * Gamma_tilde_up2kK_Born_hat(y, zZ, mu_a_hat, 0.0) +
            2.0*g2**2 * Gamma_tilde_up2kK_Born_hat(y, zW, mu_a_hat - mu_Q_hat, mu_Q_hat))

# ---------- Total Gu, Gk (dimensionless, i.e. divided by T) ----------

def Gu_Gk_hat(y, T, mu_a_hat=0.0, mu_Q_hat=0.0, nS=1.0, nG=3.0):
    Gu = Gamma_u_HTL_hat(y, T, nS=nS, nG=nG) + Gamma_u_Born_hat(y, T, mu_a_hat, mu_Q_hat)
    GupK = Gamma_up2kK_HTL_hat(y, T, nS=nS, nG=nG) + Gamma_up2kK_Born_hat(y, T, mu_a_hat, mu_Q_hat)

    # since (Gamma_u + 2k Gamma_K)/T = Gu + 2y Gk
    Gk = (GupK - Gu) / (2.0*y)
    return Gu, Gk, GupK

def gamma_indirect_m_full(y, M, T):
    Gu, Gk, GupK = Gu_Gk_hat(y, T)
    y0 = y0_onshell(y, M/T)
    anu = anu_full(y, M, T)
    bnu = bnu_full(y, M, T)
    #here v = 246 GeV, are we sure?
    return (vev(T)**2/(4*T**2)) * (ym(y, M/T)/y0) * (Gu + Gk * ym(y, M/T)) /((bnu + (1+anu) * ym(y, M/T))**2 + (Gu +Gk * ym(y, M/T))**2/4)

def gamma_indirect_m_rel(y, M, T):
    Gu, Gk, GupK = Gu_Gk_hat(y, T)
    y0 = y0_onshell(y, M/T)
    bnu = bnu_full(y, M, T)
    #here v = 246 GeV, are we sure?
    return (vev(T)**2/(4*T**2)) * (ym(y, M/T)/y0) * (Gu) /((bnu )**2 + (Gu)**2/4)

def gamma_indirect_p_full(y, M, T):
    Gu, Gk, GupK = Gu_Gk_hat(y, T)
    y0 = y0_onshell(y, M/T)
    anu = anu_full(y, M, T)
    bnu = bnu_full(y, M, T)
    #here v = 246 GeV, are we sure?
    return (vev(T)**2/(4*T**2)) * (yp(y, M/T)/y0) * (Gu + Gk * yp(y, M/T)) /((bnu + (1+anu) * yp(y, M/T))**2 + (Gu +Gk * yp(y, M/T))**2/4)

def gamma_indirect_p_rel(y, M, T):
    Gu, Gk, GupK = Gu_Gk_hat(y, T)
    y0 = y0_onshell(y, M/T)
    anu = anu_full(y, M, T)
    bnu = bnu_full(y, M, T)
    #here v = 246 GeV, are we sure?
    return (vev(T)**2/(4*T**2)) * (yp(y, M/T)/y0) * (Gu + Gk * 2 * y) /((bnu + (1+anu) * 2*y)**2 + (Gu +Gk * 2*y)**2/4)


def Sigma_indirect_m(y, M, T):
    y0 = y0_onshell(y, M/T)
    bnu = -m_ell(T)**2 /(2*y0*T)
    gamma_nu = g2**2 * T /(2*np.pi)
    return (vev(T)**2/(T)) * (gamma_nu/2) /((T*yp(y, M/T) - bnu)**2 + (gamma_nu/2)**2)

def Sigma_indirect_p(y, M, T):
    return (vev(T)**2/(T)) * SigmaA_m(y, M, T) /(SigmaA_m(y, M, T)**2 + (T*ym(y, M/T) - SigmaH_m(y, M, T))**2)

def Sigma_indirect_p_full(y, M, T):
    Gu, Gk, GupK = Gu_Gk_hat(y, T)
    y0 = y0_onshell(y, M/T)
    anu = anu_full(y, M, T)
    bnu = bnu_full(y, M, T)
    #here v = 246 GeV, are we sure?
    return (vev(T)**2/(4*T**2)) * (Gu + Gk * ym(y, M/T)) /((bnu + (1+anu) * ym(y, M/T))**2 + (Gu +Gk * ym(y, M/T))**2/4)


def Sigma_indirect_m_full(y, M, T):
    Gu, Gk, GupK = Gu_Gk_hat(y, T)
    y0 = y0_onshell(y, M/T)
    anu = anu_full(y, M, T)
    bnu = bnu_full(y, M, T)
    #here v = 246 GeV, are we sure?
    return (vev(T)**2/(4*T**2)) * (Gu + Gk * yp(y, M/T)) /((bnu + (1+anu) * yp(y, M/T))**2 + (Gu +Gk * yp(y, M/T))**2/4)

"""
Relativistic contribution extrapolated from PhysRevD.104.055010
and tables in http://www.laine.itp.unibe.ch/leptogenesis/
"""

#Take the data file and construct the arrays of data
Axes = pd.read_csv(data_dir + r'/Axes.txt', sep = " ")
Axes.drop('Unnamed: 0', inplace=True, axis=1, errors='ignore')

Qplus = pd.read_csv(data_dir + r'/Qplus.txt', sep = " ")
Qplus.drop('Unnamed: 0', inplace=True, axis=1, errors='ignore')

Qminus = pd.read_csv(data_dir + r'/Qminus.txt', sep = " ")
Qminus.drop('Unnamed: 0', inplace=True, axis=1, errors='ignore')

Rplus = pd.read_csv(data_dir + r'/Rplus.txt', sep = " ")
Rplus.drop('Unnamed: 0', inplace=True, axis=1, errors='ignore')

Rminus = pd.read_csv(data_dir + r'/Rminus.txt', sep = " ")
Rminus.drop('Unnamed: 0', inplace=True, axis=1, errors='ignore')

ks = Axes['k/T'][:249].to_numpy()
Ts = np.array([])
xs = np.array([])
for i in range(300):
    Ts = np.append(Ts, Axes['T/MeV'].iloc[249*i]*1e3) #in GeV
    xs = np.append(xs, Axes['log(Tmax/T)'].iloc[249*i])

def Qminus_interp(y, T):
    iT = find_nearest(xs, np.log(1e7/T))[0]
    ik = find_nearest(ks, y)[0]
    index = iT * len(ks) + ik
    return Qminus['Q11'].iloc[index]

def Qplus_interp(y, T):
    iT = find_nearest(xs, np.log(1e7/T))[0]
    ik = find_nearest(ks, y)[0]
    index = iT * len(ks) + ik
    return Qplus['Q11'].iloc[index]

def Rminus_interp(y, T):
    iT = find_nearest(xs, np.log(1e7/T))[0]
    ik = find_nearest(ks, y)[0]
    index = iT * len(ks) + ik
    return Rminus['R11'].iloc[index]

def Rplus_interp(y, T):
    iT = find_nearest(xs, np.log(1e7/T))[0]
    ik = find_nearest(ks, y)[0]
    index = iT * len(ks) + ik
    return Rplus['R11'].iloc[index]

#\Sigma_- / T
def Sigma_p_rel(y, T):
    iT = find_nearest(xs, np.log(1e7/T))[0]
    ik = find_nearest(ks, y)[0]
    index = iT * len(ks) + ik
    y = Axes['k/T'].iloc[index]
    return Qminus['Q11'].iloc[index] * 2*y**2 / gw #2*y**2 / gw #-- Check a potential factor of 2
    #return 1.9684e-2

#\Sigma_+ / T
def Sigma_m_rel(y, T):
    iT = find_nearest(xs, np.log(1e7/T))[0]
    ik = find_nearest(ks, y)[0]
    index = iT * len(ks) + ik
    return Qplus['Q11'].iloc[index]/(2*gw)
    #return 1.4612e-3

#\Sigma2_- / T
def Sigma2_p_rel(y, T):
    iT = find_nearest(xs, np.log(1e7/T))[0]
    ik = find_nearest(ks, y)[0]
    index = iT * len(ks) + ik
    y = Axes['k/T'].iloc[index]
    return Rminus['R11'].iloc[index] * 2*y**2 / gw #2*y**2 / gw #-- Check a potential factor of 2
    #return 1.9684e-2

#\Sigma2_+ / T
def Sigma2_m_rel(y, T):
    iT = find_nearest(xs, np.log(1e7/T))[0]
    ik = find_nearest(ks, y)[0]
    index = iT * len(ks) + ik
    return Rplus['R11'].iloc[index]/(2*gw)
    #return 1.4612e-3

#\gamma_- / T
def gamma_m_rel(y, M, T):
    y0 = y0_onshell(y, M/T)
    return Sigma_p_rel(y, T) * ym(y, M/T) * gw / y0 

#\gamma_+ / T
def gamma_p_rel(y, M, T):
    y0 = y0_onshell(y, M/T)
    return Sigma_m_rel(y, T) * yp(y, M/T) * gw / y0

#\Sigma_- / T
def Sigma_p_rel_broken(y, T):
    iT = find_nearest(xs, np.log(1e7/T))[0]
    ik = find_nearest(ks, y)[0]
    index = iT * len(ks) + ik
    y = Axes['k/T'].iloc[index]
    return Qplus['Q11'].iloc[index] * 2*y**2 / gw #2*y**2 / gw #-- Check a potential factor of 2
    #return 1.9684e-2

#\gamma_- / T
def gamma_m_rel_broken(y, M, T):
    y0 = y0_onshell(y, M/T)
    return Sigma_p_rel_broken(y, T) * ym(y, M/T) * gw / y0 

#\gamma_- / T
def gamma2_m_rel(y, M, T):
    y0 = y0_onshell(y, M/T)
    return Sigma2_p_rel(y, T) * ym(y, M/T) * gw / y0

#\gamma_+ / T
def gamma2_p_rel(y, M, T):
    y0 = y0_onshell(y, M/T)
    return Sigma2_m_rel(y, T) * yp(y, M/T) * gw / y0

def ImPi_rel(y, M, T):
    f1 = Sigma_p_rel(y, T) * ym(y, M/T) * gw
    f2 = Sigma_m_rel(y, T) * yp(y, M/T) * gw
    return f1 + f2

"""
Full Contribution
"""

def Sigma_m(y, M, T):
    z = M/T
    if T >= Tew:
        if 1-np.sqrt(xphi(z)) >= 0:
            Symmetric = np.heaviside(1-np.sqrt(xphi(z)),0.5) * (Sigma_sNm(y, z) + Sigma_sHm(y, z))
            return Sigma_m_rel(y, T) + Symmetric
        else:
            return Sigma_m_rel(y,T)
    if T < Tew:
        if 1-mW(T)/M >= 0:
            Broken = Sigma_indirect_m_full(y, M, T) + np.heaviside(1-mW(T)/M,0.5) * Sigma_bDm(y, M, T)
            return Broken + Sigma_m_rel(y, T)
        else:
            return Sigma_indirect_m_full(y, M, T) + Sigma_m_rel(y, T)

def Sigma_p(y, M, T):
    z = M/T
    if T >= Tew:
        if 1-np.sqrt(xphi(z)) >= 0:
            Symmetric = np.heaviside(1-np.sqrt(xphi(z)),0.5) * (Sigma_sNp(y, z) + Sigma_sHp(y, z))
            return Sigma_p_rel(y, T) + Symmetric
        else:
            return Sigma_p_rel(y, T)
    if T < Tew:
        if 1-mW(T)/M>= 0:
            Broken = Sigma_indirect_p_full(y, M, T) + np.heaviside(1-mW(T)/M,0.5) * Sigma_bDp(y, M, T)
            return Broken# + Sigma_p_rel(y, T)
        else:
            return  Sigma_indirect_p_full(y, M, T)# + Sigma_p_rel(y, T)

def gamma_m(y, M, T):
    z = M/T
    if T >= Tew:
        if 1-np.sqrt(xphi(z))>= 0:
            Symmetric = np.heaviside(1-np.sqrt(xphi(z)),0.5) * (gamma_sNm(y, z) + gamma_sHm(y, z))
            return gamma_m_rel(y, M, T) + Symmetric
        else:
            return gamma_m_rel(y, M, T)
    if T < Tew:
        if 1-mW(T)/M>= 0:
            Broken = np.heaviside(1-mW(T)/M,0.5) * gamma_bDm(y, M, T) + gamma_indirect_m_full(y, M, T)
            return Broken# + gamma_m_rel(y, M, T)
        else:
            return gamma_indirect_m_full(y, M, T)# + gamma_m_rel(y, M, T)

def gamma_p(y, M, T):
    z = M/T
    if T >= Tew:
        if 1-np.sqrt(xphi(z))>=0:
            Symmetric = np.heaviside(1-np.sqrt(xphi(z)),0.5) * (gamma_sNp(y, z) + gamma_sHp(y, z))
            return gamma_p_rel(y, M, T) + Symmetric
        else:
            return gamma_p_rel(y, M, T)
    if T < Tew:
        if 1-mW(T)/M >= 0:
            Broken = np.heaviside(1-mW(T)/M,0.5) * gamma_bDp(y, M, T) + gamma_indirect_p_full(y, M, T)
            return Broken + gamma_p_rel(y, M, T)
        else:
            return gamma_indirect_p_full(y, M, T) + gamma_p_rel(y, M, T)

def ImPi(y, M, T):
    z = M/T
    return ImPi_rel(y,M, T) + np.heaviside(1-np.sqrt(xphi(z)),0.5) * (ImPi_sH(y,z) + ImPi_sN(y,z))

#Correction factor for using MB statistics in the relativistic regime
corr = 3*zeta(3)/4

def neq(z):
    def integrand(x):
        return x**2 * expit(-np.sqrt(z**2 + x**2))
    r = quad(integrand, 0, np.inf, epsabs=1e-13, epsrel=1e-13)[0]
    if r !=0:
        return r
    else:
        return z**2*kn(2,z)

def I_2(z):
    def integrand(x):
        return x**2 * expit(-np.sqrt(z**2 + x**2))
    r = quad(integrand, 0, np.inf, epsabs=1e-14, epsrel=1e-6)[0]
    if r !=0:
        return r
    else:
        return z**2*kn(2,z)

def averaged_gamma_p(M, T):
    I = 0
    z = M/T
    if z <=10:
        for n in range(len(ks)-1):
            delta_a = ks[n+1] - ks[n]
            a = (ks[n+1] + ks[n])/2
            I += a**2 * gamma_p(a, M, T)*expit(-np.sqrt(z**2 + a**2)) * delta_a
        return I/(I_2(z))
    else:
        return z/(16*np.pi)

def averaged_gamma_m(M, T):
    I = 0
    z = M/T
    if z <= 10:
        for n in range(len(ks)-1):
            delta_a = ks[n+1] - ks[n]
            a = (ks[n+1] + ks[n])/2
            I += a**2 * gamma_m(a, M , T)*expit(-np.sqrt(z**2 + a**2)) * delta_a
        return I/(I_2(z))
    else:
        return z/(16*np.pi)


def averaged_gamma1_p(M, T, b = 1):
    I = 0
    z = M/T
    if z <= 10:
        for n in range(len(ks)-1):
            delta_a = ks[n+1] - ks[n]
            a = (ks[n+1] + ks[n])/2
            gamma2 = gamma2_p_rel(a, M, T) * b
            I += a**2 * expit(-np.sqrt(z**2 + a**2)) * (gamma2 + gamma_p(a, M , T) * (1 - expit(-np.sqrt(z**2 + a**2))))  * delta_a
        return I/I_2(z) #This factor enters the equations with an additional factor NNeq/Nelleq ~ z^2K2(z)/2, as in any wash-out terms in the Boltzmann equations
    else:
        return z/(16*np.pi)
        
def averaged_gamma1_m(M, T, b = 1):
    I = 0
    z = M/T
    if z <= 10:
        for n in range(len(ks)-1):
                delta_a = ks[n+1] - ks[n]
                a  = (ks[n+1] + ks[n])/2
                S2 = gamma2_m_rel(a, M, T) * b 
                I += a**2 * expit(-np.sqrt(z**2 + a**2)) * (S2 + gamma_m(a, M , T)*(1 - expit(-np.sqrt(z**2 + a**2)))) * delta_a
        return I/I_2(z)#This factor enters the equations with an additional factor NNeq/Nelleq ~ z^2K2(z)/2, as in any wash-out terms in the Boltzmann equations
    else:
        return z/(16*np.pi)
    
    
def averaged_gamma2_p(M, T):
    I = 0
    z = M/T
    for n in range(len(ks)-1):
            delta_a = ks[n+1] - ks[n]
            a = (ks[n+1] + ks[n])/2
            I += a**2 * expit(-np.sqrt(z**2 + a**2)) * gamma2_p_rel(a, M, T) * delta_a
    return I/I_2(z) #This factor enters the equations with an additional factor NNeq/Nelleq ~ z^2K2(z)/2, as in any wash-out terms in the Boltzmann equations
    
def averaged_gamma2_m(M, T):
    I = 0
    z = M/T
    for n in range(len(ks)-1):
            delta_a = ks[n+1] - ks[n]
            a  = (ks[n+1] + ks[n])/2
            I += a**2 * expit(-np.sqrt(z**2 + a**2)) * gamma2_m_rel(a, M, T) * delta_a
    return I/I_2(z)#This factor enters the equations with an additional factor NNeq/Nelleq ~ z^2K2(z)/2, as in any wash-out terms in the Boltzmann equations



"""
HAMILTONIAN
"""

"""
Zero-temperature coefficient
"""

from scipy.optimize import curve_fit
from scipy.interpolate import interp1d


def av_one_over_y0(z):
    def integrand(x):
        return np.sqrt(x**2-1) * expit(-x*z)
    if z <= 100:
        return z**2*quad(integrand, 1, np.inf)[0]/(I_2(z))
    else: 
        return 8/(15+8*z)

def av_one_over_y(z):
    def integrand(x):
        return x * expit(-x*z)
    return z**2*quad(integrand, 1, np.inf)[0]/(I_2(z))

# Define the exponential function for fitting
def exp_func(x, a, b):
    return a * np.exp(b * x)

# Define the power-law function for fitting
def power_law_func(x, a, b):
    return a * x**(b)

# Define the power-law function for fitting
@njit
def linear_func(x, a, b):
    return a * x + b

@njit
def interp_linear(x, xgrid, ygrid):
    i = np.searchsorted(xgrid, x) - 1

    if i < 0:
        i = 0
    elif i >= len(xgrid) - 1:
        i = len(xgrid) - 2

    t = (x - xgrid[i]) / (xgrid[i+1] - xgrid[i])
    return ygrid[i] * (1 - t) + ygrid[i+1] * t


zeds = np.logspace(-3, 2, 100)
one_over_y = [av_one_over_y(z) for z in zeds]
one_over_y0 = np.array([av_one_over_y0(z) for z in zeds])
interp_one_over_y_ = interp1d(zeds, one_over_y, kind='cubic')
# interp_one_over_y0_ = interp_linear(zeds, one_over_y0, kind='cubic')

# Extract the last 10 points
zeds_last_10 = zeds[-10:]
one_over_y_last_10 = one_over_y[-10:]
one_over_y0_last_10 = one_over_y0[-10:]

# Logarithmic transformation
log_zeds_last_10 = np.log10(zeds_last_10)
log_one_over_y_last_10 = np.log10(one_over_y_last_10)
log_one_over_y0_last_10 = np.log10(one_over_y0_last_10)

# Fit the log with a linear function -- it corresponds to a power-law in the original space
popty, _ = curve_fit(linear_func, log_zeds_last_10, log_one_over_y_last_10)
popty0, _ = curve_fit(linear_func, log_zeds_last_10, log_one_over_y0_last_10)

#Interpolated function allow for faster evaluation
# Define the extended function
def interp_one_over_y(z):
    if z < 1e-3:
        return interp_one_over_y_(1e-3)
    elif 1e-3 <= z <= 100:
        return interp_one_over_y_(z)
    else:
        return 10**(linear_func(np.log10(z), *popty))
    
# Define the extended function
@njit
def interp_one_over_y0(z):
    if z < 1e-3:
        return interp_linear(1e-3, zeds, one_over_y0)
    elif z <= 100:
        return interp_linear(z, zeds, one_over_y0)
    else:
        return 8/(15+8*z)


"""
Contribution from Higgs'VEV evolution
"""

#see 1710.03744 for details on all these definitions
Tsph = 131.7#GeV

def hEV(M, T):
    return vev(T)**2/(T**2) * (Tsph/T) * av_one_over_y0(M/T)

"""
LNC and LNV contributions
"""

#Approximation given in 1710.03744, see Eq. C.48 - divided by 2 to match the term that enters the DMEs
#and an additional factor of 0.5 to match the rates used by Juraj
@njit
def hm_rel(z):
    return 0.25*(3.50-0.47*np.log(z**2)+3.47*np.log(z)*np.log(z)) * 1e-2 * z**2
    
def hm_rel_num(z):
    def integrand(y, z):
        y0 = y0_onshell(y, z)
        log_term = np.log(np.abs((y0 + y) / (y0 - y)))
        return (1 / y0) * (-1 + log_term) * expit(-y0)

    integral_result, _ = quad(integrand, 0, np.inf, args=(z,))
    return 2 / (3 * zeta(3)) * (z**2 / 16) * integral_result * 0.5

@njit
def hp_rel(z):
    return np.pi**2/(144*1.20205)

@njit
def h_non_rel(z):
    return .5/(15+8*z)

#Limited functions are then interpolated between 1e-3 and 1e2 for faster evaluation
def hp_limited(z):
    def integrand(t, z):
        # Change of variables to tan(t)
        y = np.tan(t)
        y0 = y0_onshell(y, z)
        log_term = 2 * np.arctanh(y / y0)
        return (y**2 / y0) * (-log_term * (z**2/y**2) + 2 * (1 + y0/y)) * expit(-y0) * (1 / np.cos(t)**2)
    integral_result, _ = quad(integrand, 0, np.pi/2, args=(z,), epsabs=1e-12, epsrel=1e-12)
    return (1 / I_2(z)) * (1 / 16) * integral_result *0.5 #We add an additional factor of 0.5 to match hp_rel used in DMEs of ULYSSESv2

def hm_limited(z):
    def integrand(t, z):
        y = np.tan(t)
        y0 = y0_onshell(y, z)
        log_term = 2 * np.arctanh(y / y0)
        return (y**2 / y0) * (log_term * (z**2/y**2) + 2 * (1 - y0/y)) * expit(-y0) * (1 / np.cos(t)**2)
    integral_result, _ = quad(integrand, 0, np.pi/2, args=(z,), epsabs=1e-13, epsrel=1e-13)
    return (1 / I_2(z)) * (1 / 16) * integral_result * 0.5 



zeds                =  np.logspace(-3, 2, 100)
hp_array            = np.array([hp_limited(z) for z in zeds])
hm_array            = np.array([hm_limited(z) for z in zeds])
interp_hp_          = interp1d(zeds, hp_array, kind='cubic')
interp_hm_          = interp1d(zeds, hm_array, kind='cubic')

# Extract the last 10 points
zeds_last_10 = zeds[-10:] 
hm_last_10    = hm_array[-10:]  
hp_last_10    = hp_array[-10:]  

# Logarithmic transformation
log_zeds_last_10 = np.log10(zeds_last_10)
log_hm_last_10      = np.log10(hm_last_10)
log_hp_last_10      = np.log10(hp_last_10)

# Fit the exponential function to the transformed points
popt_hm, _ = curve_fit(linear_func, log_zeds_last_10, log_hm_last_10)
popt_hp, _ = curve_fit(linear_func, log_zeds_last_10, log_hp_last_10)
a_hp = popt_hp[0]
b_hp = popt_hp[1]


@njit
def interp_hp(z):
    #zm = 100
    #zp = 1000
    #c = 0.01

    if z <= 1e-3:
        return hp_rel(1e-3)
    elif 1e-3 < z <= 30:
        return interp_linear(z,zeds,hp_array)
    elif z > 30:
        return 10**(a_hp * np.log10(z) + b_hp)
        
    """
    elif zm < z <= zp:
        # Interpolate between the result at z = 100 and h_non_rel(z)
        hp_end = interp_hp_(zm)
        hNR_start = h_non_rel(zp)
        A = (hp_end**(zp**c/zm**c) / hNR_start)**(zm**c/(zp**c-zm**c))
        b = (1/(zp**c-zm**c)) * np.log(hp_end / hNR_start)
        return A*np.exp(-b*z**c)
    elif z > zp:
        return h_non_rel(z)
    """
    
@njit
def interp_hm(z):
    zm = 50
    zp = 1000
    c = 0.01

    if z <= 1e-3:
        return hm_rel(z)
    elif 1e-3 < z <= zm:
        return interp_linear(z,zeds,hm_array)
    
    #elif z > 30:
    #    return 10**(linear_func(np.log10(z), *popt_hm))


    elif zm < z <= zp:
        hm_end = interp_linear(zm,zeds,hm_array)
        hNR_start = h_non_rel(zp)
        A = (hm_end**(zp**c/zm**c) / hNR_start)**(zm**c/(zp**c-zm**c))
        b = (1/(zp**c-zm**c)) * np.log(hm_end / hNR_start)
        return A*np.exp(-b*z**c)
    elif z > zp:
        return h_non_rel(z)

""" INDIRECT CONTRIBUTION IN THE BROKEN PHASE"""

# h^ind_+ = v^2/(2T^2) * [b_nu + (1+a_nu)(y0+y)] / {[b_nu + (1+a_nu)(y0+y)]^2 + (1/4)[Gu + Gk*(y0+y)]^2}
def h_indirect_p(y, M, T):
    Gu, Gk, _ = Gu_Gk_hat(y, T)
    anu = anu_full(y, M, T)
    bnu = bnu_full(y, M, T)
    yp_val = yp(y, M/T)
    num = bnu + (1 + anu) * yp_val
    den = num**2 + 0.25 * (Gu + Gk * yp_val)**2
    return (vev(T)**2 / (2 * T**2)) * num / den

# h^ind_- = v^2/(2T^2) * [b_nu + (1+a_nu)(y0-y)] / {[b_nu + (1+a_nu)(y0-y)]^2 + (1/4)[Gu + Gk*(y0-y)]^2}
def h_indirect_m(y, M, T):
    Gu, Gk, _ = Gu_Gk_hat(y, T)
    anu = anu_full(y, M, T)
    bnu = bnu_full(y, M, T)
    ym_val = ym(y, M/T)
    num = bnu + (1 + anu) * ym_val
    den = num**2 + 0.25 * (Gu + Gk * ym_val)**2
    return (vev(T)**2 / (2 * T**2)) * num / den


def _stable_occ_weight(y, z):
    if not np.isfinite(y) or not np.isfinite(z):
        return 0.0

    e = np.hypot(z, y)
    # For large z use a rescaled MB-like weight exp(-(E-z)) so normalization
    # and numerator stay in a safe dynamic range; the common exp(-z) cancels.
    if z > 10:
        # Use E-z = y^2/(E+z) to avoid catastrophic cancellation at large z.
        denom = e + z
        if denom <= 0 or not np.isfinite(denom):
            return 0.0
        delta = (y * y) / denom
        if not np.isfinite(delta):
            return 0.0
        return np.exp(-delta)
    return expit(-e)

def averaged_h_indirect_p(M, T):
    I = 0
    N = 0
    z = M/T
    for n in range(len(ks)-1):
        delta_a = ks[n+1] - ks[n]
        a = (ks[n+1] + ks[n])/2
        w = _stable_occ_weight(a, z)
        h = h_indirect_p(a, M, T)
        if not np.isfinite(w) or not np.isfinite(h):
            continue
        I += a**2 * h * w * delta_a
        N += a**2 * w * delta_a
    return I / max(N, 1e-300)

def averaged_h_indirect_m(M, T):
    I = 0
    N = 0
    z = M/T
    for n in range(len(ks)-1):
        delta_a = ks[n+1] - ks[n]
        a = (ks[n+1] + ks[n])/2
        w = _stable_occ_weight(a, z)
        h = h_indirect_m(a, M, T)
        if not np.isfinite(w) or not np.isfinite(h):
            continue
        I += a**2 * h * w * delta_a
        N += a**2 * w * delta_a
    return I / max(N, 1e-300)
