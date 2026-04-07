import pandas as pd
import os
import numpy as np
from scipy.integrate import quad
from mpmath import polylog
from scipy import interpolate
from scipy.special import expit
from scipy.special import kn
from scipy.special import zeta
from numba import jit, njit


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
v   = 246# GeV


#vev as in PhysRevLett.117.091801
def vev(T):
    return np.sqrt((1 - T**2/Tew**2)*np.heaviside(Tew-T, 0.5)*v**2)

def mH(T):
    return np.sqrt(2*lam*vev(T)**2)

def mW(T):
    return vev(T)*g2/2

def mZ(T):
    return vev(T)*np.sqrt(g1**2 + g2**2)/2

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
        if (ym + y0) >= 500:
            l1 = -(ym + y0)#It avoids the overflow of the exponential
        elif (ym + y0) < 500:
            l1 = -np.log(np.exp(ym+y0)-1)
        if (yp + y0) >= 500:
            l2 = yp + y0#It avoids the overflow of the exponential
        elif (yp + y0) < 500:
            l2 = np.log(np.exp(yp+y0)-1)
        l3 = np.logaddexp(0, ym) - np.logaddexp(0, yp)
        return l1 + l2 + l3
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
def gamma_indirect_m(y, M, T):
    y0 = y0_onshell(y, M/T)
    return (vev(T)**2/(2*T)) * (ym(y, M/T)/y0) * SigmaA_m(y, M, T) /(SigmaA_m(y, M, T)**2 + (T*ym(y, M/T) - SigmaH_m(y, M, T))**2)

def gamma_indirect_p(y, M, T):
    y0 = y0_onshell(y, M/T)
    bnu = -m_ell(T)**2 /(2*y0*T)
    gamma_nu = g2**2 * T /(2*np.pi)
    return (vev(T)**2/(2*T)) * (yp(y, M/T)/y0) * (gamma_nu/2) /((T*yp(y, M/T) - bnu)**2 + (gamma_nu/2)**2)

def Sigma_indirect_m(y, M, T):
    y0 = y0_onshell(y, M/T)
    bnu = -m_ell(T)**2 /(2*y0*T)
    gamma_nu = g2**2 * T /(2*np.pi)
    return (vev(T)**2/(2*T)) * (gamma_nu/2) /((T*yp(y, M/T) - bnu)**2 + (gamma_nu/2)**2)

def Sigma_indirect_p(y, M, T):
    return (vev(T)**2/(2*T)) * SigmaA_m(y, M, T) /(SigmaA_m(y, M, T)**2 + (T*ym(y, M/T) - SigmaH_m(y, M, T))**2)

"""
Relativistic contribution extrapolated from PhysRevD.104.055010
and tables in http://www.laine.itp.unibe.ch/leptogenesis/
"""

#Take the data file and construct the arrays of data
Axes = pd.read_csv(data_dir + r'/Axes.txt', sep = " ")
Axes.drop('Unnamed: 0', inplace  = True, axis = 1)

Qplus = pd.read_csv(data_dir + r'/Qplus.txt', sep = " ")
Qplus.drop('Unnamed: 0', inplace  = True, axis = 1)

Qminus = pd.read_csv(data_dir + r'/Qminus.txt', sep = " ")
Qminus.drop('Unnamed: 0', inplace  = True, axis = 1)

Rplus = pd.read_csv(data_dir + r'/Rplus.txt', sep = " ")
Rplus.drop('Unnamed: 0', inplace  = True, axis = 1)

Rminus = pd.read_csv(data_dir + r'/Rminus.txt', sep = " ")
Rminus.drop('Unnamed: 0', inplace  = True, axis = 1)

ks = Axes['k/T'][:249].to_numpy()
Ts = np.array([])
xs = np.array([])
for i in range(300):
    Ts = np.append(Ts, Axes['T/MeV'][249*i]*1e3) #in GeV
    xs = np.append(xs, Axes['log(Tmax/T)'][249*i])

def Qminus_interp(y, T):
    iT = find_nearest(xs, np.log(1e7/T))[0]
    ik = find_nearest(ks, y)[0]
    index = iT * len(ks) + ik
    return Qminus['Q11'][index]

def Qplus_interp(y, T):
    iT = find_nearest(xs, np.log(1e7/T))[0]
    ik = find_nearest(ks, y)[0]
    index = iT * len(ks) + ik
    return Qplus['Q11'][index]

def Rminus_interp(y, T):
    iT = find_nearest(xs, np.log(1e7/T))[0]
    ik = find_nearest(ks, y)[0]
    index = iT * len(ks) + ik
    return Rminus['R11'][index]

def Rplus_interp(y, T):
    iT = find_nearest(xs, np.log(1e7/T))[0]
    ik = find_nearest(ks, y)[0]
    index = iT * len(ks) + ik
    return Rplus['R11'][index]

#\Sigma_- / T
def Sigma_p_rel(y, T):
    iT = find_nearest(xs, np.log(1e7/T))[0]
    ik = find_nearest(ks, y)[0]
    index = iT * len(ks) + ik
    y = Axes['k/T'][index]
    return Qminus['Q11'][index] * 2*y**2 / gw #2*y**2 / gw #-- Check a potential factor of 2
    #return 1.9684e-2

#\Sigma_+ / T
def Sigma_m_rel(y, T):
    iT = find_nearest(xs, np.log(1e7/T))[0]
    ik = find_nearest(ks, y)[0]
    index = iT * len(ks) + ik
    return Qplus['Q11'][index]/(2*gw)
    #return 1.4612e-3

#\Sigma2_- / T
def Sigma2_p_rel(y, T):
    iT = find_nearest(xs, np.log(1e7/T))[0]
    ik = find_nearest(ks, y)[0]
    index = iT * len(ks) + ik
    y = Axes['k/T'][index]
    return Rminus['R11'][index] * 2*y**2 / gw #2*y**2 / gw #-- Check a potential factor of 2
    #return 1.9684e-2

#\Sigma2_+ / T
def Sigma2_m_rel(y, T):
    iT = find_nearest(xs, np.log(1e7/T))[0]
    ik = find_nearest(ks, y)[0]
    index = iT * len(ks) + ik
    return Rplus['R11'][index]/(2*gw)
    #return 1.4612e-3

#\gamma_- / T
def gamma_m_rel(y, M, T):
    y0 = y0_onshell(y, M/T)
    return Sigma_p_rel(y, T) * ym(y, M/T) * gw / y0

#\gamma_+ / T
def gamma_p_rel(y, M, T):
    y0 = y0_onshell(y, M/T)
    return Sigma_m_rel(y, T) * yp(y, M/T) * gw / y0

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
    if T >= 160:
        if 1-np.sqrt(xphi(z)) >= 0:
            Symmetric = np.heaviside(1-np.sqrt(xphi(z)),0.5) * (Sigma_sNm(y, z) + Sigma_sHm(y, z))
            return Sigma_m_rel(y, T) + Symmetric
        else:
            return Sigma_m_rel(y,T)
    if T < 160:
        if 1-mW(T)/M >= 0:
            Broken = Sigma_indirect_m(y, M, T) + np.heaviside(1-mW(T)/M,0.5) * Sigma_bDm(y, M, T)
            return Sigma_m_rel(y, T) + Broken
        else:
            return Sigma_m_rel(y, T)

def Sigma_p(y, M, T):
    z = M/T
    if T >= 160:
        if 1-np.sqrt(xphi(z)) >= 0:
            Symmetric = np.heaviside(1-np.sqrt(xphi(z)),0.5) * (Sigma_sNp(y, z) + Sigma_sHp(y, z))
            return Sigma_p_rel(y, T) + Symmetric
    if T < 160:
        if 1-mW(T)/M>= 0:
            Broken = Sigma_indirect_p(y, M, T) + np.heaviside(1-mW(T)/M,0.5) * Sigma_bDp(y, M, T)
            return Sigma_p_rel(y, T) + Broken
        else:
            return Sigma_p_rel(y, T)
    
def gamma_m(y, M, T):
    z = M/T
    if T >= 160:
        if 1-np.sqrt(xphi(z))>= 0:
            Symmetric = np.heaviside(1-np.sqrt(xphi(z)),0.5) * (gamma_sNm(y, z) + gamma_sHm(y, z))
            return gamma_m_rel(y, M, T) + Symmetric
        else:
            return gamma_m_rel(y, M, T)
    if T < 160:
        if 1-mW(T)/M>= 0:
            Broken = gamma_indirect_m(y, M, T) + np.heaviside(1-mW(T)/M,0.5) * gamma_bDm(y, M, T)
            return gamma_m_rel(y, M, T) + Broken
        else:
            return gamma_m_rel(y, M, T)

def gamma_p(y, M, T):
    z = M/T
    if T >= 160:
        if 1-np.sqrt(xphi(z))>=0:
            Symmetric = np.heaviside(1-np.sqrt(xphi(z)),0.5) * (gamma_sNp(y, z) + gamma_sHp(y, z))
            return gamma_p_rel(y, M, T) + Symmetric
        else:
            return gamma_p_rel(y, M, T)
    if T < 160:
        if 1-mW(T)/M >= 0:
            Broken = gamma_indirect_p(y, M, T) + np.heaviside(1-mW(T)/M,0.5) * gamma_bDp(y, M, T)
            return gamma_p_rel(y, M, T) + Broken
        else:
            return gamma_p_rel(y, M, T)
    
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
    return vev(T)**2/(2*T**2) * (Tsph/T) * av_one_over_y0(M/T)

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

    