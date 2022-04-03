# ARS leptogenesis
from scipy.special import kn
import ulysses
import numpy as np
from odeintw import odeintw


# global constants (masses in GeV)
Tew   = 131.7
gss   = 107.75
M0    = 7.1e17
zeta3 = 1.2020569

def f_TSM(z):
    return Tew/z
    
def f_ss(z):
    return (2 * np.pi  * np.pi * gss * f_TSM(z)**3)/ 45.
    
def f_HH(z):
    return (f_TSM(z) * f_TSM(z))/M0
    
def f_nphieqSM(z):
    return f_TSM(z)**3/(np.pi * np.pi)
    
def f_YHeqSM(z): # this is simply a constant ask Brian
    return (2 * f_nphieqSM(z) )/ f_ss(z)
    
def f_nNeq(M ,z):
    return (M * M * Tew * kn(2, M * z /Tew)) / (2. * np.pi * np.pi * z )
    
def f_YNeq(M z):
    return f_nNeq(M ,z)/ f_ss(z)
    
def f_Yieldeq(M, z):
    return (45. / (4. * np.pi**4 * gss) * (M * z / Tew) * (M * z / Tew)) * kn(2, M * z /Tew)
    
def f_convertmutoY(x):
    return (x * 90.) /(12 * np.pi**2 * gss)

def f_convertYBLtoYB(x):
    return x * 28./79.


