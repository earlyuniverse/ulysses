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
    
from ulysses.numba import jit
@jit
def fast_RHS(y, z, Fmat,  d, w1, n1eq, epstt,epsmm,epsee):

# constants
    c0LPM  = 4.22
    c1aLPM = 3.56
    c1bLPM = 4.77
    cQ0    = 2.57
    cQ1a   = 3.10
    cQ1b   = 2.27
    cV0    = 3.17
    cV1a   = 3.83
    cV1b   = 2.89
    g      = 0.652
    gp     = 0.357
    ht     = 0.9888
    phi0   = 0.106482
    phi1a  = 0.114281
    phi1b  = 0.0525642
    phit0  = (0.00855458 * z * z)/ (Tew * Tew)
    phit1a = (0.202923 * z * z)/ (Tew * Tew)
    phitb  = (0.101461 * z * z)/ (Tew * Tew)

  
# RHS matrices
    chi_mat     = -1./711. * np.matrix([[257,  20,  20], [20,  257,  20], [20, 20, 257 ]], dtype=np.real128)
    RN_mat      = np.matrix([[y[0], y[1]], [y[2], y[3]]], dtype=np.complex128)
    RNb_mat     = np.matrix([[y[4], y[5]], [y[6], y[7]]], dtype=np.complex128)
    mud_mat     = np.matrix([[y[8],  0,  0], [0,  y[9],  0], [0, 0, y[10] ]], dtype=np.complex128)
    mu_mat      = 2 * chi_mat @ mud_mat
    WN_mat      = (0.057018 * M0)/Tew * FdF
    WNLNV_mat   = (0.057018 * M0)/Tew * M1 * M1 * FdF

    rhs2 =    (epstt+epsmm+epsee)*d*(N1-n1eq)-w1*NBL

    return [rhs1, rhs2]


