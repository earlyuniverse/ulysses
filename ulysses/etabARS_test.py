# ARS leptogenesis
from scipy.special import kn
import ulysses
import numpy as np
from odeintw import odeintw
import matplotlib.pyplot as plt


# global constants (masses in GeV)
Tew   = 131.7
gss   = 107.75
M0    = 7.1e17
zeta3 = 1.20206

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
    temp = M * z /Tew
    return (M * M * Tew * kn(2, temp.real)) / (2. * np.pi * np.pi * z )
    
def f_YNeq(M, z):
    return f_nNeq(M ,z)/ f_ss(z)
    
def f_Yieldeq(M, z):
    temp = M * z /Tew
    return (45. / (4. * np.pi**4 * gss) * (M * z / Tew) * (M * z / Tew)) * kn(2, temp.real)
    
def f_convertmutoY(z):
    return (z * 90.) /(12 * np.pi**2 * gss)

def f_convertYBLtoYB(z):
    return z * 28./79.

def f_DYNeq(M, x):
    temp = M * x /Tew
    mathematicaoutput =  (45. * M**2 * x * kn(2, temp.real))/(2. * gss * np.pi**4 * Tew**2)  - (45. * M**3 * x**2 * (-kn(1 , temp.real) - kn( 3 , temp.real)))/  -   (8. * gss * np.pi**4 * Tew**3)
    return mathematicaoutput

def commutator(X, Y):
    return X @ Y - Y @ X

def anticommutator(X, Y):
    return X @ Y + Y @ X

def matrix_diag3(d1,d2,d3):
    return np.array([[d1, 0.0, 0.0], [0.0, d2, 0.0], [0.0, 0.0, d3]])

def matrix_diag2(d1,d2):
    return np.array([[d1, 0.0], [0.0, d2]])


# Generic Rotations #
def matrix_rot23(th23):
    return np.array([[1.0,          0.0 , 0.0],
                    [0.0,  np.cos(th23), np.sin(th23)],
                    [0.0, -np.sin(th23), np.cos(th23)]])

def matrix_rot12(th12):
    return np.array([[ np.cos(th12), np.sin(th12), 0.0],
                    [-np.sin(th12), np.cos(th12), 0.0],
                    [          0.0,  0.0,         1.0]])

def matrix_rot13(th13):
    return np.array([[                     np.cos(th13), 0.0, np.sin(th13) ],
                    [                     0.0         , 1.0, 0.0   ],
                    [-np.sin(th13), 0.0, np.cos(th13)]],
                    dtype=np.complex64)

def matrix_pmns(th12, th13, th23, delta, alpha1):
    return matrix_rot23(th23) @ matrix_diag3(np.exp(-1j * delta/2),  1, np.exp(1j * delta/2)) @ matrix_rot13(th13) @  matrix_diag3(np.exp(1j * delta/2),  1, np.exp(-1j * delta/2)) @ matrix_rot12(th12) @ matrix_diag3(1, np.exp(-1j * alpha1), 1)

    

def fast_RHS(y, z, Fmat11, Fmat12,Fmat21,Fmat22,Fmat31,Fmat32, M1, deltaM):
#    Fmat11, Fmat12,Fmat21,Fmat22,Fmat31,Fmat32, M1, deltaM = params
    Fmat    = np.matrix([[Fmat11, Fmat12],[Fmat21, Fmat22],[Fmat31, Fmat32]])
# constants
    c0LPM   = 4.22
    c1aLPM  = 3.56
    c1bLPM  = 4.77
    cQ0     = 2.57
    cQ1a    = 3.10
    cQ1b    = 2.27
    cV0     = 3.17
    cV1a    = 3.83
    cV1b    = 2.89
    g       = 0.652
    gp      = 0.357
    ht      = 0.9888
    phi0    = 0.106482
    phi1a   = 0.114281
    phi1b   = 0.0525642
    phit0   = (0.00855458 * z * z)/ (Tew * Tew)
    phit1a  = (0.202923 * z * z)/ (Tew * Tew)
    phit1b  = (0.101461 * z * z)/ (Tew * Tew)
    FdF     = Fmat.H @ Fmat
    deltaM2 = 2 * M1 * deltaM + deltaM * deltaM
    TL2     = 377923. * deltaM2**(1./3.)
    r2      = TL2/Tew

# matrices
    r_mat       = np.matrix([[0,0], [0, r2**3]], dtype=np.complex128)
  
# RHS matrices
    
    RN_mat      =  np.matrix([[1, 2], [1, 1]], dtype=np.complex128)
    RNb_mat     =  np.matrix([[1, 1], [1, 1]], dtype=np.complex128)
    
    mud_mat     =  np.matrix([[3,  0,  0], [0, 1,  0], [0, 0, 1 ]], np.dtype(np.int32) )
    chi_mat     =  -1./711. * np.matrix([[257,  20,  20], [20,  257,  20], [20, 20, 257 ]], dtype=np.complex128)
 
    mu_mat      = 2 * chi_mat @ mud_mat

    WN_mat      =  (np.pi * np.pi * M0)/(144 * zeta3 * Tew) * FdF
    


    WNLNV_mat   =  (0.057018 * M0)/Tew * M1 * M1 * FdF

    
    omu_mat     =  (4.048280300459774e16 / Tew) * Fmat.H @ mu_mat @ Fmat
    

    omub_mat    = -(4.048280300459774e16 / Tew) * Fmat.T @ mu_mat @ np.conjugate(Fmat)
    
    omuLNV_mat  =  (np.pi * np.pi * M0 * M1 * M1)/(144 * zeta3 * Tew) * Fmat.H @ mu_mat @ Fmat
    
    omubLNV_mat =  omub_mat * M1 * M1
    
# matrices for muDeltaRHS
    FRNFdagger_mat     = Fmat @ RN_mat @ Fmat.H
    FstarRNbFtrans_mat = np.conjugate(Fmat) @ RNb_mat @ Fmat.T
    FFdagger           = Fmat @ Fmat.H
    

    RNRHS_mat   = 1j * commutator(RN_mat, WN_mat)     + 3j * z * z * commutator(RN_mat,r_mat)    - phit0  * anticommutator(RN_mat, WNLNV_mat)          - phi0 * anticommutator(RN_mat, WN_mat)    + 2 * phi0 * WN_mat      + 2 * phit0 * WNLNV_mat    - RN_mat/f_YNeq(M1, z) * f_DYNeq(M1, z)   + phi1a * omu_mat    - phit1a * omuLNV_mat       + 0.5 * phi1b * anticommutator(omu_mat, RN_mat)   - 0.5 * phit1b * anticommutator(omuLNV_mat, RN_mat)
    

    RNbRHS_mat  = 1j * commutator(RNb_mat, WN_mat.T)  + 3j * z * z * commutator(RNb_mat,r_mat)  - phit0  * anticommutator(RNb_mat, WNLNV_mat.T)   - RNb_mat/f_YNeq(M1, z) * f_DYNeq(M1, z)     - phi0 * anticommutator(RNb_mat, WN_mat.T)  + 2 * phi0 * WN_mat.T    + 2 * phit0 * WNLNV_mat.T   + phi1a * omub_mat    - phit1a * omubLNV_mat    - 0.5 * phit1b * anticommutator(omubLNV_mat, RNb_mat)    + 0.5 * phi1b * anticommutator(omub_mat, RNb_mat)

#    muDeltaRHS  = M0/(32 * Tew) * (  - phi0  * (FRNFdagger_mat - FstarRNbFtrans_mat)                                                                                +  M1 * M1 * phit0 * (FRNFdagger_mat - FstarRNbFtrans_mat) + phi1a *    FFdagger @ mu_mat ).diagonal()
##      + phi1a *   FFdagger @ mu_mat  + phi1a *  mu_mat @  FFdagger   + phi1b/2.* mu_mat @ (FRNFdagger_mat + FstarRNbFtrans_mat) - M1 * M1 * phit1a * mu_mat @ FFdagger            - 0.5 * M1 * M1 * phit1b  * mu_mat @ (FRNFdagger_mat + FstarRNbFtrans_mat)

    muDeltaRHS  = M0/(32 * Tew) * (   -phi0  * (FRNFdagger_mat - FstarRNbFtrans_mat).diagonal()    + phi1a * (np.diag(np.diag(FFdagger))  @ np.diag(np.diag( mu_mat)).diagonal() )             + phi1b/2. * (np.diag(np.diag(FRNFdagger_mat + FstarRNbFtrans_mat))  @ np.diag(np.diag( mu_mat))).diagonal() + M1 * M1 * phit0 * ( FRNFdagger_mat - FstarRNbFtrans_mat).diagonal() - M1 * M1 * phit1a *  (np.diag(np.diag(FFdagger))  @ np.diag(np.diag( mu_mat)).diagonal() )   - 0.5 * M1 * M1 * phit1b  * (np.diag(np.diag(FRNFdagger_mat + FstarRNbFtrans_mat))  @ np.diag(np.diag( mu_mat))).diagonal())


    
    
    return [RNRHS_mat,RNbRHS_mat,muDeltaRHS,mu_mat]
    
class EtaB_ARS(ulysses.ULSBase):
    """
    add description of where to find BEs
    """

    def RHS(self, y, z, Fmat11, Fmat12,Fmat21,Fmat22,Fmat31,Fmat32, M1, deltaM ):

        return fast_RHS(y, z, Fmat11, Fmat12,Fmat21,Fmat22,Fmat31,Fmat32, M1, deltaM)

    @property
    def EtaB(self):
        y0       = np.array([1+0j,0+0j, 0+0j, 1+0j, 1+0j, 0+0j, 0+0j, 1+0j, 0+0j, 0+0j, 0+0j], dtype=np.complex128)

#        params   = np.array([epstt,epsmm,epsee,k], dtype=np.complex128)
        
        # CI parameters for test
        th12val  = np.arcsin(0.557)
        th13val  = np.arcsin(0.1497)
        th23val  = np.arcsin(0.75)
        omegaval = np.pi/4 - 0.7 * 1j
        m1val    = 0.
        m2val    = 8.6e-12
        m3val    = 58e-12
        M1val    = 40.0
        dMval    = 0.2e-10
        M2val    = M1val + dMval
        deltaval = 221. * np.pi/180.
        alpha1val= np.pi/3 - deltaval
        vev      = 246.
        

        
        # construct the CI parametrsition for internal tests
        
        mnu      = matrix_diag3(m1val, m2val, m3val)
        mM       = matrix_diag2(M1val, M2val)
        R_mat    = np.matrix([ [0,0] , [np.cos(omegaval), np.sin(omegaval)] , [-np.sin(omegaval), np.cos(omegaval)] ])
        
        Fmat     = (np.sqrt(2)/ vev) * matrix_pmns(th12val, th13val, th23val, deltaval, alpha1val).conjugate() @ np.sqrt( mnu ) @ R_mat.conjugate() @ np.sqrt(mM)
        Fmat11    = Fmat[0, 0]
        Fmat12    = Fmat[0, 1]
        Fmat21    = Fmat[1, 0]
        Fmat22    = Fmat[1, 1]
        Fmat31    = Fmat[2, 0]
        Fmat32    = Fmat[2, 1]
        Fmat      = np.matrix([[Fmat11, Fmat12],[Fmat21, Fmat22],[Fmat31, Fmat32]])
        FdF       = Fmat.H @ Fmat


        from IPython import embed
        embed()
        fast_RHS(y0, 1, Fmat11, Fmat12,Fmat21,Fmat22,Fmat31,Fmat32, M1val, dMval)
#        fast_RHS(y0, 1, Fmat, M1val, dMval)
        
