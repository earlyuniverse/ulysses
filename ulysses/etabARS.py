# ARS leptogenesis
from scipy.special import kn
import ulysses
import numpy as np
from odeintw import odeintw
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy.interpolate import interp1d, RectBivariateSpline
from scipy.integrate import quad, ode, solve_ivp, odeint
from scipy.special import zeta
import math
import os
plt.rcParams['text.usetex'] = True
from termcolor import colored
#import numba
#---------------------------------------------#
#              Auxiliary Functions            #
#---------------------------------------------#


def f_TSM(z, Tew):
    return Tew/z

def f_ss(z, Tew, gss):
    return (2 * np.pi  * np.pi * gss * f_TSM(z, Tew)**3)/ 45.

def f_HH(z, Tew, M0):
    return (f_TSM(z, Tew) * f_TSM(z, Tew))/M0

def f_nphieqSM(z, Tew):
    return f_TSM(z, Tew)**3/(np.pi * np.pi)

def f_YHeqSM(z, Tew, gss):
    return (2 * f_nphieqSM(z, Tew) )/f_ss(z, Tew, gss)

def f_nNeq(M, z, Tew):
    temp = M * z /Tew
    return (M * M * Tew * kn(2, temp.real)) / (2. * np.pi * np.pi * z )

def f_YNeq(M, z, Tew, gss):
    return f_nNeq(M, z, Tew)/ f_ss(z, Tew, gss)

def f_Yieldeq(M, z, Tew, gss):
    temp = M * z /Tew
    return (45. / (4. * np.pi**4 * gss) * (M * z / Tew) * (M * z / Tew)) * kn(2, temp.real)

def f_convertmutoY(z, gss):
    return (z * 90.) /(12 * np.pi**2 * gss)

def f_convertYBLtoYB(z):
    return z * 28./79.

def f_DYNeq(M, x, Tew, gss):
    temp = M * x /Tew
    const = gss * np.pi*np.pi*np.pi*np.pi * Tew*Tew
    kn1, kn2, kn3 = kn([1,2,3], temp.real) # bessel values, a bit faster to compute them like this
    mathematicaoutput =  (45. * M*M * x * kn2 )/(2. * const)  - (45. * M**3 * x*x * (- kn1 - kn3   ))/  -   (8. * const *Tew)
    # mathematicaoutput =  (45. * M*M * x * kn(2, temp.real))/(2. * gss * np.pi**4 * Tew*Tew)  - (45. * M**3 * x*x * (-kn(1 , temp.real) - kn(3 , temp.real)))/  -   (8. * gss * np.pi**4 * Tew*Tew*Tew)
    return mathematicaoutput

def g_run(T):

    g  = 0.652
    MZ = 91.19

    inv_g2 = 1./(g*g) + (19./(48.*np.pi**2))*math.log(np.pi*T/MZ)

    return 1/np.sqrt(inv_g2)

def gp_run(T):

    gp = 0.357
    MZ = 91.19

    inv_gp2 = 1./(gp * gp) - 0.08654517769449685 * math.log(np.pi*T/MZ)

    return 1/math.sqrt(inv_gp2)

def commutator(X, Y):
    return X @ Y - Y @ X

def anticommutator(X, Y):
    return X @ Y + Y @ X


from numba import jit, njit

@njit
def explicit_anticommutator(X, Y, R):
    """
    Compute anti commutator of 2x2 matrices X and Y, store result in 2x2 matrix R
    """
    R[0][0] = X[0][0]*Y[0][0] + X[0][1]*Y[1][0] + Y[0][0]*X[0][0] + Y[0][1]*X[1][0]
    R[0][1] = X[0][0]*Y[0][1] + X[0][1]*Y[1][1] + Y[0][0]*X[0][1] + Y[0][1]*X[1][1]
    R[1][0] = X[1][0]*Y[0][0] + X[1][1]*Y[1][0] + Y[1][0]*X[0][0] + Y[1][1]*X[1][0]
    R[1][1] = X[1][0]*Y[0][1] + X[1][1]*Y[1][1] + Y[1][0]*X[0][1] + Y[1][1]*X[1][1]

@njit
def explicit_anticommutator_array(L, R, r):
    """
    Compute 8 anticommutators of 2x2 matrics stores in 16x2 arrays
    """
    for i in range(8):
        explicit_anticommutator(L[2*i:2*(i+1)], R[2*i:2*(i+1)], r[2*i:2*(i+1)])

def diagdiag(mat):
    arr        = np.identity(3, dtype=np.complex128)
    arr[0][0]  = mat.item(0,0)
    arr[1][1]  = mat.item(1,1)
    arr[2][2]  = mat.item(2,2)
    return  arr

#-------------------------------------------------------------------------------------#
#                    Quantum Kinetic Equations for ARS Leptogenesis                   #
#-------------------------------------------------------------------------------------#

def fast_RHS(z, y, Fmat, M2, deltaM, Tew, gss, M0, M_mat, Dm2_mat, chi_mat, Lvec, Rvec, acr, funcs, use_interpolation=True):

    if use_interpolation:

        G0_fun, G1_fun, G2_fun, S0_M_fun, S1_M_fun, S2_M_fun = funcs

        G0 = G0_fun(z)
        G1 = G1_fun(z)
        G2 = G2_fun(z)
        S0 = S0_M_fun(M2, z) * (z*z/(Tew*Tew))
        S1 = S1_M_fun(M2, z) * (z*z/(Tew*Tew))
        S2 = S2_M_fun(M2, z) * (z*z/(Tew*Tew))
    else:
        G0 = 0.01007
        G1 = 0.00547
        G2 = -0.00252
        S0 = 0.0402 *   (z*z/(Tew*Tew))
        S1 = 0.00735 *  (z*z/(Tew*Tew))
        S2 = -0.01616 * (z*z/(Tew*Tew))


    T = Tew/z

    cons_1 = 0.05701803240084191
    cons_2 = 0.5480722270510788

    # RHS matrices
    RN_mat      =  np.array([[y[0], y[1]], [y[2], y[3]]], dtype=np.complex128)
    RNb_mat     =  np.array([[y[4], y[5]], [y[6], y[7]]], dtype=np.complex128)
    
    #Vector of the lepton chemical potentials
    mud_vec     =  np.array([y[8], y[9], y[10]], dtype =np.complex128)
    mu_vec      = 2 * chi_mat @ mud_vec
    #Diagonal matrix of the lepton chemical potentials
    mu_mat     =  np.diag([mu_vec[0], mu_vec[1], mu_vec[2]])

    # matrices appearing in Eqs
    FmatH               = np.transpose(np.conjugate(Fmat))
    FmatT               = np.transpose(Fmat)
    FmatC               = np.conjugate(Fmat)
    FdF                 = np.transpose(np.conjugate(Fmat)) @ Fmat

    Ham_RN              = cons_1 * (FdF + 4.* (z*z/(Tew*Tew)) * Dm2_mat)
    Fdagger_mu_F        = FmatH @ mu_mat @ Fmat
    M_Ftrans_Fstar_M    = M_mat @ FmatT @ FmatC @ M_mat
    M_Ftrans_mu_Fstar_M = M_mat @ FmatT @ mu_mat @ FmatC @ M_mat

    Ham_RNb             = cons_1 * (np.conjugate(FdF) + 4.* (z*z/(Tew*Tew)) * Dm2_mat)
    Ftrans_mu_Fstar     = FmatT @ mu_mat @ FmatC
    M_Fdagger_F_M       = M_mat @ FmatH @ Fmat @ M_mat
    M_Fdagger_mu_F_M    = M_mat @ FmatH @ mu_mat @ Fmat @ M_mat

    F_RN_Fdagger        = Fmat @ RN_mat @ FmatH
    Fstar_RNb_Ftrans    = FmatC @ RNb_mat @ FmatT
    F_Fdagger           = Fmat @ FmatH
    Fstar_M_RN_M_Ftrans = FmatC @ M_mat @ RN_mat @ M_mat @ FmatT
    F_M_RNb_M_Fdagger   = Fmat @ M_mat @ RNb_mat @ M_mat @ FmatH
    F_M_M_Fdagger       = Fmat @ M_mat @ M_mat @ FmatH
    RNmatmId2           = RN_mat  - np.identity(2)
    RNbmatmId2          = RNb_mat - np.identity(2)

    # Structurally, we have 8 calls here to anticommutator with 2x2 matrices as arguments
    # Let's  define left hand, right hand and results vectors of shape (16,2) once and explicitly compute elements:
    # Build the arrays first

    Lvec[0:2]   = FdF
    Lvec[2:4]   = Fdagger_mu_F
    Lvec[4:6]   = M_Ftrans_Fstar_M
    Lvec[6:8]   = M_Ftrans_mu_Fstar_M
    Lvec[8:10]  = np.conjugate(Lvec[0:2])
    Lvec[10:12] = Ftrans_mu_Fstar
    Lvec[12:14] = M_Fdagger_F_M
    Lvec[14:16] = M_Fdagger_mu_F_M

    Rvec[0:2]   = RNmatmId2
    Rvec[2:4]   = RN_mat
    Rvec[4:6]   = RNmatmId2
    Rvec[6:8]   = RN_mat
    Rvec[8:10]  = RNbmatmId2
    Rvec[10:12] = RNb_mat
    Rvec[12:14] = RNbmatmId2
    Rvec[14:16] = RNb_mat

    # Compute anticommutators
    explicit_anticommutator_array(Lvec, Rvec, acr)

    # This is expensive so cache it
    fdyneq = f_DYNeq(M2, z, Tew, gss)
    fyneq = f_YNeq(M2, z, Tew, gss)

    # ARS Equations
    RNRHS_mat   = (M0/Tew) * (- 1j * commutator(Ham_RN, RN_mat) - 0.5 * G0 * acr[0:2]
                              + G1 * Fdagger_mu_F
                              - 0.5 * G2 * acr[2:4]
                              - 0.5 * S0 * acr[4:6]
                              - S1 * M_Ftrans_mu_Fstar_M
                              + 0.5 * S2 * acr[6:8]
                              - (Tew/M0) * (RN_mat/ fyneq) * fdyneq
                              )


    RNbRHS_mat  = (M0/Tew) * (- 1j * commutator(Ham_RNb, RNb_mat)
                              - 0.5 * G0 * acr[8:10]
                              - G1 * Ftrans_mu_Fstar
                              + 0.5 * G2 * acr[10:12]
                              - 0.5 * S0 * acr[12:14]
                              + S1 * M_Fdagger_mu_F_M
                              - 0.5 * S2 * acr[14:16]
                              - (Tew/M0) * (RNb_mat/ fyneq) * fdyneq
                              )


    eqtns = np.zeros(11, dtype=np.complex128)

    eqtns[0]  = RNRHS_mat[0][0]
    eqtns[1]  = RNRHS_mat[0][1]
    eqtns[2]  = RNRHS_mat[1][0]
    eqtns[3]  = RNRHS_mat[1][1]

    eqtns[4]  = RNbRHS_mat[0][0]
    eqtns[5]  = RNbRHS_mat[0][1]
    eqtns[6]  = RNbRHS_mat[1][0]
    eqtns[7]  = RNbRHS_mat[1][1]

    # This is just a 3 vector and we only care about diagonal elements of the 3x3 matrices
    for i in range(3):
        eqtns[8+i] = cons_2 * (M0/Tew) * (
                - 0.5 * G0 * (F_RN_Fdagger[i][i] - Fstar_RNb_Ftrans[i][i])
                +       G1 * F_Fdagger[i][i]* mu_vec[i]
                - 0.5 * G2 * (F_RN_Fdagger[i][i] + Fstar_RNb_Ftrans[i][i]) * mu_vec[i]
                + 0.5 * S0 * (Fstar_M_RN_M_Ftrans[i][i] - F_M_RNb_M_Fdagger[i][i])
                +       S1 * F_M_M_Fdagger[i][i] * mu_vec[i]
                - 0.5 * S2 * (Fstar_M_RN_M_Ftrans[i][i] + F_M_RNb_M_Fdagger[i][i]) * mu_vec[i]
            )


    return eqtns

def fast_averaged_RHS(z, y, Fmat, M2, deltaM, Tew, gss, M0, M_mat, Dm2_mat, chi_mat, Lvec, Rvec, acr, funcs, use_interpolation=True):

    if use_interpolation:

        G0_fun, G1_fun, G2_fun, S0_M_fun, S1_M_fun, S2_M_fun = funcs

        G0 = G0_fun(z)
        G1 = G1_fun(z)
        G2 = G2_fun(z)
        S0 = S0_M_fun(M2, z) * (z*z/(Tew*Tew))
        S1 = S1_M_fun(M2, z) * (z*z/(Tew*Tew))
        S2 = S2_M_fun(M2, z) * (z*z/(Tew*Tew))
    else:
        G0 = 0.01007
        G1 = 0.00547
        G2 = -0.00252
        S0 = 0.0402 * (z*z/(Tew*Tew))
        S1 = 0.00735 * (z*z/(Tew*Tew))
        S2 = -0.01616 * (z*z/(Tew*Tew))

    M_mat[0][0]   = M2
    M_mat[1][1]   = M2 + deltaM

    deltaM2 = 2 * M2 * deltaM + deltaM * deltaM

    cons_1 = 0.05701803240084191
    cons_2 = 0.5480722270510788

    RN_mat      =  np.diag([y[0], y[1]])
    RNb_mat     =  np.diag([y[2], y[3]])
    
    #Vector of the lepton chemical potentials
    mud_vec     =  np.array([y[8], y[9], y[10]], dtype =np.complex128)
    mu_vec      = 2 * chi_mat @ mud_vec
    #Diagonal matrix of the lepton chemical potentials
    mu_mat     =  np.diag([mu_vec[0], mu_vec[1], mu_vec[2]])
    

    # matrices appearing in Eqs
    FmatH               = np.transpose(np.conjugate(Fmat))
    FmatT               = np.transpose(Fmat)
    FmatC               = np.conjugate(Fmat)
    FdF                 = np.transpose(np.conjugate(Fmat)) @ Fmat

    Fdagger_mu_F        = FmatH @ mu_mat @ Fmat
    M_Ftrans_Fstar_M    = M_mat @ FmatT @ np.conjugate(Fmat) @ M_mat
    M_Ftrans_mu_Fstar_M = M_mat @ FmatT @ mu_mat @ np.conjugate(Fmat) @ M_mat

    Ftrans_mu_Fstar     = FmatT @ mu_mat @ np.conjugate(Fmat)
    M_Fdagger_F_M       = M_mat @ FmatH @ Fmat @ M_mat
    M_Fdagger_mu_F_M    = M_mat @ FmatH @ mu_mat @ Fmat @ M_mat

    F_RN_Fdagger        = Fmat @ RN_mat @ FmatH
    Fstar_RNb_Ftrans    = np.conjugate(Fmat) @ RNb_mat @ FmatT
    F_Fdagger           = Fmat @ FmatH
    Fstar_M_RN_M_Ftrans = np.conjugate(Fmat) @ M_mat @ RN_mat @ M_mat @ FmatT
    F_M_RNb_M_Fdagger   = Fmat @ M_mat @ RNb_mat @ M_mat @ FmatH
    F_M_M_Fdagger       = Fmat @ M_mat @ M_mat @ FmatH
    RNmatmId2           = RN_mat  - np.identity(2)
    RNbmatmId2          = RNb_mat - np.identity(2)


    Lvec[0:2]   = FdF
    Lvec[2:4]   = Fdagger_mu_F
    Lvec[4:6]   = M_Ftrans_Fstar_M
    Lvec[6:8]   = M_Ftrans_mu_Fstar_M
    Lvec[8:10]  = np.conjugate(Lvec[0:2])
    Lvec[10:12] = Ftrans_mu_Fstar
    Lvec[12:14] = M_Fdagger_F_M
    Lvec[14:16] = M_Fdagger_mu_F_M

    Rvec[0:2]   = RNmatmId2
    Rvec[2:4]   = RN_mat
    Rvec[4:6]   = RNmatmId2
    Rvec[6:8]   = RN_mat
    Rvec[8:10]  = RNbmatmId2
    Rvec[10:12] = RNb_mat
    Rvec[12:14] = RNbmatmId2
    Rvec[14:16] = RNb_mat

    # Compute anticommutators
    explicit_anticommutator_array(Lvec, Rvec, acr)

    # This is expensive so cache it
    fdyneq = f_DYNeq(M2, z, Tew, gss)
    fyneq = f_YNeq(M2, z, Tew, gss)

    # ARS Equations
    RNRHS_mat   = (M0/Tew) * (- 0.5 * G0 * acr[0:2]
                              + G1 * Fdagger_mu_F
                              - 0.5 * G2 * acr[2:4]
                              - 0.5 * S0 * acr[4:6]
                              - S1 * M_Ftrans_mu_Fstar_M
                              + 0.5 * S2 * acr[6:8]
                              - (Tew/M0) * (RN_mat/ fyneq) * fdyneq
                              )

    RNbRHS_mat  = (M0/Tew) * (- 0.5 * G0 * acr[8:10]
                              - G1 * Ftrans_mu_Fstar
                              + 0.5 * G2 * acr[10:12]
                              - 0.5 * S0 * acr[12:14]
                              + S1 * M_Fdagger_mu_F_M
                              - 0.5 * S2 * acr[14:16]
                              - (Tew/M0) * (RNb_mat/ fyneq) * fdyneq
                              )

    eqtns = np.zeros(7, dtype=np.complex128)

    eqtns[0] = RNRHS_mat[0,0]
    eqtns[1] = RNRHS_mat[1,1]

    eqtns[2] = RNbRHS_mat[0,0]
    eqtns[3] = RNbRHS_mat[1,1]

    for i in range(3):
        eqtns[4+i] = cons_2 * (M0/Tew) * (
                - 0.5 * G0 * (F_RN_Fdagger[i][i] - Fstar_RNb_Ftrans[i][i])
                +       G1 * F_Fdagger[i][i]* mu_vec[i]
                - 0.5 * G2 * (F_RN_Fdagger[i][i] + Fstar_RNb_Ftrans[i][i]) * mu_vec[i]
                + 0.5 * S0 * (Fstar_M_RN_M_Ftrans[i][i] - F_M_RNb_M_Fdagger[i][i])
                +       S1 * F_M_M_Fdagger[i][i] * mu_vec[i]
                - 0.5 * S2 * (Fstar_M_RN_M_Ftrans[i][i] + F_M_RNb_M_Fdagger[i][i]) * mu_vec[i]
            )

    return eqtns


class EtaB_ARS(ulysses.ULSBase):
    """
    add description of where to find BEs
    """

    def shortname(self): return "BEARS"

    def flavourindices(self): return [1]

    def flavourlabels(self): return ["$NBL$"]

    def RHS(self, z, y, Fmat, M2, deltaM, Tew, gss, M0, M_mat, Dm2_mat, chi_mat, Lvec, Rvec, acr):

        funcs = []
        return fast_RHS(z, y, Fmat, M2, deltaM, Tew, gss, M0, M_mat, Dm2_mat, chi_mat, Lvec, Rvec, acr, funcs, use_interpolation=False)

    def RHS_averaged(self, z, y, Fmat, M2, deltaM, Tew, gss, M0, M_mat, Dm2_mat, chi_mat, Lvec, Rvec, acr):

        funcs = []
        return fast_averaged_RHS(z, y, Fmat, M2, deltaM, Tew, gss, M0, M_mat, Dm2_mat, chi_mat, Lvec, Rvec, acr, funcs, use_interpolation=False)

    @property
    def EtaB(self):
        
        # global constants (masses in GeV)
        Tew    = 131.7
        gss    = 106.75
        M0     = 7.e17
        zeta3  = zeta(3.)
        gstaro = 43/11 #entropic effective degrees of freedom at present

        mp       = 1.672621898*1e-24  #proton mass in g
        ngamma   = 410.7              #present photon number density in cm^-3
        rhoc     = 1.87840*1e-29      #critical density of the Universe in h^2 g cm^-3
        gstaro   = 43/11              #entropic effective degrees of freedom at present
        ToYb     = 45 * zeta3 /(gstaro * np.pi**4)
        ToOmegab = mp * ngamma/rhoc

        self.M1 = 0.
        self.m  = 0.

        Fmat     = np.delete(self.h, 0, 1)

        if 30. <= self.M2 and self.M2 < 100.:
            print("\n")
            print(colored("Warning: RHN masses sufficiently heavy that non-relativistic corrections to lepton-number-violating rates can be important at the time of sphaleron decoupling.\n", 'red'))

        if self.M2 >= 100.:
            print("\n")
            print(colored("Warning: Our implementation of the kinetic equations assumes that RHNs are at least somewhat relativistic. \n\t Results are not validated and may be incorrect for this set of RHN masses. \n", 'red'))


        # FFd = Fmat @ Fmat.H
        FFd = Fmat @ np.transpose(np.conjugate( Fmat))

        G1_t = 0.00547 
        S1_t = 0.00735/(Tew*Tew)

        G1_Hub = 0.5*(M0/Tew) * G1_t * np.diag(FFd).real
        G2_Hub = 0.5*(M0/Tew) * S1_t * np.diag(FFd).real

        if np.min(G1_Hub) >= 1. and  np.min(G2_Hub) >= 1.:
            print("\n")
            print(colored("Warning: This set of parameters is in the strong washout regime for lepton asymmetries. \n", 'red'))

        dMval = self.M3 - self.M2
       

        params  = np.array([Fmat[0, 0], Fmat[0, 1], Fmat[1, 0], Fmat[1, 1], Fmat[2, 0], Fmat[2, 1], self.M2, dMval, Tew, gss, M0],
                           dtype=np.complex128)

        # initial conditions in the order RN11, RN12, RN21, RN22, RNb11, RNb12, RNb21, RNb22, mudelta1, mudelta2,  mudelta3
        YNin = 0.+0j
        y0 = np.array([YNin, 0+0j, 0+0j, YNin, YNin, 0+0j, 0+0j, YNin, 0+0j, 0+0j, 0+0j], dtype=np.complex128)

        M_mat       = np.zeros([2,2])
        M_mat[0][0] = self.M2
        M_mat[1][1] = self.M2 + dMval
        deltaM2     = 2 * self.M2 * dMval + dMval * dMval
        Dm2_mat     = np.array([[0 , 0], [0 , deltaM2 ]])

        chi_mat  =  -1./711. * np.array([[257.,  20.,  20.], [20.,  257.,  20.], [20., 20., 257. ]], dtype=np.complex128)
        Lvec     = np.zeros((16,2), dtype=np.complex128)
        Rvec     = np.zeros((16,2), dtype=np.complex128)
        acr      = np.zeros((16,2), dtype=np.complex128)

        #Dimensionless oscillation time, see 2109.10908, Eq 6.
        zosc = np.cbrt(12.*Tew**3/(deltaM2 * M0)) 

        # if zosc > 0.1, oscillations start late so the full integration [1e-6,1] can be performed without problems
        if zosc > 0.1:
        
            print(colored("Info: zosc = {}".format(zosc), 'red'))

            ys = solve_ivp(lambda t, z: self.RHS(t, z, Fmat, self.M2, dMval, Tew, gss, M0, M_mat, Dm2_mat, chi_mat, Lvec, Rvec, acr),
                           [1.e-6, 1], y0, method='BDF', rtol=1e-7, atol=1e-10)
            
            t, muD1, muD2, muD3 = [ys.t,ys.y[8], ys.y[9], ys.y[10]]
        
        # if zosc < 0.1, oscillations start early. The numerical integration could be slow. The default behaviour is that the full integration is performed
        # unless the user provides a zcut  value (this is the point i the z-integration where the stitching occurs and the off diagonal terms are ignored)
        # see below for more details
        else:

            print(colored("Info: zosc = {}, stitching occurs at {}.".format(zosc,self._zcut), 'red'))
            
            if self._zcut == 1:
            
                    print(colored("Info: zcut is at default value 1, consider inputting  0.5 < zcut < 1 if integration is slow.".format(zosc,self._zcut), 'red'))

            ys_1 = solve_ivp(lambda t, z: self.RHS(t, z, Fmat, self.M2, dMval, Tew, gss, M0, M_mat, Dm2_mat, chi_mat, Lvec, Rvec, acr),
                                 [1.e-6, self._zcut], y0, method='BDF', rtol=1.e-5, atol=1.e-5)
            
            t_1, muD1_1, muD2_1, muD3_1 = [ys_1.t,ys_1.y[8], ys_1.y[9], ys_1.y[10]]

                #  this is where the stitch occurs if the user provides a zcut
                # For the stitching case, we ignore the evolution of the off-diagonal terms in the RN and RNb matrices.
                # Thus, we only take the solutions related to RN11, RN22, RNb11, RNb22 for the averaged equations

            y0_2  = np.array([np.abs(ys_1.y[0,-1]),  np.abs(ys_1.y[3,-1]),  np.abs(ys_1.y[4,-1]), np.abs(ys_1.y[7,-1]),
                                  np.real(ys_1.y[8,-1]), np.real(ys_1.y[9,-1]), np.real(ys_1.y[10,-1])], dtype=np.complex128)

                

            solARS = solve_ivp(lambda t, z: self.RHS_averaged(t, z, Fmat, self.M2, dMval, Tew, gss, M0, M_mat, Dm2_mat, chi_mat, Lvec, Rvec, acr), [self._zcut, 1], y0_2, method='BDF', rtol=1.e-5, atol=1.e-10)

            t_2, muD1_2, muD2_2, muD3_2 = [solARS.t, solARS.y[4], solARS.y[5], solARS.y[6]]
            t,   muD1,    muD2, muD3 = [np.concatenate((t_1,t_2), axis=0), np.concatenate((muD1_1, muD1_2) , axis=0), np.concatenate((muD2_1, muD2_2), axis=0), np.concatenate((muD3_1, muD3_2), axis=0)]
                
        YB_sol  = np.real(muD1[-1] + muD2[-1] + muD3[-1])
        
        plt.plot(t, np.abs(muD1), label=r"$|\mu_e|$")
        plt.plot(t, np.abs(muD2), label=r"$|\mu_\mu|$")
        plt.plot(t, np.abs(muD3), label=r"$|\mu_\tau|$")
        plt.xlabel(r"$T_{\rm{ew}}/T$", fontsize=16)#z=T_{ew}/T in this module
        plt.legend(loc='lower right', fontsize=16)
        plt.ylabel(r"$|\mu|$",  fontsize=16)
        plt.show()
        YB  = f_convertYBLtoYB(f_convertmutoY(YB_sol, gss))
        etaB= YB/ToYb
    
        return etaB


class EtaB_ARS_INTERP(EtaB_ARS):

    def shortname(self): return "BEARS_INTERP"

    def flavourindices(self): return [1]

    def flavourlabels(self): return ["$NBL$"]

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        #---------------------------------------------#
        #                 Integrated Rates            #
        #---------------------------------------------#

        
        data_dir = os.path.dirname(ulysses.__file__)
      
        G0_f = os.path.join(data_dir, "./data/g0log.dat")
        G1_f = os.path.join(data_dir, "./data/g1log.dat")
        G2_f = os.path.join(data_dir, "./data/g2log.dat")
        

        G0Tab = np.loadtxt(G0_f, skiprows=0)
        G1Tab = np.loadtxt(G1_f, skiprows=0)
        G2Tab = np.loadtxt(G2_f, skiprows=0)

        self.G0Int_ = interpolate.splrep(G0Tab[:,0], G0Tab[:,1], s=0)
        self.G1Int_ = interpolate.splrep(G1Tab[:,0], G1Tab[:,1], s=0)
        self.G2Int_ = interpolate.splrep(G2Tab[:,0], G2Tab[:,1], s=0)

        M2_f  = os.path.join(data_dir, "./data/Log_M.txt")
        z2_f  = os.path.join(data_dir, "./data/Log_z.txt")
        
        S0M_f = os.path.join(data_dir, "./data/s0log_massdep.txt")
        S1M_f = os.path.join(data_dir, "./data/s1log_massdep.txt")
        S2M_f = os.path.join(data_dir, "./data/s2log_massdep.txt")

        M2tab  = np.loadtxt(M2_f, skiprows=0)
        z2tab  = np.loadtxt(z2_f, skiprows=0)
        
        #Mass-dependent G0, G1, G2 functions. Use these to include the sub-leading non-relativistic contributions to G0, G1 and G2.
        

        G0M_f = os.path.join(data_dir, "./data/g0log_massdep.txt")
        G1M_f = os.path.join(data_dir, "./data/g1log_massdep.txt")
        G2M_f = os.path.join(data_dir, "./data/g2log_massdep.txt")
        
      
     
        G0Mtab = np.loadtxt(G0M_f, skiprows=0)
        G1Mtab = np.loadtxt(G1M_f, skiprows=0)
        G2Mtab = np.loadtxt(G2M_f, skiprows=0)
        
        S0Mtab = np.loadtxt(S0M_f, skiprows=0)
        S1Mtab = np.loadtxt(S1M_f, skiprows=0)
        S2Mtab = np.loadtxt(S2M_f, skiprows=0)

        self.S0MInt_ = RectBivariateSpline(M2tab, z2tab, S0Mtab) # 2-D Interpolation
        self.S1MInt_ = RectBivariateSpline(M2tab, z2tab, S1Mtab) # 2-D Interpolation
        self.S2MInt_ = RectBivariateSpline(M2tab, z2tab, S2Mtab) # 2-D Interpolation
        
        self.G0MInt_ = RectBivariateSpline(M2tab, z2tab, G0Mtab) # 2-D Interpolation
        self.G1MInt_ = RectBivariateSpline(M2tab, z2tab, G1Mtab) # 2-D Interpolation
        self.G2MInt_ = RectBivariateSpline(M2tab, z2tab, G2Mtab) # 2-D Interpolation
 
    def G0_fun(self,z): return interpolate.splev(math.log(z), self.G0Int_, der=0)

    def G1_fun(self,z): return interpolate.splev(math.log(z), self.G1Int_, der=0)

    def G2_fun(self,z): return interpolate.splev(math.log(z), self.G2Int_, der=0)

    def S0_M_fun(self,M, z):

        if M < 0.1:
            M = 0.1
        elif M > 100.:
            M = 100.

        return self.S0MInt_(np.log(M), np.log(z))[0,0]

    def S1_M_fun(self,M, z):

        if M < 0.1:
            M = 0.1
        elif M > 100.:
            M = 100.

        return self.S1MInt_(np.log(M), np.log(z))[0,0]

    def S2_M_fun(self,M, z):

        if M < 0.1:
            M = 0.1
        elif M > 100.:
            M = 100.

        return self.S2MInt_(np.log(M), np.log(z))[0,0]
    
      
    def G0_M_fun(self,M, z):

        if M < 0.1:
            M = 0.1
        elif M > 100.:
            M = 100.

        return self.G0MInt_(np.log(M), np.log(z))[0,0]

    def G1_M_fun(self,M, z):

        if M < 0.1:
            M = 0.1
        elif M > 100.:
            M = 100.

        return self.G1MInt_(np.log(M), np.log(z))[0,0]

    def G2_M_fun(self,M, z):

        if M < 0.1:
            M = 0.1
        elif M > 100.:
            M = 100.

        return self.G2MInt_(np.log(M), np.log(z))[0,0]


    def RHS(self, z, y, Fmat, M2, deltaM, Tew, gss, M0, M_mat, Dm2_mat, chi_mat, Lvec, Rvec, acr):

        funcs = [self.G0_fun, self.G1_fun, self.G2_fun, self.S0_M_fun, self.S1_M_fun, self.S2_M_fun]

        return fast_RHS(z, y, Fmat, M2, deltaM, Tew, gss, M0, M_mat, Dm2_mat, chi_mat, Lvec, Rvec, acr, funcs, use_interpolation=True)

    def RHS_averaged(self, z, y, Fmat, M2, deltaM, Tew, gss, M0, M_mat, Dm2_mat, chi_mat, Lvec, Rvec, acr):

        funcs = [self.G0_fun, self.G1_fun, self.G2_fun, self.S0_M_fun, self.S1_M_fun, self.S2_M_fun]

        return fast_averaged_RHS(z, y, Fmat, M2, deltaM, Tew, gss, M0, M_mat, Dm2_mat, chi_mat, Lvec, Rvec, acr, funcs, use_interpolation=True)
