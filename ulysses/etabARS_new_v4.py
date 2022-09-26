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
from termcolor import colored
import numba
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

#---------------------------------------------#
#                 Integrated Rates            #
#---------------------------------------------#

G0TabAux = np.loadtxt("./Data/g0log.dat", skiprows=0)
G1TabAux = np.loadtxt("./Data/g1log.dat", skiprows=0)
G2TabAux = np.loadtxt("./Data/g2log.dat", skiprows=0)

S0TabAux = np.loadtxt("./Data/s0log.dat", skiprows=0)
S1TabAux = np.loadtxt("./Data/s1log.dat", skiprows=0)
S2TabAux = np.loadtxt("./Data/s2log.dat", skiprows=0)

ztab  = G0TabAux[:,0]

G0tab = G0TabAux[:,1]
G0Int = interpolate.splrep(ztab, G0tab, s=0)
def G0_fun(z): return interpolate.splev(math.log(z), G0Int, der=0)

G1tab = G1TabAux[:,1]
G1Int = interpolate.splrep(ztab, G1tab, s=0)
def G1_fun(z): return interpolate.splev(math.log(z), G1Int, der=0)

G2tab = G2TabAux[:,1]
G2Int = interpolate.splrep(ztab, G2tab, s=0)
def G2_fun(z): return interpolate.splev(math.log(z), G2Int, der=0)

S0tab = S0TabAux[:,1]
S0Int = interpolate.splrep(ztab, S0tab, s=0)
def S0_fun(z): return interpolate.splev(math.log(z), S0Int, der=0)

S1tab = S1TabAux[:,1]
S1Int = interpolate.splrep(ztab, S1tab, s=0)
def S1_fun(z): return interpolate.splev(math.log(z), S1Int, der=0)

S2tab = S2TabAux[:,1]
S2Int = interpolate.splrep(ztab, S2tab, s=0)
def S2_fun(z): return interpolate.splev(math.log(z), S2Int, der=0)


M1tab  = np.loadtxt("./Data/Log_M.txt", skiprows=0)
z2tab  = np.loadtxt("./Data/Log_z.txt", skiprows=0)
S0Mtab = np.loadtxt("./Data/s0log_massdep.txt", skiprows=0)
S1Mtab = np.loadtxt("./Data/s1log_massdep.txt", skiprows=0)
S2Mtab = np.loadtxt("./Data/s2log_massdep.txt", skiprows=0)

S0MInt = RectBivariateSpline(M1tab, z2tab, S0Mtab) # 2-D Interpolation
S1MInt = RectBivariateSpline(M1tab, z2tab, S1Mtab) # 2-D Interpolation
S2MInt = RectBivariateSpline(M1tab, z2tab, S2Mtab) # 2-D Interpolation

def S0_M_fun(M, z):

    if M < 0.1:
        M = 0.1
    elif M > 100.:
        M = 100.

    return S0MInt(np.log(M), np.log(z))[0,0]

def S1_M_fun(M, z):

    if M < 0.1:
        M = 0.1
    elif M > 100.:
        M = 100.

    return S1MInt(np.log(M), np.log(z))[0,0]

def S2_M_fun(M, z):

    if M < 0.1:
        M = 0.1
    elif M > 100.:
        M = 100.

    return S2MInt(np.log(M), np.log(z))[0,0]

#-------------------------------------------------------------------------------------#
#                    Quantum Kinetic Equations for ARS Leptogenesis                   #
#-------------------------------------------------------------------------------------#

#
# NOTE: the mass labelling is off. What is called M1 here is called self.M2 in EtaB
#
def fast_RHS(z, y, Fmat, M1, deltaM, Tew, gss, M0, M_mat, Dm2_mat, chi_mat, Lvec, Rvec, acr, use_interpolation=True):
    if use_interpolation:

        S0 = S0_M_fun(M1, z) * (z*z/(Tew*Tew))
        S1 = S1_M_fun(M1, z) * (z*z/(Tew*Tew))
        S2 = S2_M_fun(M1, z) * (z*z/(Tew*Tew))
        G0 = G0_fun(z)
        G1 = G1_fun(z)
        G2 = G2_fun(z)
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
    mud_mat     =  np.diag([y[8], y[9], y[10]])
    mu_mat      = 2 * chi_mat @ mud_mat

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
    fdyneq = f_DYNeq(M1, z, Tew, gss)
    fyneq = f_YNeq(M1, z, Tew, gss)

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


    # This is just a 3 vector and we only care about diagonal elements of the 3x3 matrices
    eqtns = np.array([0+0j,0+0j, 0+0j, 0+0j, 0+0j, 0+0j, 0+0j, 0+0j, 0+0j, 0+0j, 0+0j], dtype=np.complex128)

    for i in range(3):
        eqtns[8+i] = cons_2 * (M0/Tew) * (
                - 0.5 * G0 * (F_RN_Fdagger[i][i] - Fstar_RNb_Ftrans[i][i])
                +       G1 * F_Fdagger[i][i]* mu_mat[i][i]
                - 0.5 * G2 * (F_RN_Fdagger[i][i] + Fstar_RNb_Ftrans[i][i]) * mu_mat[i][i]
                + 0.5 * S0 * (Fstar_M_RN_M_Ftrans[i][i] - F_M_RNb_M_Fdagger[i][i])
                +       S1 * F_M_M_Fdagger[i][i] * mu_mat[i][i]
                - 0.5 * S2 * (Fstar_M_RN_M_Ftrans[i][i] + F_M_RNb_M_Fdagger[i][i]) * mu_mat[i][i]
            )

    eqtns[0]  = RNRHS_mat[0][0]
    eqtns[1]  = RNRHS_mat[0][1]
    eqtns[2]  = RNRHS_mat[1][0]
    eqtns[3]  = RNRHS_mat[1][1]

    eqtns[4]  = RNbRHS_mat[0][0]
    eqtns[5]  = RNbRHS_mat[0][1]
    eqtns[6]  = RNbRHS_mat[1][0]
    eqtns[7]  = RNbRHS_mat[1][1]

    return eqtns

def fast_averaged_RHS(z, y, Fmat, M1, deltaM, Tew, gss, M0, M_mat,  use_interpolation=True):
    if use_interpolation:

        S0 = S0_M_fun(M1, z) * (z*z/(Tew*Tew))
        S1 = S1_M_fun(M1, z) * (z*z/(Tew*Tew))
        S2 = S2_M_fun(M1, z) * (z*z/(Tew*Tew))
        G0 = G0_fun(z)
        G1 = G1_fun(z)
        G2 = G2_fun(z)
    else:
        G0 = 0.01007
        G1 = 0.00547
        G2 = -0.00252
        S0 = 0.0402 * (z*z/(Tew*Tew))
        S1 = 0.00735 * (z*z/(Tew*Tew))
        S2 = -0.01616 * (z*z/(Tew*Tew))

    M_mat[0][0]   = M1
    M_mat[1][1]   = M1 + deltaM
    FdF     = Fmat.H @ Fmat
    deltaM2 = 2 * M1 * deltaM + deltaM * deltaM

    cons_1 = 0.05701803240084191
    cons_2 = 0.5480722270510788

    # RHS matrices

    RN_mat      =  np.diag([y[0], y[1]])
    RNb_mat     =  np.diag([y[2], y[3]])
    mud_mat     =  np.diag([y[4], y[5], y[6]])

    chi_mat     =  -1./711. * np.matrix([[257.,  20.,  20.], [20.,  257.,  20.], [20., 20., 257. ]], dtype=np.complex128)
    mu_mat      =  2. * chi_mat @ mud_mat


    # matrices appearing in Eqs

    Fdagger_mu_F        = Fmat.H @ mu_mat @ Fmat
    M_Ftrans_Fstar_M    = M_mat @ Fmat.T @ np.conjugate(Fmat) @ M_mat
    M_Ftrans_mu_Fstar_M = M_mat @ Fmat.T @ mu_mat @ np.conjugate(Fmat) @ M_mat

    Ftrans_mu_Fstar     = Fmat.T @ mu_mat @ np.conjugate(Fmat)
    M_Fdagger_F_M       = M_mat @ Fmat.H @ Fmat @ M_mat
    M_Fdagger_mu_F_M    = M_mat @ Fmat.H @ mu_mat @ Fmat @ M_mat

    F_RN_Fdagger        = Fmat @ RN_mat @ Fmat.H
    Fstar_RNb_Ftrans    = np.conjugate(Fmat) @ RNb_mat @ Fmat.T
    F_Fdagger           = Fmat @ Fmat.H
    Fstar_M_RN_M_Ftrans = np.conjugate(Fmat) @ M_mat @ RN_mat @ M_mat @ Fmat.T
    F_M_RNb_M_Fdagger   = Fmat @ M_mat @ RNb_mat @ M_mat @ Fmat.H
    F_M_M_Fdagger       = Fmat @ M_mat @ M_mat @ Fmat.H
    RNmatmId2           = RN_mat  - np.identity(2)
    RNbmatmId2          = RNb_mat - np.identity(2)



    # ARS Equations

    RNRHS_mat   = (M0/Tew) * (- 0.5 * G0 * anticommutator(FdF, RNmatmId2)
                              + G1 * Fdagger_mu_F
                              - 0.5 * G2 * anticommutator(Fdagger_mu_F, RN_mat)
                              - 0.5 * S0 * anticommutator(M_Ftrans_Fstar_M, RNmatmId2)
                              - S1 * M_Ftrans_mu_Fstar_M
                              + 0.5 * S2 * anticommutator(M_Ftrans_mu_Fstar_M, RN_mat)
                              - (Tew/M0) * (RN_mat/f_YNeq(M1, z, Tew, gss)) * f_DYNeq(M1, z, Tew, gss)
                              )

    RNbRHS_mat  = (M0/Tew) * (- 0.5 * G0 * anticommutator(np.conjugate(FdF), RNbmatmId2)
                              - G1 * Ftrans_mu_Fstar
                              + 0.5 * G2 * anticommutator(Ftrans_mu_Fstar, RNb_mat)
                              - 0.5 * S0 * anticommutator(M_Fdagger_F_M,  RNbmatmId2) + S1 * M_Fdagger_mu_F_M
                              - 0.5 * S2 * anticommutator(M_Fdagger_mu_F_M, RNb_mat)
                              - (Tew/M0) * (RNb_mat/f_YNeq(M1, z, Tew, gss)) * f_DYNeq(M1, z, Tew, gss)
                              )

    muDeltaRHS  = cons_2 * (M0/Tew) * (- 0.5 * G0 * (F_RN_Fdagger - Fstar_RNb_Ftrans).diagonal()
                                       + G1 * (diagdiag(F_Fdagger) @ np.diag(mu_mat) )
                                       - 0.5 * G2 * (diagdiag(F_RN_Fdagger + Fstar_RNb_Ftrans) @ diagdiag(mu_mat) ).diagonal()
                                       + 0.5 * S0 * (Fstar_M_RN_M_Ftrans - F_M_RNb_M_Fdagger).diagonal()
                                       + S1 * (diagdiag(F_M_M_Fdagger)  @ np.diag(mu_mat) )
                                       - 0.5 * S2 * (diagdiag(Fstar_M_RN_M_Ftrans + F_M_RNb_M_Fdagger) @ np.diag(mu_mat) ).diagonal()
                                        )

    eqtns = np.array([0+0j, 0+0j, 0+0j, 0+0j, 0+0j, 0+0j, 0+0j], dtype=np.complex128)

    eqtns[0] = RNRHS_mat[0,0]
    eqtns[1] = RNRHS_mat[1,1]

    eqtns[2] = RNbRHS_mat[0,0]
    eqtns[3] = RNbRHS_mat[1,1]

    eqtns[4] = muDeltaRHS[0,0]
    eqtns[5] = muDeltaRHS[0,1]
    eqtns[6] = muDeltaRHS[0,2]


    return eqtns


class EtaB_ARS(ulysses.ULSBase):
    """
    add description of where to find BEs
    """

    def shortname(self): return "BEARS"

    def flavourindices(self): return [1]

    def flavourlabels(self): return ["$NBL$"]

    def RHS(self, z, y, Fmat, M1, deltaM, Tew, gss, M0, M_mat, Dm2_mat, chi_mat, Lvec, Rvec, acr):

        return fast_RHS(z, y, Fmat, M1, deltaM, Tew, gss, M0, M_mat, Dm2_mat, chi_mat, Lvec, Rvec, acr, use_interpolation=False)

    def RHS_averaged(self, z, y, Fmat, M1, deltaM, Tew, gss, M0, M_mat):

        return fast_averaged_RHS(z, y, Fmat, M1, deltaM, Tew, gss, M0, M_mat, use_interpolation=False)

    @property
    def EtaB(self):

        # global constants (masses in GeV)
        Tew    = 131.7
        gss    = 106.75
        M0     = 7.e17
        zeta3  = zeta(3.)
        gstaro = 43/11 #entropic effective degrees of freedom at present

        mp       = 1.672621898*1e-24 #proton mass in g
        ngamma   = 410.7 #present photon number density in cm^-3
        rhoc     = 1.87840*1e-29 #critical density of the Universe in h^2 g cm^-3
        gstaro   = 43/11 #entropic effective degrees of freedom at present
        ToYb     = 45 * zeta3 /(gstaro * np.pi**4)
        ToOmegab = mp * ngamma/rhoc

        self.M1 = 0.
        self.m  = 0.

        self.x1 *= -1.
        self.y1 *= -1.

        self.delta *= -1.
        self.a21   *= -1.

        Fmat     = np.delete(self.h, 0, 1)

        if 30. <= self.M2 and self.M2 < 100.:
            print("\n")
            print(colored("Warning: RHN masses sufficiently heavy that non-relativistic corrections to lepton-number-violating rates can be important at the time of sphaleron decoupling.\n", 'red'))

        if self.M2 >= 100.:
            print("\n")
            print(colored("Warning: Our implementation of the kinetic equations assumes that RHNs are at least somewhat relativistic. \n\t Results are not validated and may be incorrect for this set of RHN masses. \n", 'red'))


        # FFd = Fmat @ Fmat.H
        FFd = Fmat @ np.transpose(np.conjugate( Fmat))

        G1_t = G1_fun(1.)
        # TODO inconsistent/confusing labelling --- also this recalculates self.M2 which is then used as input parameter to RHS (but called M1)
        S1_t = S1_M_fun(self.M2, 1.)


        G1_Hub = 0.5*(M0/Tew) * G1_t * np.diag(FFd).real
        G2_Hub = 0.5*(M0/Tew) * S1_t * np.diag(FFd).real

        if np.min(G1_Hub) >= 1. and  np.min(G2_Hub) >= 1.:
            print("\n")
            print(colored("Warning: This set of parameters is in the strong washout regime for lepton asymmetries. \n", 'red'))

        dMval = self.M3 - self.M2

        zcut = 0.3

        params  = np.array([Fmat[0, 0], Fmat[0, 1], Fmat[1, 0], Fmat[1, 1], Fmat[2, 0], Fmat[2, 1], self.M2, dMval, Tew, gss, M0],
                           dtype=np.complex128)

        # initial conditions in the order RN11, RN12, RN21, RN22, RNb11, RNb12, RNb21, RNb22, mudelta1, mudelta2,  mudelta3

        YNin = 0.+0j
        y0 = np.array([YNin, 0+0j, 0+0j, YNin, YNin, 0+0j, 0+0j, YNin, 0+0j, 0+0j, 0+0j], dtype=np.complex128)


        M_mat    = np.zeros([2,2])
        M_mat[0][0]   = self.M2
        M_mat[1][1]   = self.M2 + dMval
        deltaM2 = 2 * self.M2 * dMval + dMval * dMval
        Dm2_mat = np.array([[0 , 0], [0 , deltaM2 ]])

        chi_mat     =  -1./711. * np.array([[257.,  20.,  20.], [20.,  257.,  20.], [20., 20., 257. ]], dtype=np.complex128)
        Lvec = np.zeros((16,2), dtype=np.complex128)
        Rvec = np.zeros((16,2), dtype=np.complex128)
        acr = np.zeros((16,2), dtype=np.complex128)

        if dMval <= 1.e-8:

            ys = solve_ivp(lambda t, z: self.RHS(t, z, Fmat, self.M2, dMval, Tew, gss, M0, M_mat, Dm2_mat, chi_mat, Lvec, Rvec, acr), [1.e-6, 1], y0,
                           method='BDF', rtol=1.e-7, atol=1.e-10)

            muD1, muD2, muD3 = [ys.y[8,-1], ys.y[9,-1], ys.y[10,-1]]


        # TODO propagate changes from if branch and make corresponding changes for the RHS_averaged bit
        else:
            print(self.M2, dMval)
            print("stitching\n")

            ys_1 = solve_ivp(lambda t, z: self.RHS(t, z, Fmat, self.M2, dMval, Tew, gss, M0, M_mat, Dm2_mat, chi_mat), [1.e-6, zcut], y0,
                             method='BDF', rtol=1.e-7, atol=1.e-10)


            y0_2  = np.array([np.abs(ys_1.y[0,-1]),  np.abs(ys_1.y[3,-1]),  np.abs(ys_1.y[4,-1]), np.abs(ys_1.y[7,-1]),
                              np.real(ys_1.y[8,-1]), np.real(ys_1.y[9,-1]), np.real(ys_1.y[10,-1])], dtype=np.complex128)

            ys = solve_ivp(lambda t, z: self.RHS_averaged(t, z, Fmat, self.M2, dMval, Tew, gss, M0, M_mat),
                             [zcut, 1], y0_2, method='BDF', rtol=1.e-5, atol=1.e-10)

            muD1, muD2, muD3 = [ys.y[4,-1], ys.y[5,-1], ys.y[6,-1]]

        YB_sol  = np.abs(muD1 + muD2 + muD3)
        YB  = f_convertYBLtoYB(f_convertmutoY(YB_sol, gss))
        etaB= YB/ToYb

        return etaB

class EtaB_ARS_INTERP(EtaB_ARS):

    def shortname(self): return "BEARS_INTERP"

    def RHS(self, z, y, Fmat, M1, deltaM, Tew, gss, M0, M_mat, Dm2_mat, chi_mat, Lvec, Rvec, acr):

        return fast_RHS(z, y, Fmat, M1, deltaM, Tew, gss, M0, M_mat, Dm2_mat, chi_mat, Lvec, Rvec, acr, use_interpolation=True)

    def RHS_averaged(self, z, y, Fmat, M1, deltaM, Tew, gss, M0, M_mat):

        return fast_averaged_RHS(z, y, Fmat, M1, deltaM, Tew, gss, M0, M_mat, use_interpolation=True)
