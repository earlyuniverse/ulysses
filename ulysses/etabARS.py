# ARS leptogenesis
from scipy.special import kn
import ulysses
import numpy as np
from odeintw import odeintw
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy.interpolate import interp1d, RectBivariateSpline
from scipy.integrate import quad, ode, solve_ivp, odeint, LSODA
from scipy.special import zeta, kn
import math
plt.rcParams['text.usetex'] = True
from termcolor import colored
#import numba
#---------------------------------------------#
#              Auxiliary Functions            #
#---------------------------------------------#

def ticks_log(i0, i1):
    ticks = []
    minor = []
    for i in range(i0, i1):
        for j in range(0, 9):
            if j == 0:
                minor.append('major')
            else:
                minor.append('minor')
            ticks.append((1+j)*10**i)
    ticks2 = np.append(ticks,10**i1)
    minor.append('major')
    return ticks2.astype(float), minor

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
    return (z * 15.) /(2 * np.pi**2 * gss)

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
    Compute 8 anticommutators of 2x2 matrices stores in 16x2 arrays
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

        G0_M_fun, G1_M_fun, G2_M_fun, S0_M_fun, S1_M_fun, S2_M_fun = funcs

        G0 = G0_M_fun(M2,z)
        G1 = G1_M_fun(M2,z)
        G2 = G2_M_fun(M2,z)
        S0 = S0_M_fun(M2, z) * (z*z/(Tew*Tew))
        S1 = S1_M_fun(M2, z) * (z*z/(Tew*Tew))
        S2 = S2_M_fun(M2, z) * (z*z/(Tew*Tew))
    else:
        G0 = 9.15e-3#0.01007
        G1 = 5.1e-3#0.00547
        G2 = -2.19e-3#-0.00252
        S0 = 0.04337 * (z*z/(Tew*Tew))
        S1 = 0.00855 * (z*z/(Tew*Tew))
        S2 = -0.01651 * (z*z/(Tew*Tew))



    T = Tew/z

    cons_1 = 0.05701803240084191
    cons_2 = 0.5480722270510788

    # RHS matrices
    RN_mat      =  np.array([[y[0], y[1]], [y[2], y[3]]], dtype=np.complex128)
    RNb_mat     =  np.array([[y[4], y[5]], [y[6], y[7]]], dtype=np.complex128)
    mud_vec     =  np.array([y[8], y[9], y[10]], dtype=np.complex128)
    mud_mat     =  np.diag([y[8], y[9], y[10]])
    mu_mat_old  = 2. * chi_mat @ mud_mat#THIS IS WRONG!
    mu_vec      = 2. * chi_mat @ mud_vec
    mu_mat      = np.diag([mu_vec[0], mu_vec[1], mu_vec[2]])
    
    
    """
    #Kill the off-diagonal terms
    RN_mat[0][1] = 0
    RN_mat[1][0] = 0
    RNb_mat[0][1] = 0
    RNb_mat[1][0] = 0
    """

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
    
    #Additional factors for checks
    F_RNb_Fdagger        = Fmat @ RNb_mat @ FmatH
    Fstar_M_RN_M_Ftrans = np.conjugate(Fmat) @ M_mat @ RN_mat @ M_mat @ FmatT
    Fstar_M_RNb_M_Ftrans = np.conjugate(Fmat) @ M_mat @ RNb_mat @ M_mat @ FmatT
    F_M_RN_M_Fdagger   = Fmat @ M_mat @ RN_mat @ M_mat @ FmatH
    Fstar_M_M_FT       = FmatC @ M_mat @ M_mat @ FmatT


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
    #fdyneq = f_DYNeq(M2, z, Tew, gss)
    #fyneq = f_YNeq(M2, z, Tew, gss)
    temp_z = M2*z/Tew
    zz = temp_z.real
    kn1, kn2, kn3 = kn([1,2,3], zz)
    if kn2 == 0:
        fact = -(M2/Tew)  * zz/2
        fact2 = cons_2
    else:
        fact = -(M2/Tew)  * (kn1/kn2)
        fact2 = cons_2 *2/(zz**2 * kn2)
    
    
    
    # ARS Equations
    RNRHS_mat   = (M0/Tew) * (- 1j * commutator(Ham_RN, RN_mat) 
                              - 0.5 * G0 * acr[0:2]
                              + G1 * Fdagger_mu_F
                              - 0.5 * G2 * acr[2:4]
                              - 0.5 * S0 * acr[4:6]
                              - S1 * M_Ftrans_mu_Fstar_M
                              + 0.5 * S2 * acr[6:8]
                              + (Tew/M0) * RN_mat * fact
                              )
    

    RNbRHS_mat  = (M0/Tew) * (- 1j * commutator(Ham_RNb, RNb_mat)
                              - 0.5 * G0 * acr[8:10]
                              - G1 * Ftrans_mu_Fstar
                              + 0.5 * G2 * acr[10:12]
                              - 0.5 * S0 * acr[12:14]
                              + S1 * M_Fdagger_mu_F_M
                              - 0.5 * S2 * acr[14:16]
                              + (Tew/M0) * (RNb_mat) * fact
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
        eqtns[8+i] = fact2 * (M0/Tew) * (
                - 0.5 * G0 * (F_RN_Fdagger[i][i] - Fstar_RNb_Ftrans[i][i])
                +       G1 * F_Fdagger[i][i]* mu_vec[i]
                - 0.5 * G2 * (F_RN_Fdagger[i][i] + Fstar_RNb_Ftrans[i][i]) * mu_vec[i]
                + 0.5 * S0 * (Fstar_M_RN_M_Ftrans[i][i] - F_M_RNb_M_Fdagger[i][i])
                +       S1 * Fstar_M_M_FT[i][i] * mu_vec[i]
                -0.5 * S2 * (Fstar_M_RN_M_Ftrans[i][i] + F_M_RNb_M_Fdagger[i][i]) * mu_vec[i]
            )
        
    return eqtns

def real_RHS(z, y, Fmat, M2, deltaM, Tew, gss, M0, M_mat, Dm2_mat, chi_mat, Lvec, Rvec, acr, funcs, use_interpolation=True):

    if use_interpolation:

        G0_M_fun, G1_M_fun, G2_M_fun, S0_M_fun, S1_M_fun, S2_M_fun = funcs

        G0 = G0_M_fun(M2,z)
        G1 = G1_M_fun(M2,z)
        G2 = G2_M_fun(M2,z)
        S0 = S0_M_fun(M2, z) * (z*z/(Tew*Tew))
        S1 = S1_M_fun(M2, z) * (z*z/(Tew*Tew))
        S2 = S2_M_fun(M2, z) * (z*z/(Tew*Tew))
    else:
        G0 = 9.15e-3#0.01007
        G1 = 5.1e-3#0.00547
        G2 = -2.19e-3#-0.00252
        S0 = 0.04337 * (z*z/(Tew*Tew))
        S1 = 0.00855 * (z*z/(Tew*Tew))
        S2 = -0.01651 * (z*z/(Tew*Tew))



    T = Tew/z

    cons_1 = 0.05701803240084191
    cons_2 = 0.5480722270510788

    # RHS matrices
    RN_mat      =  np.array([[y[0]+y[3], y[1]-1j*y[2]], [y[1]+1j*y[2], y[0]-y[3]]], dtype=np.complex128)
    RNb_mat     =  np.array([[y[4]+y[7], y[5]-1j*y[6]], [y[5]+1j*y[6], y[4]-y[7]]], dtype=np.complex128)
    r0          =  y[0]
    r1          =  y[1]
    r2          =  y[2]
    r3          =  y[3]
    rb0         =  y[4]
    rb1         =  y[5]
    rb2         =  y[6]
    rb3         =  y[7]
    mud_vec     =  np.array([y[8], y[9], y[10]], dtype=np.float64)
    mud_mat     =  np.diag([y[8], y[9], y[10]])
    #mu_mat_old  = 2. * chi_mat @ mud_mat#THIS IS WRONG!
    mu_vec      = 2. * chi_mat @ mud_vec
    mu_mat      = np.diag([mu_vec[0], mu_vec[1], mu_vec[2]])
    
    #Pauli matrices
    s0          = np.array([[1, 0], [0, 1]], dtype=np.complex128)
    s1          = np.array([[0, 1], [1, 0]], dtype=np.complex128)
    s2          = np.array([[0, -1j], [1j, 0]], dtype=np.complex128)
    s3          = np.array([[1, 0], [0, -1]], dtype=np.complex128)
    
    
    """
    #Kill the off-diagonal terms
    RN_mat[0][1] = 0
    RN_mat[1][0] = 0
    RNb_mat[0][1] = 0
    RNb_mat[1][0] = 0
    """

    # matrices appearing in Eqs
    FmatH               = np.transpose(np.conjugate(Fmat))
    FmatT               = np.transpose(Fmat)
    FmatC               = np.conjugate(Fmat)
    FdF                 = np.transpose(np.conjugate(Fmat)) @ Fmat

    Ham_RN              = cons_1 * (FdF + 4.* (z*z/(Tew*Tew)) * Dm2_mat)
    Fdagger_mu_F        = FmatH @ mu_mat @ Fmat
    M_Ftrans_Fstar_M    = M_mat @ FmatT @ FmatC @ M_mat
    M_Ftrans_mu_Fstar_M = M_mat @ FmatT @ mu_mat @ FmatC @ M_mat
    Ham_s0              = Ham_RN@s0 - s0@Ham_RN
    Ham_s1              = Ham_RN@s1 - s1@Ham_RN
    Ham_s2              = Ham_RN@s2 - s2@Ham_RN
    Ham_s3              = Ham_RN@s3 - s3@Ham_RN
    FdF_s0              = FdF@s0    + s0@FdF
    FdF_s1              = FdF@s1    + s1@FdF
    FdF_s2              = FdF@s2    + s2@FdF
    FdF_s3              = FdF@s3    + s3@FdF

    Ham_RNb             = cons_1 * (np.conjugate(FdF) + 4.* (z*z/(Tew*Tew)) * Dm2_mat)
    Ftrans_mu_Fstar     = FmatT @ mu_mat @ FmatC
    M_Fdagger_F_M       = M_mat @ FmatH @ Fmat @ M_mat
    M_Fdagger_mu_F_M    = M_mat @ FmatH @ mu_mat @ Fmat @ M_mat
    Hamb_s0             = Ham_RNb@s0 - s0@Ham_RNb
    Hamb_s1             = Ham_RNb@s1 - s1@Ham_RNb
    Hamb_s2             = Ham_RNb@s2 - s2@Ham_RNb
    Hamb_s3             = Ham_RNb@s3 - s3@Ham_RNb
    MFTFsM_s0           = Ftrans_mu_Fstar@s0    + s0@Ftrans_mu_Fstar 
    MFTFsM_s1           = Ftrans_mu_Fstar@s1    + s1@Ftrans_mu_Fstar 
    MFTFsM_s2           = Ftrans_mu_Fstar@s2    + s2@Ftrans_mu_Fstar 
    MFTFsM_s3           = Ftrans_mu_Fstar@s3    + s3@Ftrans_mu_Fstar 
    

    F_RN_Fdagger        = Fmat @ RN_mat @ FmatH
    Fstar_RNb_Ftrans    = FmatC @ RNb_mat @ FmatT
    F_Fdagger           = Fmat @ FmatH
    Fstar_M_RN_M_Ftrans = FmatC @ M_mat @ RN_mat @ M_mat @ FmatT
    F_M_RNb_M_Fdagger   = Fmat @ M_mat @ RNb_mat @ M_mat @ FmatH
    F_M_M_Fdagger       = Fmat @ M_mat @ M_mat @ FmatH
    RNmatmId2           = RN_mat  - np.identity(2)
    RNbmatmId2          = RNb_mat - np.identity(2)
    
    #Additional factors for checks
    F_RNb_Fdagger        = Fmat @ RNb_mat @ FmatH
    Fstar_M_RN_M_Ftrans = np.conjugate(Fmat) @ M_mat @ RN_mat @ M_mat @ FmatT
    Fstar_M_RNb_M_Ftrans = np.conjugate(Fmat) @ M_mat @ RNb_mat @ M_mat @ FmatT
    F_M_RN_M_Fdagger   = Fmat @ M_mat @ RN_mat @ M_mat @ FmatH
    Fstar_M_M_FT       = FmatC @ M_mat @ M_mat @ FmatT


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
    #fdyneq = f_DYNeq(M2, z, Tew, gss)
    #fyneq = f_YNeq(M2, z, Tew, gss)
    temp_z = M2*z/Tew
    zz = temp_z.real
    kn1, kn2, kn3 = kn([1,2,3], zz)
    if kn2 == 0:
        fact = -(M2/Tew)  * zz/2
        fact2 = cons_2
    else:
        fact = -(M2/Tew)  * (kn1/kn2)
        fact2 = cons_2 *2/(zz**2 * kn2)
    
    

    
    # ARS Equations
    RNRHS_mat   = (M0/Tew) * (- 1j * commutator(Ham_RN, RN_mat) 
                              - 0.5 * G0 * acr[0:2]
                              + G1 * Fdagger_mu_F
                              - 0.5 * G2 * acr[2:4]
                              - 0.5 * S0 * acr[4:6]
                              - S1 * M_Ftrans_mu_Fstar_M
                              + 0.5 * S2 * acr[6:8]
                              + (Tew/M0) * RN_mat * fact
                              )
    

    RNbRHS_mat  = (M0/Tew) * (- 1j * commutator(Ham_RNb, RNb_mat)
                              - 0.5 * G0 * acr[8:10]
                              - G1 * Ftrans_mu_Fstar
                              + 0.5 * G2 * acr[10:12]
                              - 0.5 * S0 * acr[12:14]
                              + S1 * M_Fdagger_mu_F_M
                              - 0.5 * S2 * acr[14:16]
                              + (Tew/M0) * (RNb_mat) * fact
                              )
  
   
    #Define the real equations
    eqtns = np.zeros(11, dtype=np.float64)

    eqtns[0]  = np.real(RNRHS_mat[0][0]+RNRHS_mat[1][1])/2
    eqtns[1]  = np.real(RNRHS_mat[0][1]+RNRHS_mat[1][0])/2
    eqtns[2]  = np.real((-RNRHS_mat[0][1]+RNRHS_mat[1][0])/2j)
    eqtns[3]  = np.real(RNRHS_mat[0][0]-RNRHS_mat[1][1])/2

    eqtns[4]  = np.real(RNbRHS_mat[0][0]+RNbRHS_mat[1][1])/2
    eqtns[5]  = np.real(RNbRHS_mat[0][1]+RNbRHS_mat[1][0])/2
    eqtns[6]  = np.real((-RNbRHS_mat[0][1]+RNbRHS_mat[1][0])/2j)
    eqtns[7]  = np.real(RNbRHS_mat[0][0]-RNbRHS_mat[1][1])/2
    
    Fs0Fdag     = Fmat@s0@FmatH
    Fs1Fdag     = Fmat@s1@FmatH
    Fs2Fdag     = Fmat@s2@FmatH
    Fs3Fdag     = Fmat@s3@FmatH
    FstarMs0MFT = np.conjugate(Fmat)@M_mat@s0@M_mat@FmatT
    FstarMs1MFT = np.conjugate(Fmat)@M_mat@s1@M_mat@FmatT
    FstarMs2MFT = np.conjugate(Fmat)@M_mat@s2@M_mat@FmatT
    FstarMs3MFT = np.conjugate(Fmat)@M_mat@s3@M_mat@FmatT
    
    Fstars0FT   = np.conjugate(Fmat)@s0@FmatT
    Fstars1FT   = np.conjugate(Fmat)@s1@FmatT
    Fstars2FT   = np.conjugate(Fmat)@s2@FmatT
    Fstars3FT   = np.conjugate(Fmat)@s3@FmatT
    FMs0MFdag   = Fmat@M_mat@s0@M_mat@FmatH
    FMs1MFdag   = Fmat@M_mat@s1@M_mat@FmatH
    FMs2MFdag   = Fmat@M_mat@s2@M_mat@FmatH
    FMs3MFdag   = Fmat@M_mat@s3@M_mat@FmatH

    # This is just a 3 vector and we only care about diagonal elements of the 3x3 matrices
    for i in range(3):
        
        eqtns[8+i] = np.real(fact2 * (M0/Tew) * (
                - 0.5 * G0 * (Fs0Fdag[i][i]*r0 - Fstars0FT[i][i]*rb0)
                - 0.5 * G0 * (Fs1Fdag[i][i]*r1 - Fstars1FT[i][i]*rb1)
                - 0.5 * G0 * (Fs2Fdag[i][i]*r2 - Fstars2FT[i][i]*rb2)
                - 0.5 * G0 * (Fs3Fdag[i][i]*r3 - Fstars3FT[i][i]*rb3)
                +       G1 * F_Fdagger[i][i]* mu_vec[i]
                - 0.5 * G2 * (Fs0Fdag[i][i]*r0 + Fstars0FT[i][i]*rb0) * mu_vec[i]
                - 0.5 * G2 * (Fs1Fdag[i][i]*r1 + Fstars1FT[i][i]*rb1) * mu_vec[i]
                - 0.5 * G2 * (Fs2Fdag[i][i]*r2 + Fstars2FT[i][i]*rb2) * mu_vec[i]
                - 0.5 * G2 * (Fs3Fdag[i][i]*r3 + Fstars3FT[i][i]*rb3) * mu_vec[i]
                + 0.5 * S0 * (FstarMs0MFT[i][i]*r0 - FMs0MFdag[i][i]*rb0)
                + 0.5 * S0 * (FstarMs1MFT[i][i]*r1 - FMs1MFdag[i][i]*rb1)                
                + 0.5 * S0 * (FstarMs2MFT[i][i]*r2 - FMs2MFdag[i][i]*rb2)
                + 0.5 * S0 * (FstarMs3MFT[i][i]*r3 - FMs3MFdag[i][i]*rb3)
                +       S1 * Fstar_M_M_FT[i][i] * mu_vec[i]
                -0.5 * S2 * (FstarMs0MFT[i][i]*r0 + FMs0MFdag[i][i]*rb0) * mu_vec[i]
                -0.5 * S2 * (FstarMs1MFT[i][i]*r1 + FMs1MFdag[i][i]*rb1) * mu_vec[i]
                -0.5 * S2 * (FstarMs2MFT[i][i]*r2 + FMs2MFdag[i][i]*rb2) * mu_vec[i]
                -0.5 * S2 * (FstarMs3MFT[i][i]*r3 + FMs3MFdag[i][i]*rb3) * mu_vec[i]
            ))


    return eqtns

def fast_averaged_RHS(z, y, Fmat, M2, deltaM, Tew, gss, M0, M_mat, Dm2_mat, chi_mat, Lvec, Rvec, acr, funcs, use_interpolation=True):

    if use_interpolation:

        G0_M_fun, G1_M_fun, G2_M_fun, S0_M_fun, S1_M_fun, S2_M_fun = funcs

        G0 = G0_M_fun(M2,z)
        G1 = G1_M_fun(M2,z)
        G2 = G2_M_fun(M2,z)
        S0 = S0_M_fun(M2, z) * (z*z/(Tew*Tew))
        S1 = S1_M_fun(M2, z) * (z*z/(Tew*Tew))
        S2 = S2_M_fun(M2, z) * (z*z/(Tew*Tew))
    else:
        G0 = 9.15e-3#0.01007
        G1 = 5.1e-3#0.00547
        G2 = -2.19e-3#-0.00252
        S0 = 0.04337 * (z*z/(Tew*Tew))
        S1 = 0.00855 * (z*z/(Tew*Tew))
        S2 = -0.01651 * (z*z/(Tew*Tew))

    M_mat[0][0]   = M2
    M_mat[1][1]   = M2 + deltaM

    deltaM2 = 2 * M2 * deltaM + deltaM * deltaM

    cons_1 = 0.05701803240084191
    cons_2 = 0.5480722270510788

    RN_mat      =  np.diag([y[0], y[1]])
    RNb_mat     =  np.diag([y[2], y[3]])
    mud_mat     =  np.diag([y[4], y[5], y[6]])
    mud_vec     =  np.array([y[4], y[5], y[6]])

    mu_mat_old  = 2. * chi_mat @ mud_mat#THIS IS WRONG!
    mu_vec      = 2. * chi_mat @ mud_vec
    mu_mat      = np.diag([mu_vec[0], mu_vec[1], mu_vec[2]])

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
    
    Fstar_M_RNb_M_Ftrans = np.conjugate(Fmat) @ M_mat @ RNb_mat @ M_mat @ FmatT
    F_M_RN_M_Fdagger   = Fmat @ M_mat @ RN_mat @ M_mat @ FmatH
    Fstar_M_M_FT       = FmatC @ M_mat @ M_mat @ FmatT


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
    #fdyneq = f_DYNeq(M2, z, Tew, gss)
    #fyneq = f_YNeq(M2, z, Tew, gss)
    temp_z = M2*z/Tew
    zz = temp_z.real
    kn1, kn2, kn3 = kn([1,2,3], zz)
    if kn2 == 0:
        fact = -(M2/Tew)  * zz/2
        fact2 = cons_2
    else:
        fact = -(M2/Tew)  * (kn1/kn2)
        fact2 = cons_2 *2/(zz**2 * kn2)

    # ARS Equations
    RNRHS_mat   = (M0/Tew) * (- 0.5 * G0 * acr[0:2]
                              + G1 * Fdagger_mu_F
                              - 0.5 * G2 * acr[2:4]
                              - 0.5 * S0 * acr[4:6]
                              - S1 * M_Ftrans_mu_Fstar_M
                              + 0.5 * S2 * acr[6:8]
                              - (Tew/M0) * (RN_mat) * fact
                              )
    
    RNbRHS_mat  = (M0/Tew) * (- 0.5 * G0 * acr[8:10]
                              - G1 * Ftrans_mu_Fstar
                              + 0.5 * G2 * acr[10:12]
                              - 0.5 * S0 * acr[12:14]
                              + S1 * M_Fdagger_mu_F_M
                              - 0.5 * S2 * acr[14:16]
                              - (Tew/M0) * (RNb_mat) * fact #/ fyneq) * fdyneq
                              )

    eqtns = np.zeros(7, dtype=np.complex128)

    eqtns[0] = RNRHS_mat[0,0]
    eqtns[1] = RNRHS_mat[1,1]

    eqtns[2] = RNbRHS_mat[0,0]
    eqtns[3] = RNbRHS_mat[1,1]

    for i in range(3):
        eqtns[4+i] = fact2 * (M0/Tew) * (
                - 0.5 * G0 * (F_RN_Fdagger[i][i] - Fstar_RNb_Ftrans[i][i])
                +       G1 * F_Fdagger[i][i]* mu_vec[i]
                - 0.5 * G2 * (F_RN_Fdagger[i][i] + Fstar_RNb_Ftrans[i][i]) * mu_vec[i]
                + 0.5 * S0 * (Fstar_M_RN_M_Ftrans[i][i] - F_M_RNb_M_Fdagger[i][i])
                +       S1 * Fstar_M_M_FT[i][i] * mu_vec[i]
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
    
    def Re_RHS(self, z, y, Fmat, M2, deltaM, Tew, gss, M0, M_mat, Dm2_mat, chi_mat, Lvec, Rvec, acr):

        funcs = []
        return real_RHS(z, y, Fmat, M2, deltaM, Tew, gss, M0, M_mat, Dm2_mat, chi_mat, Lvec, Rvec, acr, funcs, use_interpolation=False)


    def RHS_averaged(self, z, y, Fmat, M2, deltaM, Tew, gss, M0, M_mat, Dm2_mat, chi_mat, Lvec, Rvec, acr):

        funcs = []
        return fast_averaged_RHS(z, y, Fmat, M2, deltaM, Tew, gss, M0, M_mat, Dm2_mat, chi_mat, Lvec, Rvec, acr, funcs, use_interpolation=False)

    @property
    def EtaB(self):
        
        # global constants (masses in GeV)
        Tew    = 131.7
        gss    = 106.75
        M0     = 7.112582895088419e+17
        zeta3  = zeta(3.)

        mp       = 1.672621898*1e-24  #proton mass in g
        ngamma   = 410.7              #present photon number density in cm^-3
        rhoc     = 1.87840*1e-29      #critical density of the Universe in h^2 g cm^-3
        gstaro   = 43/11              #entropic effective degrees of freedom at present
        ToYb     = 45 * zeta3 /(gstaro * np.pi**4)
        ToOmegab = mp * ngamma/rhoc
        self.plot = False
        self.save_plot = False

        self.m  = 0.

        #Selects only the first two columns of the Yukawa matrix corresponding to N1 and N2. N3 is decoupled.
        Fmat     =  np.delete(self.h, 2, 1)
        
        if 30. <= self.M1 and self.M1 < 100.:
            print("\n")
            print(colored("Warning: RHN masses sufficiently heavy that non-relativistic corrections to lepton-number-violating rates can be important at the time of sphaleron decoupling.\n", 'red'))

        if self.M1 >= 100.:
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

        dMval = self.M2 - self.M1
    
        params  = np.array([Fmat[0, 0], Fmat[0, 1], Fmat[1, 0], Fmat[1, 1], Fmat[2, 0], Fmat[2, 1], self.M1, dMval, Tew, gss, M0],
                           dtype=np.complex128)

        # initial conditions in the order RN11, RN12, RN21, RN22, RNb11, RNb12, RNb21, RNb22, mudelta1, mudelta2,  mudelta3
        y0 = np.array([np.real(self.initial_abundance), 0., 0., 0., np.real(self.initial_abundance), 0., 0., 0., 0., 0., 0.], dtype=np.float64)

        M_mat       = np.zeros([2,2])
        M_mat[0][0] = self.M1
        M_mat[1][1] = self.M1 + dMval
        deltaM2     = 2 * self.M1 * dMval + dMval * dMval
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

            ys = solve_ivp(lambda t, z: self.Re_RHS(t, z, Fmat, self.M1, dMval, Tew, gss, M0, M_mat, Dm2_mat, chi_mat, Lvec, Rvec, acr),
                           [1e-6, 1], y0, method='BDF', atol=1e-13, rtol = 1e-13)
            
            t, rN, rNb, muD1, muD2, muD3 = [ys.t, (ys.y[0]+ys.y[3]), (ys.y[4]+ ys.y[7]), ys.y[8], ys.y[9], ys.y[10]]

        
        # if zosc < 0.1, oscillations start early. The numerical integration could be slow. The default behaviour is that the full integration is performed
        # unless the user provides a zcut  value (this is the point i the z-integration where the stitching occurs and the off diagonal terms are ignored)
        # see below for more details
        else:

            print(colored("Info: zosc = {}, stitching occurs at {}.".format(zosc,self._zcut), 'red'))
            
            if self._zcut == 1:
            
                    print(colored("Info: zcut is at default value 1, consider inputting  0.5 < zcut < 1 if integration is slow.".format(zosc,self._zcut), 'red'))

            ys_1 = solve_ivp(lambda t, z: self.Re_RHS(t, z, Fmat, self.M1, dMval, Tew, gss, M0, M_mat, Dm2_mat, chi_mat, Lvec, Rvec, acr),
                                 [1e-6, self._zcut], y0, method='BDF', atol=1e-9, rtol = 1e-10)
            
            t_1, rN_1, rNb_1, muD1_1, muD2_1, muD3_1 = [ys_1.t, (ys_1.y[0]+ys_1.y[3]), (ys_1.y[4]+ys_1.y[7]), ys_1.y[8], ys_1.y[9], ys_1.y[10]]

                #  this is where the stitch occurs if the user provides a zcut
                # For the stitching case, we ignore the evolution of the off-diagonal terms in the RN and RNb matrices.
                # Thus, we only take the solutions related to RN11, RN22, RNb11, RNb22 for the averaged equations

            y0_2  = np.array([np.abs(ys_1.y[0,-1]+ys_1.y[3,-1]),  np.abs(ys_1.y[0,-1]- ys_1.y[3,-1]),  np.abs(ys_1.y[4,-1]+ys_1.y[7,-1]), np.abs(ys_1.y[4,-1]-ys_1.y[7,-1]),
                              np.real(ys_1.y[8,-1]), np.real(ys_1.y[9,-1]), np.real(ys_1.y[10,-1])], dtype=np.complex128)

                

            solARS = solve_ivp(lambda t, z: self.RHS_averaged(t, z, Fmat, self.M1, dMval, Tew, gss, M0, M_mat, Dm2_mat, chi_mat, Lvec, Rvec, acr), [self._zcut, 1], y0_2, method='BDF', rtol=1.e-10, atol=1.e-10)

            t_2, rN_2, rNb_2, muD1_2, muD2_2, muD3_2 = [solARS.t, solARS.y[0], solARS.y[2], solARS.y[4], solARS.y[5], solARS.y[6]]
            
            t, muD1, muD2, muD3 = [np.concatenate((t_1,t_2), axis=0), np.concatenate((muD1_1, muD1_2) , axis=0), np.concatenate((muD2_1, muD2_2), axis=0), np.concatenate((muD3_1, muD3_2), axis=0)]
            rN, rNb = [np.concatenate((rN_1, rN_2), axis = 0), np.concatenate((rNb_1, rNb_2), axis = 0)  ]
        
        YB_sol_e = np.real(muD1[-1])
        YB_sol_mu = np.real(muD2[-1])
        YB_sol_tau = np.real(muD3[-1])
        YB_sol  = np.real(muD1[-1] + muD2[-1] + muD3[-1])
        #print(np.real(muD1[-1]),  np.real(muD2[-1]), np.real(muD3[-1]))
        
        
        cf = (28/79) * np.pi**2 /(27*6* zeta(3)) 
        YBe = f_convertYBLtoYB(f_convertmutoY(YB_sol_e, gss))
        YBmu = f_convertYBLtoYB(f_convertmutoY(YB_sol_mu, gss))
        YBtau = f_convertYBLtoYB(f_convertmutoY(YB_sol_tau, gss))
        YB  = f_convertYBLtoYB(f_convertmutoY(YB_sol, gss))
        
        etaBe = YBe/ToYb
        etaBmu = YBmu/ToYb
        etaBtau = YBtau/ToYb
        etaB= YB/ToYb
                
        Neq, Neq2 = [], []
        for x in t:
            temp = self.M1*x/Tew
            Neq.append((3/(16)) * temp**2 * kn(2,temp))
        
        lw = 0.9
        x1, x2 = -6, 0
        y1, y2 = -13, 1
        xticks, xminor = ticks_log(x1, x2)
        yticks, yminor = ticks_log(y1, y2)
        prop_cycle = plt.rcParams['axes.prop_cycle']
        colors = prop_cycle.by_key()['color']

        if self.plot == True:
            fig, ax = plt.subplots()
            ax.axhline(6e-10, color = 'grey', linewidth = 0.6)
            #ax.axvline(zosc, color = 'grey', linewidth = 0.6)
            if self._zcut < 1:
                ax.axvline(self._zcut, color = 'grey', linewidth = 0.6, linestyle = '-.')
            ax.plot(t, np.real(rN*Neq), label=r"$N_{N_1}$", linewidth = lw, color = colors[4])
            ax.plot(t, np.real(muD1), label=r"$\mu_{\Delta_e}$", linewidth = lw, color = colors[0])
            ax.plot(t, np.real(muD2), label=r"$\mu_{\Delta_\mu}$", linewidth = lw, color = colors[1])
            ax.plot(t, np.real(muD3), label=r"$\mu_{\Delta_\tau}$", linewidth = lw, color = colors[2])
            ax.plot(t, cf*np.real(muD1 + muD2 + muD3), label=r"$\eta_{B}$", linewidth = lw, color = colors[3])
            #ax.plot(t, np.real((rN-rNb)*Neq), label=r"$N_N-N_{\bar{N}}$", linewidth = lw)
            ax.plot(t, Neq, linewidth = lw, color = 'k', linestyle = ':')
            ax.plot(t, -np.real(muD1), linestyle = '--', color = colors[0], linewidth = lw)
            ax.plot(t, -np.real(muD2), linestyle = '--', color = colors[1], linewidth = lw)
            ax.plot(t, -np.real(muD3), linestyle = '--', color = colors[2], linewidth = lw)
            ax.plot(t, -cf*np.real(muD1 + muD2 + muD3), linestyle = '--', linewidth = lw, color = colors[3])
            #ax.plot(t, -np.real((rN-rNb)*Neq), linestyle = '--', linewidth = lw, color = colors[4])
            #ax.plot(t, np.real((rN-1)*Neq), linewidth = lw, linestyle = '--', color = colors[5])
            ax.set_xscale('log')
            ax.set_yscale('log')
            ax.set_xticks(xticks, minor = xminor)
            ax.set_yticks(yticks, minor = yminor)
            ax.set_xlim([10**x1,10**x2])
            ax.set_ylim([10**y1,10**y2])
            ax.axvline(1, lw = .8, ls = ':', color = 'k')
            ax.xaxis.set_ticks_position('both')
            ax.yaxis.set_ticks_position('both')
            ax.tick_params(direction='in',  which = 'both', labelsize = 13, width = 0.6)
            ax.set_xlabel(r"$T_{\rm{ew}}/T$", fontsize=13)#x=T_{ew}/T
            ax.set_ylabel(r"$|\eta_B|, \, |\mu_{\Delta_\alpha}|, \, N_{N_1}$", fontsize=13)
            ax.legend(loc = 'lower left', shadow = True, ncol = 1, prop = {'size' : 12})

            ax.set_title(r'\texttt{etabARS.py}', fontsize=13, loc = 'right')

            ax.text(1.5e-6, 8.5e-1, r'$N_{N}^{\rm{eq}}$', fontsize=11, color='k')

            ax.grid(which='both', linestyle='-', linewidth='0.6', color='grey', alpha=0.2)

            if self.save_plot == True:
                fig.savefig(self.path, dpi=300, bbox_inches='tight')
                print(f"Plot saved to {self.path}")

            plt.show()
                
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

        import os
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
        G0M_f = os.path.join(data_dir, "./data/G0log_massdep.txt")
        G1M_f = os.path.join(data_dir, "./data/G1log_massdep.txt")
        G2M_f = os.path.join(data_dir, "./data/G2log_massdep.txt")
        S0M_f = os.path.join(data_dir, "./data/S0log_massdep.txt")
        S1M_f = os.path.join(data_dir, "./data/S1log_massdep.txt")
        S2M_f = os.path.join(data_dir, "./data/S2log_massdep.txt")
        
        #S0M_f = os.path.join(data_dir, "./data/s0log_massdep.txt")
        #S1M_f = os.path.join(data_dir, "./data/s1log_massdep.txt")
        #S2M_f = os.path.join(data_dir, "./data/s2log_massdep.txt")

        M2tab  = np.loadtxt(M2_f, skiprows=0)
        z2tab  = np.loadtxt(z2_f, skiprows=0)
        G0Mtab = np.loadtxt(G0M_f, skiprows=0)
        G1Mtab = np.loadtxt(G1M_f, skiprows=0)
        G2Mtab = np.loadtxt(G2M_f, skiprows=0)
        S0Mtab = np.loadtxt(S0M_f, skiprows=0)
        S1Mtab = np.loadtxt(S1M_f, skiprows=0)
        S2Mtab = np.loadtxt(S2M_f, skiprows=0)
        
        self.G0MInt_ = RectBivariateSpline(M2tab, z2tab, G0Mtab) # 2-D Interpolation
        self.G1MInt_ = RectBivariateSpline(M2tab, z2tab, G1Mtab) # 2-D Interpolation
        self.G2MInt_ = RectBivariateSpline(M2tab, z2tab, G2Mtab) # 2-D Interpolation

        self.S0MInt_ = RectBivariateSpline(M2tab, z2tab, S0Mtab) # 2-D Interpolation
        self.S1MInt_ = RectBivariateSpline(M2tab, z2tab, S1Mtab) # 2-D Interpolation
        self.S2MInt_ = RectBivariateSpline(M2tab, z2tab, S2Mtab) # 2-D Interpolation

    
    def G0_fun(self,z): return interpolate.splev(math.log(z), self.G0Int_, der=0)

    def G1_fun(self,z): return interpolate.splev(math.log(z), self.G1Int_, der=0)

    def G2_fun(self,z): return interpolate.splev(math.log(z), self.G2Int_, der=0)

    

    def G0_M_fun(self, M, z):

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

    def RHS(self, z, y, Fmat, M2, deltaM, Tew, gss, M0, M_mat, Dm2_mat, chi_mat, Lvec, Rvec, acr):

        funcs = [self.G0_M_fun, self.G1_M_fun, self.G2_M_fun, self.S0_M_fun, self.S1_M_fun, self.S2_M_fun]

        return fast_RHS(z, y, Fmat, M2, deltaM, Tew, gss, M0, M_mat, Dm2_mat, chi_mat, Lvec, Rvec, acr, funcs, use_interpolation=True)

    def Re_RHS(self, z, y, Fmat, M2, deltaM, Tew, gss, M0, M_mat, Dm2_mat, chi_mat, Lvec, Rvec, acr):

        funcs = [self.G0_M_fun, self.G1_M_fun, self.G2_M_fun, self.S0_M_fun, self.S1_M_fun, self.S2_M_fun]

        return real_RHS(z, y, Fmat, M2, deltaM, Tew, gss, M0, M_mat, Dm2_mat, chi_mat, Lvec, Rvec, acr, funcs, use_interpolation=True)


    def RHS_averaged(self, z, y, Fmat, M2, deltaM, Tew, gss, M0, M_mat, Dm2_mat, chi_mat, Lvec, Rvec, acr):

        funcs = [self.G0_M_fun, self.G1_M_fun, self.G2_M_fun, self.S0_M_fun, self.S1_M_fun, self.S2_M_fun]

        return fast_averaged_RHS(z, y, Fmat, M2, deltaM, Tew, gss, M0, M_mat, Dm2_mat, chi_mat, Lvec, Rvec, acr, funcs, use_interpolation=True)
