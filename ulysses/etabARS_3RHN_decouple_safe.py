# ARS leptogenesis
from scipy.special import kn
import ulysses
import numpy as np
from odeintw import odeintw
from numba import jit, njit
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy.interpolate import interp1d, RectBivariateSpline
from scipy.integrate import quad, ode, solve_ivp, odeint, LSODA
from scipy.special import zeta, kn
import math
plt.rcParams['text.usetex'] = True
from termcolor import colored
import numba
from ulysses import ProdRates
from ulysses.ProdRates import interp_one_over_y0, interp_hp, interp_hm



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

def f_TSM(x, Tew):
    return Tew/x

def f_ss(x, Tew, gss):
    return (2 * np.pi  * np.pi * gss * f_TSM(x, Tew)**3)/ 45.

def f_HH(x, Tew, M0):
    return (f_TSM(x, Tew) * f_TSM(x, Tew))/M0

def f_nphieqSM(x, Tew):
    return f_TSM(x, Tew)**3/(np.pi * np.pi)

def f_YHeqSM(x, Tew, gss):
    return (2 * f_nphieqSM(x, Tew) )/f_ss(x, Tew, gss)

def f_nNeq(M, x, Tew):
    temp = M * x /Tew
    return (M * M * Tew * kn(2, temp.real)) / (2. * np.pi * np.pi * x )

def f_YNeq(M, x, Tew, gss):
    return f_nNeq(M, x, Tew)/ f_ss(x, Tew, gss)

def f_Yieldeq(M, x, Tew, gss):
    temp = M * x /Tew
    return (45. / (4. * np.pi**4 * gss) * (M * x / Tew) * (M * x / Tew)) * kn(2, temp.real)

def f_convertmutoY(x, gss):
    return (x * 15.) /(2 * np.pi**2 * gss)

def f_convertYBLtoYB(x):
    return x * 28./79.

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




@njit
def explicit_anticommutator(X, Y, R):
    """
    Compute anti commutator of 3x3 matrices X and Y, store result in 3x3 matrix R
    """
    for i in range(3):
        for j in range(3):
            R[i][j] = X[i][0]*Y[0][j] + X[i][1]*Y[1][j] + X[i][2]*Y[2][j] + Y[i][0]*X[0][j] + Y[i][1]*X[1][j] + Y[i][2]*X[2][j]

@njit
def explicit_anticommutator_array(L, R, r):
    """
    Compute 8 anticommutators of 3x3 matrices stored in 24x3 arrays
    """
    for i in range(8):
        explicit_anticommutator(L[3*i:3*(i+1)], R[3*i:3*(i+1)], r[3*i:3*(i+1)])

@njit
def diagdiag(mat):
    arr        = np.identity(3, dtype=np.complex128)
    arr[0][0]  = mat.item(0,0)
    arr[1][1]  = mat.item(1,1)
    arr[2][2]  = mat.item(2,2)
    return  arr

#The Gel-Mann matrices and the 3x3 identity
# @njit
def gellmann_matrices():
    lambda0 = np.array([[1, 0, 0],
                        [0, 1, 0],
                        [0, 0, 1]], dtype=np.complex128)
    
    lambda1 = np.array([[0, 1, 0],
                        [1, 0, 0],
                        [0, 0, 0]], dtype=np.complex128)

    lambda2 = np.array([[0, -1j, 0],
                        [1j, 0, 0],
                        [0, 0, 0]], dtype=np.complex128)

    lambda3 = np.array([[1, 0, 0],
                        [0, -1, 0],
                        [0, 0, 0]], dtype=np.complex128)

    lambda4 = np.array([[0, 0, 1],
                        [0, 0, 0],
                        [1, 0, 0]], dtype=np.complex128)

    lambda5 = np.array([[0, 0, -1j],
                        [0, 0, 0],
                        [1j, 0, 0]], dtype=np.complex128)

    lambda6 = np.array([[0, 0, 0],
                        [0, 0, 1],
                        [0, 1, 0]], dtype=np.complex128)

    lambda7 = np.array([[0, 0, 0],
                        [0, 0, -1j],
                        [0, 1j, 0]], dtype=np.complex128)

    lambda8 = np.array([[1/np.sqrt(3), 0, 0],
                        [0, 1/np.sqrt(3), 0],
                        [0, 0, -2/np.sqrt(3)]], dtype=np.complex128)

    return [lambda0, lambda1, lambda2, lambda3, lambda4, lambda5, lambda6, lambda7, lambda8]


#this function takes an hermitian matrix and give the coeffients of the decomposition in terms of the Gell-Mann matrices
@njit
def GellMann_coefficients(M):
    c0 = (1/3) * (M[0][0] + M[1][1] + M[2][2])
    c1 = (1/2) * (M[0][1] + M[1][0])
    c2 = (1j/2) * (M[0][1] - M[1][0])
    c3 = (1/2) * (M[0][0] - M[1][1])
    c4 = (1/2) * (M[0][2] + M[2][0])
    c5 = (1j/2) * (M[0][2] - M[2][0])
    c6 = (1/2) * (M[1][2] + M[2][1])
    c7 = (1j/2) * (M[1][2] - M[2][1])
    c8 = (np.sqrt(3)/6) * (M[0][0] + M[1][1] - 2*M[2][2])
    return [c0, c1, c2, c3, c4, c5, c6, c7, c8]


#use this if n3 is decoupled
def GellMann_coefficients_decoupled(M):
    c0 = (1/3) * (M[0][0] + M[1][1])   # + M[2][2])
    c1 = (1/2) * (M[0][1] + M[1][0])
    c2 = (1j/2) * (M[0][1] - M[1][0])
    c3 = (1/2) * (M[0][0] - M[1][1])
    c4 = 0   #(1/2) * (M[0][2] + M[2][0])
    c5 = 0   #(1j/2) * (M[0][2] - M[2][0])
    c6 = 0   #(1/2) * (M[1][2] + M[2][1])
    c7 = 0   #(1j/2) * (M[1][2] - M[2][1])
    c8 = (np.sqrt(3)/6) * (M[0][0] + M[1][1] )   #- 2*M[2][2])
    return [c0, c1, c2, c3, c4, c5, c6, c7, c8]


"""
ATTEMPT: Define all the quantities that do not depend on x outside of the RHS. This makes the code faster.
"""

#import time

#This function computes all the relevant combination that enters the coefficients of the 
#density matrix equations given the matrix of the Yukawas, the matrix of the RHN masses and that of the spectator processes, 
def compute_linearised_coefficients(Fmat, M_mat, chi_mat, Dm2_mat, n3_decoupled=False):
    #start_time = time.time()

    FmatH   = np.transpose(np.conjugate(Fmat))
    FmatT   = np.transpose(Fmat)
    FmatC   = np.conjugate(Fmat)
    FdF     = FmatH @ Fmat

    M_Ftrans_Fstar_M = M_mat @ FmatT @ FmatC @ M_mat
    M_Fdagger_F_M = M_mat @ FmatH @ Fmat @ M_mat

    L = gellmann_matrices()
    
    F_L_Fdag = [Fmat@l@FmatH for l in L]
    FstarM_L_MFtrans = [FmatC@M_mat@l@M_mat@FmatT for l in L]
    Fstar_L_Ftrans = [FmatC@l@FmatT for l in L]
    FM_L_MFdag = [Fmat@M_mat@l@M_mat@FmatH for l in L]

    chi_mat_diag_e = np.diag([chi_mat[0,0], chi_mat[1,0], chi_mat[2,0]])
    chi_mat_diag_mu = np.diag([chi_mat[0,1], chi_mat[1,1], chi_mat[2,1]])
    chi_mat_diag_tau = np.diag([chi_mat[0,2], chi_mat[1,2], chi_mat[2,2]])
    
    Fd_chie_F = FmatH@chi_mat_diag_e@Fmat
    Fd_chimu_F = FmatH@chi_mat_diag_mu@Fmat
    Fd_chitau_F = FmatH@chi_mat_diag_tau@Fmat
    MF_chie_FstarM = M_mat@FmatT@chi_mat_diag_e@FmatC@M_mat
    MF_chimu_FstarM = M_mat@FmatT@chi_mat_diag_mu@FmatC@M_mat
    MF_chitau_FstarM = M_mat@FmatT@chi_mat_diag_tau@FmatC@M_mat

    FFd = Fmat@FmatH 
    FMMFd = Fmat@M_mat@M_mat@FmatH

    #Terms in the equations for RN and RbN that are linear in RN-1 and RbN-1. Will need to multily by G0 and S0, and (M0/Tew)
    R_G0 = [- 0.5 * anticommutator(FdF, l) for l in L]
    R_S0 = [- 0.5 * anticommutator(M_Ftrans_Fstar_M, l) for l in L]

    Rb_G0 = [- 0.5 * anticommutator(np.conjugate(FdF), l) for l in L]
    Rb_S0 = [- 0.5 * anticommutator(M_Fdagger_F_M, l) for l in L]

    #Terms in the equations for RN and RbN that are linear in the chemical potentials, need to multiply by G1 and S1, and (M0/Tew)
    RMe_G1 = -2*Fd_chie_F
    RMe_S1 = 2*MF_chie_FstarM

    RMmu_G1 = -2*Fd_chimu_F
    RMmu_S1 = 2*MF_chimu_FstarM

    RMtau_G1 = -2*Fd_chitau_F
    RMtau_S1 = 2*MF_chitau_FstarM

    RbMe_G1 = 2*np.conjugate(Fd_chie_F)
    RbMe_S1 = -2*np.conjugate(MF_chie_FstarM)

    RbMmu_G1 = 2*np.conjugate(Fd_chimu_F)
    RbMmu_S1 = -2*np.conjugate(MF_chimu_FstarM)

    RbMtau_G1 = 2*np.conjugate(Fd_chitau_F)
    RbMtau_S1 = -2*np.conjugate(MF_chitau_FstarM)

    #Terms in the equations for the chemical potentials that are linear in RN-1 and RbN-1. Will need to multily by G0 and S0, and (M0/Tew) * fact2.
    MRe_G0 = [-0.5 * f1[0][0] for f1 in F_L_Fdag]
    MRe_S0 = [ 0.5 * f2[0][0] for f2 in FstarM_L_MFtrans]

    MRmu_G0 = [-0.5 * f1[1][1] for f1 in F_L_Fdag]
    MRmu_S0 = [ 0.5 * f2[1][1] for f2 in FstarM_L_MFtrans]

    MRtau_G0 = [-0.5 * f1[2][2] for f1 in F_L_Fdag]
    MRtau_S0 = [ 0.5 * f2[2][2] for f2 in FstarM_L_MFtrans]

    MRbe_G0 = [ 0.5 * f1[0][0] for f1 in Fstar_L_Ftrans]
    MRbe_S0 = [-0.5 * f2[0][0] for f2 in FM_L_MFdag]

    MRbmu_G0 = [0.5 * f1[1][1] for f1 in Fstar_L_Ftrans]
    MRbmu_S0 = [-0.5 * f2[1][1] for f2 in FM_L_MFdag]

    MRbtau_G0 = [0.5 * f1[2][2] for f1 in Fstar_L_Ftrans]
    MRbtau_S0 = [-0.5 * f2[2][2] for f2 in FM_L_MFdag]


    #Terms in the equations for the chemical potentials that are linear in the chemical potentials. Will need to multily by G1 and S1, and (M0/Tew) * fact2.
    Mee_G1 = -2 * FFd[0][0] * chi_mat[0][0] 
    Mee_S1 = -2 * FMMFd[0][0] * chi_mat[0][0]

    Mem_G1 = -2 * FFd[0][0] * chi_mat[0][1]
    Mem_S1 = -2 * FMMFd[0][0] * chi_mat[0][1]

    Met_G1 = -2 * FFd[0][0] * chi_mat[0][2]
    Met_S1 = -2 * FMMFd[0][0] * chi_mat[0][2]

    Mme_G1 = -2 * FFd[1][1] * chi_mat[1][0]
    Mme_S1 = -2 * FMMFd[1][1] * chi_mat[1][0]

    Mmm_G1 = -2 * FFd[1][1] * chi_mat[1][1]
    Mmm_S1 = -2 * FMMFd[1][1] * chi_mat[1][1]

    Mmt_G1 = -2 * FFd[1][1] * chi_mat[1][2]
    Mmt_S1 = -2 * FMMFd[1][1] * chi_mat[1][2]

    Mte_G1 = -2 * FFd[2][2] * chi_mat[2][0]
    Mte_S1 = -2 * FMMFd[2][2] * chi_mat[2][0]

    Mtm_G1 = -2 * FFd[2][2] * chi_mat[2][1]
    Mtm_S1 = -2 * FMMFd[2][2] * chi_mat[2][1]

    Mtt_G1 = -2 * FFd[2][2] * chi_mat[2][2]
    Mtt_S1 = -2 * FMMFd[2][2] * chi_mat[2][2]

    #Columns of the matrix A of the linearised equations related to rN-1, the sum here combines lists into a single list
    #Will need to multily by G0 and S0, and (M0/Tew)
    R_G0_combined = [GellMann_coefficients(R_G0[i]) + [0]*12 for i in range(9)]
    R_S0_combined = [GellMann_coefficients(R_S0[i]) + [0]*12 for i in range(9)]
                     
    Rb_G0_combined = [[0]*9  + GellMann_coefficients(Rb_G0[i]) + [0, 0, 0] for i in range(9)]
    Rb_S0_combined = [[0]*9  + GellMann_coefficients(Rb_S0[i]) + [0, 0, 0] for i in range(9)]

    #Will need to multily by G0 and S0, and (M0/Tew) * fact2
    MR_G0_combined = [[0] * 9 + [0] * 9 + [MRe_G0[i], MRmu_G0[i], MRtau_G0[i]] for i in range(9)]
    MR_S0_combined = [[0] * 9 + [0] * 9 + [MRe_S0[i], MRmu_S0[i], MRtau_S0[i]] for i in range(9)]

    MRbe_G0_combined = [[0]*9 + [0] * 9 + [MRbe_G0[i], MRbmu_G0[i], MRbtau_G0[i]] for i in range(9)]
    MRbe_S0_combined = [[0]*9 + [0] * 9 + [MRbe_S0[i], MRbmu_S0[i], MRbtau_S0[i]] for i in range(9)]

    # Will need to multily by G1 and S1, and (M0/Tew)
    RMe_cols_G1 = GellMann_coefficients(RMe_G1) + GellMann_coefficients(RbMe_G1) + [0,0,0]
    RMmu_cols_G1 = GellMann_coefficients(RMmu_G1) + GellMann_coefficients(RbMmu_G1) + [0,0,0]
    RMtau_cols_G1 = GellMann_coefficients(RMtau_G1) + GellMann_coefficients(RbMtau_G1) + [0,0,0]

    RMe_cols_S1 = GellMann_coefficients(RMe_S1) + GellMann_coefficients(RbMe_S1) + [0,0,0]
    RMmu_cols_S1 = GellMann_coefficients(RMmu_S1) + GellMann_coefficients(RbMmu_S1) + [0,0,0]
    RMtau_cols_S1 = GellMann_coefficients(RMtau_S1) + GellMann_coefficients(RbMtau_S1) + [0,0,0]

    #Will need to multiply by G1 and S1, and (M0/Tew) * fact2
    Me_G1_cols = [0]*9 + [0]*9 + [Mee_G1, Mme_G1, Mte_G1]
    Mmu_G1_cols = [0]*9 + [0]*9 + [Mem_G1, Mmm_G1, Mtm_G1]
    Mtau_G1_cols = [0]*9 + [0]*9 + [Met_G1, Mmt_G1, Mtt_G1]

    Me_S1_cols = [0]*9 + [0]*9 + [Mee_S1, Mme_S1, Mte_S1]
    Mmu_S1_cols = [0]*9 + [0]*9 + [Mem_S1, Mmm_S1, Mtm_S1]
    Mtau_S1_cols = [0]*9 + [0]*9 + [Met_S1, Mmt_S1, Mtt_S1]


    A_G0 = np.transpose(np.array(R_G0_combined[0]+R_G0_combined[1]+R_G0_combined[2]+R_G0_combined[3]+R_G0_combined[4]+R_G0_combined[5]+R_G0_combined[6]+R_G0_combined[7]+R_G0_combined[8] +
                                 Rb_G0_combined[0]+Rb_G0_combined[1]+Rb_G0_combined[2]+Rb_G0_combined[3]+Rb_G0_combined[4]+Rb_G0_combined[5]+Rb_G0_combined[6]+Rb_G0_combined[7]+Rb_G0_combined[8] +
                                 [0]*21 + [0]*21 + [0]*21, dtype=np.complex128).reshape(21,21))
    
    A_S0 = np.transpose(np.array(R_S0_combined[0]+R_S0_combined[1]+R_S0_combined[2]+R_S0_combined[3]+R_S0_combined[4]+R_S0_combined[5]+R_S0_combined[6]+R_S0_combined[7]+R_S0_combined[8] +
                                 Rb_S0_combined[0]+Rb_S0_combined[1]+Rb_S0_combined[2]+Rb_S0_combined[3]+Rb_S0_combined[4]+Rb_S0_combined[5]+Rb_S0_combined[6]+Rb_S0_combined[7]+Rb_S0_combined[8] +
                                 [0]*21 + [0]*21 + [0]*21, dtype=np.complex128).reshape(21,21))
    
    A_G0_fact2 = np.transpose(np.array(MR_G0_combined[0]+MR_G0_combined[1]+MR_G0_combined[2]+MR_G0_combined[3]+MR_G0_combined[4]+MR_G0_combined[5]+MR_G0_combined[6]+MR_G0_combined[7]+MR_G0_combined[8] +
                                       MRbe_G0_combined[0]+MRbe_G0_combined[1]+MRbe_G0_combined[2]+MRbe_G0_combined[3]+MRbe_G0_combined[4]+MRbe_G0_combined[5]+MRbe_G0_combined[6]+MRbe_G0_combined[7]+MRbe_G0_combined[8] +
                                       [0]*21 + [0]*21 + [0]*21, dtype=np.complex128).reshape(21,21))
    
    A_S0_fact2 = np.transpose(np.array(MR_S0_combined[0]+MR_S0_combined[1]+MR_S0_combined[2]+MR_S0_combined[3]+MR_S0_combined[4]+MR_S0_combined[5]+MR_S0_combined[6]+MR_S0_combined[7]+MR_S0_combined[8] +
                                       MRbe_S0_combined[0]+MRbe_S0_combined[1]+MRbe_S0_combined[2]+MRbe_S0_combined[3]+MRbe_S0_combined[4]+MRbe_S0_combined[5]+MRbe_S0_combined[6]+MRbe_S0_combined[7]+MRbe_S0_combined[8] +
                                       [0]*21 + [0]*21 + [0]*21, dtype=np.complex128).reshape(21,21))
    
    A_G1 = np.transpose(np.array([0]*21 + [0]*21 + [0]*21 + [0]*21 + [0]*21 + [0]*21 + [0]*21 + [0]*21 +[0]*21 + 
                                 [0]*21 + [0]*21 + [0]*21 + [0]*21 + [0]*21 + [0]*21 + [0]*21 + [0]*21 +[0]*21 + 
                                 RMe_cols_G1 + RMmu_cols_G1 + RMtau_cols_G1, dtype=np.complex128).reshape(21,21))
    
    A_S1 = np.transpose(np.array([0]*21 + [0]*21 + [0]*21 + [0]*21 + [0]*21 + [0]*21 + [0]*21 + [0]*21 +[0]*21 +
                                 [0]*21 + [0]*21 + [0]*21 + [0]*21 + [0]*21 + [0]*21 + [0]*21 + [0]*21 +[0]*21 +
                                 RMe_cols_S1 + RMmu_cols_S1 + RMtau_cols_S1, dtype=np.complex128).reshape(21,21))
    
    A_G1_fact2 = np.transpose(np.array([0]*21 + [0]*21 + [0]*21 + [0]*21 + [0]*21 + [0]*21 + [0]*21 + [0]*21 +[0]*21 + 
                                       [0]*21 + [0]*21 + [0]*21 + [0]*21 + [0]*21 + [0]*21 + [0]*21 + [0]*21 +[0]*21 + 
                                       Me_G1_cols + Mmu_G1_cols + Mtau_G1_cols, dtype=np.complex128).reshape(21,21))
    
    A_S1_fact2 = np.transpose(np.array([0]*21 + [0]*21 + [0]*21 + [0]*21 + [0]*21 + [0]*21 + [0]*21 + [0]*21 +[0]*21 + 
                                       [0]*21 + [0]*21 + [0]*21 + [0]*21 + [0]*21 + [0]*21 + [0]*21 + [0]*21 +[0]*21 + 
                                       Me_S1_cols + Mmu_S1_cols + Mtau_S1_cols, dtype=np.complex128).reshape(21,21))

    """
    Hamiltionian
    """

    #Zeroth-order contribution
    #Needs to be multiplied by h0 and x*x/(Tew*Tew)
    Ham_RN_0 = 0.5*Dm2_mat
    Ham_RNb_0 = 0.5*Dm2_mat

    #LNC contribution
    #Needs to be multiplied by hLNC
    Ham_RN_LNC = FdF
    Ham_RNb_LNC = np.conjugate(FdF)

    #LNC contribution
    #Needs to be multiplied by hLNV/M1**2
    Ham_RN_LNV = M_Ftrans_Fstar_M
    Ham_RNb_LNV = np.conjugate(M_Ftrans_Fstar_M)

    """
    We need to compute 18 3x3 matrices for the linear terms in cN and cbN, 3 3x3 matrices for the chemical potential, one 3x3 matrix for the constant term
    and 27 3x3 matrices for the non-linear terms. We compute each column and then assemble the full matrix.
    For the time being we neglect the non-lienar terms.
    """
    
    #Terms in the equations for RN and RbN that are linear in RN and RbN that depend on x and cannot be calculated outside this function
    R_ham_0   = [- 1j * commutator(Ham_RN_0, l) for l in L]
    Rb_ham_0  = [- 1j * commutator(Ham_RNb_0, l) for l in L]

    R_ham_LNC   = [- 1j * commutator(Ham_RN_LNC, l) for l in L]
    Rb_ham_LNC  = [- 1j * commutator(Ham_RNb_LNC, l) for l in L]

    R_ham_LNV   = [- 1j * commutator(Ham_RN_LNV, l) for l in L]
    Rb_ham_LNV  = [- 1j * commutator(Ham_RNb_LNV, l) for l in L]


    #Columns of the matrix A of the linearised equations related to rN, the sum here combines lists into a single list
    R_H_0 = [GellMann_coefficients(R_ham_0[i]) + [0]*12 for i in range(9)]
    Rb_H_0 = [[0]*9 + GellMann_coefficients(Rb_ham_0[i]) + [0]*3 for i in range(9)]

    R_H_LNC = [GellMann_coefficients(R_ham_LNC[i]) + [0]*12 for i in range(9)]
    Rb_H_LNC = [[0]*9 + GellMann_coefficients(Rb_ham_LNC[i]) + [0]*3 for i in range(9)]
    
    R_H_LNV = [GellMann_coefficients(R_ham_LNV[i]) + [0]*12 for i in range(9)]
    Rb_H_LNV = [[0]*9 + GellMann_coefficients(Rb_ham_LNV[i]) + [0]*3 for i in range(9)]
    
    A_ham_0 = np.transpose(np.array(R_H_0[0]+R_H_0[1]+R_H_0[2]+R_H_0[3]+R_H_0[4]+R_H_0[5]+R_H_0[6]+R_H_0[7]+R_H_0[8] +
                                    Rb_H_0[0]+Rb_H_0[1]+Rb_H_0[2]+Rb_H_0[3]+Rb_H_0[4]+Rb_H_0[5]+Rb_H_0[6]+Rb_H_0[7]+Rb_H_0[8] +
                                    [0]*21 + [0]*21 + [0]*21, dtype=np.complex128).reshape(21,21))
    

    A_ham_LNC = np.transpose(np.array(R_H_LNC[0]+R_H_LNC[1]+R_H_LNC[2]+R_H_LNC[3]+R_H_LNC[4]+R_H_LNC[5]+R_H_LNC[6]+R_H_LNC[7]+R_H_LNC[8] +
                                  Rb_H_LNC[0]+Rb_H_LNC[1]+Rb_H_LNC[2]+Rb_H_LNC[3]+Rb_H_LNC[4]+Rb_H_LNC[5]+Rb_H_LNC[6]+Rb_H_LNC[7]+Rb_H_LNC[8] +
                                  [0]*21 + [0]*21 + [0]*21, dtype=np.complex128).reshape(21,21))
    
    A_ham_LNV = np.transpose(np.array([R_H_LNV[0]+R_H_LNV[1]+R_H_LNV[2]+R_H_LNV[3]+R_H_LNV[4]+R_H_LNV[5]+R_H_LNV[6]+R_H_LNV[7]+R_H_LNV[8] +
                                       Rb_H_LNV[0]+Rb_H_LNV[1]+Rb_H_LNV[2]+Rb_H_LNV[3]+Rb_H_LNV[4]+Rb_H_LNV[5]+Rb_H_LNV[6]+Rb_H_LNV[7]+Rb_H_LNV[8] +
                                       [0]*21 + [0]*21 + [0]*21], dtype=np.complex128).reshape(21,21))
    

    #We are writing the equations in terms of rn-1, therefore, there is an additional matrix, which is simply L, and needs to be multiplied by fact only
    if n3_decoupled:
        L_gm = [GellMann_coefficients_decoupled(l) + GellMann_coefficients_decoupled(l) + [0]*3 for l in L]#[_gm_proj(l) + _gm_proj(l) + [0]*3 for l in L]
    else:
        L_gm = [GellMann_coefficients(l) + GellMann_coefficients(l) + [0]*3 for l in L]

    A_L  = np.transpose(np.array(L_gm[0]+L_gm[1]+L_gm[2]+L_gm[3]+L_gm[4]+L_gm[5]+L_gm[6]+L_gm[7]+L_gm[8] +
                                  L_gm[0]+L_gm[1]+L_gm[2]+L_gm[3]+L_gm[4]+L_gm[5]+L_gm[6]+L_gm[7]+L_gm[8] +
                                  [0]*21 + [0]*21 + [0]*21, dtype=np.complex128).reshape(21,21))
    

    """
    #Non-linear terms
    """
    
    #Terms non linear in RN-1 and chemical potential, in the equation for RN-1. Will need to multiply by G2 or S2 and M0/Tew
    NL_RN_G2_e    = [anticommutator(Fd_chie_F, l) for l in L]
    NL_RN_G2_mu   = [anticommutator(Fd_chimu_F, l) for l in L]
    NL_RN_G2_tau  = [anticommutator(Fd_chitau_F, l) for l in L]

    NL_RN_S2_e    = [-anticommutator(MF_chie_FstarM, l) for l in L]
    NL_RN_S2_mu   = [-anticommutator(MF_chimu_FstarM, l) for l in L]
    NL_RN_S2_tau  = [-anticommutator(MF_chitau_FstarM, l) for l in L]

    #Terms non linear in RNb-1 and chemical potential, in the equation for RNb-1. Will need to multiply by G2 or S2 and M0/Tew
    NL_RNb_G2_e    = [-anticommutator(np.conjugate(Fd_chie_F), l) for l in L]
    NL_RNb_G2_mu   = [-anticommutator(np.conjugate(Fd_chimu_F), l) for l in L]
    NL_RNb_G2_tau  = [-anticommutator(np.conjugate(Fd_chitau_F), l) for l in L]

    NL_RNb_S2_e    = [anticommutator(np.conjugate(MF_chie_FstarM), l) for l in L]
    NL_RNb_S2_mu   = [anticommutator(np.conjugate(MF_chimu_FstarM), l) for l in L]
    NL_RNb_S2_tau  = [anticommutator(np.conjugate(MF_chitau_FstarM), l) for l in L]

    #We perform the Gell-Mann transformation of all these matrices
    Bvec_G2_e      = [GellMann_coefficients(r) for r in NL_RN_G2_e]
    Bvec_G2_mu     = [GellMann_coefficients(r) for r in NL_RN_G2_mu]
    Bvec_G2_tau    = [GellMann_coefficients(r) for r in NL_RN_G2_tau]

    Bvec_S2_e      = [GellMann_coefficients(r) for r in NL_RN_S2_e]
    Bvec_S2_mu     = [GellMann_coefficients(r) for r in NL_RN_S2_mu]
    Bvec_S2_tau    = [GellMann_coefficients(r) for r in NL_RN_S2_tau]

    Bbvec_G2_e     = [GellMann_coefficients(r) for r in NL_RNb_G2_e]
    Bbvec_G2_mu    = [GellMann_coefficients(r) for r in NL_RNb_G2_mu]
    Bbvec_G2_tau   = [GellMann_coefficients(r) for r in NL_RNb_G2_tau]

    Bbvec_S2_e     = [GellMann_coefficients(r) for r in NL_RNb_S2_e]
    Bbvec_S2_mu    = [GellMann_coefficients(r) for r in NL_RNb_S2_mu]
    Bbvec_S2_tau   = [GellMann_coefficients(r) for r in NL_RNb_S2_tau]

    #At this point Bkae_G2 = Bvec_G2_e[a][k] etc., we define Bvec as a list of matrices B_\kappa, with entries (B_\kappa)_{\rho, \sigma}
    Bvec_G2, Bvec_S2 = list(), list()

    for k in range(0, 9):
        b_G2    = np.array([0]*18 + [Bvec_G2_e[0][k]] + [Bvec_G2_mu[0][k]] + [Bvec_G2_tau[0][k]]+
                           [0]*18 + [Bvec_G2_e[1][k]] + [Bvec_G2_mu[1][k]] + [Bvec_G2_tau[1][k]]+
                           [0]*18 + [Bvec_G2_e[2][k]] + [Bvec_G2_mu[2][k]] + [Bvec_G2_tau[2][k]]+
                           [0]*18 + [Bvec_G2_e[3][k]] + [Bvec_G2_mu[3][k]] + [Bvec_G2_tau[3][k]]+
                           [0]*18 + [Bvec_G2_e[4][k]] + [Bvec_G2_mu[4][k]] + [Bvec_G2_tau[4][k]]+
                           [0]*18 + [Bvec_G2_e[5][k]] + [Bvec_G2_mu[5][k]] + [Bvec_G2_tau[5][k]]+
                           [0]*18 + [Bvec_G2_e[6][k]] + [Bvec_G2_mu[6][k]] + [Bvec_G2_tau[6][k]]+
                           [0]*18 + [Bvec_G2_e[7][k]] + [Bvec_G2_mu[7][k]] + [Bvec_G2_tau[7][k]]+
                           [0]*18 + [Bvec_G2_e[8][k]] + [Bvec_G2_mu[8][k]] + [Bvec_G2_tau[8][k]]+
                           [0]*18 + [0]*3+
                           [0]*18 + [0]*3+
                           [0]*18 + [0]*3+
                           [0]*18 + [0]*3+
                           [0]*18 + [0]*3+
                           [0]*18 + [0]*3+
                           [0]*18 + [0]*3+
                           [0]*18 + [0]*3+
                           [0]*18 + [0]*3+
                           [0]*21 + [0]*21 + [0] * 21, dtype=np.complex128).reshape(21,21)
        b_S2    = np.array([0]*18 + [Bvec_S2_e[0][k]] + [Bvec_S2_mu[0][k]] + [Bvec_S2_tau[0][k]]+
                           [0]*18 + [Bvec_S2_e[1][k]] + [Bvec_S2_mu[1][k]] + [Bvec_S2_tau[1][k]]+
                           [0]*18 + [Bvec_S2_e[2][k]] + [Bvec_S2_mu[2][k]] + [Bvec_S2_tau[2][k]]+
                           [0]*18 + [Bvec_S2_e[3][k]] + [Bvec_S2_mu[3][k]] + [Bvec_S2_tau[3][k]]+
                           [0]*18 + [Bvec_S2_e[4][k]] + [Bvec_S2_mu[4][k]] + [Bvec_S2_tau[4][k]]+
                           [0]*18 + [Bvec_S2_e[5][k]] + [Bvec_S2_mu[5][k]] + [Bvec_S2_tau[5][k]]+
                           [0]*18 + [Bvec_S2_e[6][k]] + [Bvec_S2_mu[6][k]] + [Bvec_S2_tau[6][k]]+
                           [0]*18 + [Bvec_S2_e[7][k]] + [Bvec_S2_mu[7][k]] + [Bvec_S2_tau[7][k]]+
                           [0]*18 + [Bvec_S2_e[8][k]] + [Bvec_S2_mu[8][k]] + [Bvec_S2_tau[8][k]]+
                           [0]*18 + [0]*3+
                           [0]*18 + [0]*3+
                           [0]*18 + [0]*3+
                           [0]*18 + [0]*3+
                           [0]*18 + [0]*3+
                           [0]*18 + [0]*3+
                           [0]*18 + [0]*3+
                           [0]*18 + [0]*3+
                           [0]*18 + [0]*3+
                           [0]*21 + [0]*21 + [0] * 21, dtype=np.complex128).reshape(21,21)
        Bvec_G2.append(b_G2 + np.transpose(b_G2))
        Bvec_S2.append(b_S2 + np.transpose(b_S2))

    for k in range(0, 9):
        bb_G2    = np.array([0]*18 + [0]*3+
                            [0]*18 + [0]*3+
                            [0]*18 + [0]*3+
                            [0]*18 + [0]*3+
                            [0]*18 + [0]*3+
                            [0]*18 + [0]*3+
                            [0]*18 + [0]*3+
                            [0]*18 + [0]*3+
                            [0]*18 + [0]*3+ 
                            [0]*18 + [Bbvec_G2_e[0][k]] + [Bbvec_G2_mu[0][k]] + [Bbvec_G2_tau[0][k]]+
                            [0]*18 + [Bbvec_G2_e[1][k]] + [Bbvec_G2_mu[1][k]] + [Bbvec_G2_tau[1][k]]+
                            [0]*18 + [Bbvec_G2_e[2][k]] + [Bbvec_G2_mu[2][k]] + [Bbvec_G2_tau[2][k]]+
                            [0]*18 + [Bbvec_G2_e[3][k]] + [Bbvec_G2_mu[3][k]] + [Bbvec_G2_tau[3][k]]+
                            [0]*18 + [Bbvec_G2_e[4][k]] + [Bbvec_G2_mu[4][k]] + [Bbvec_G2_tau[4][k]]+
                            [0]*18 + [Bbvec_G2_e[5][k]] + [Bbvec_G2_mu[5][k]] + [Bbvec_G2_tau[5][k]]+
                            [0]*18 + [Bbvec_G2_e[6][k]] + [Bbvec_G2_mu[6][k]] + [Bbvec_G2_tau[6][k]]+
                            [0]*18 + [Bbvec_G2_e[7][k]] + [Bbvec_G2_mu[7][k]] + [Bbvec_G2_tau[7][k]]+
                            [0]*18 + [Bbvec_G2_e[8][k]] + [Bbvec_G2_mu[8][k]] + [Bbvec_G2_tau[8][k]]+
                            [0]*21 + [0]*21 + [0]*21, dtype=np.complex128).reshape(21,21)
        bb_S2    = np.array([0]*18 + [0]*3+
                            [0]*18 + [0]*3+
                            [0]*18 + [0]*3+
                            [0]*18 + [0]*3+
                            [0]*18 + [0]*3+
                            [0]*18 + [0]*3+
                            [0]*18 + [0]*3+
                            [0]*18 + [0]*3+
                            [0]*18 + [0]*3+
                            [0]*18 + [Bbvec_S2_e[0][k]] + [Bbvec_S2_mu[0][k]] + [Bbvec_S2_tau[0][k]]+
                            [0]*18 + [Bbvec_S2_e[1][k]] + [Bbvec_S2_mu[1][k]] + [Bbvec_S2_tau[1][k]]+
                            [0]*18 + [Bbvec_S2_e[2][k]] + [Bbvec_S2_mu[2][k]] + [Bbvec_S2_tau[2][k]]+
                            [0]*18 + [Bbvec_S2_e[3][k]] + [Bbvec_S2_mu[3][k]] + [Bbvec_S2_tau[3][k]]+
                            [0]*18 + [Bbvec_S2_e[4][k]] + [Bbvec_S2_mu[4][k]] + [Bbvec_S2_tau[4][k]]+
                            [0]*18 + [Bbvec_S2_e[5][k]] + [Bbvec_S2_mu[5][k]] + [Bbvec_S2_tau[5][k]]+
                            [0]*18 + [Bbvec_S2_e[6][k]] + [Bbvec_S2_mu[6][k]] + [Bbvec_S2_tau[6][k]]+
                            [0]*18 + [Bbvec_S2_e[7][k]] + [Bbvec_S2_mu[7][k]] + [Bbvec_S2_tau[7][k]]+
                            [0]*18 + [Bbvec_S2_e[8][k]] + [Bbvec_S2_mu[8][k]] + [Bbvec_S2_tau[8][k]]+
                            [0]*21 + [0]*21 + [0]*21, dtype=np.complex128).reshape(21,21)
        Bvec_G2.append(bb_G2 + np.transpose(bb_G2))
        Bvec_S2.append(bb_S2 + np.transpose(bb_S2))

    for k in range(0, 3):
        bmu_G2 = np.zeros(21*21).reshape(21,21)
        bmu_S2 = np.zeros(21*21).reshape(21,21)
        Bvec_G2.append(bmu_G2)
        Bvec_S2.append(bmu_S2)


    #Terms non linear in RN-1 and chemical potential, in the equation for the chemical potential. Will need to multiply by G2 or S2, and M0/Tew * fact2
    NL_chem_G2_ee = [f1[0][0]*chi_mat[0, 0] for f1 in F_L_Fdag]
    NL_chem_G2_emu = [f1[0][0]*chi_mat[0, 1] for f1 in F_L_Fdag]
    NL_chem_G2_etau = [f1[0][0]*chi_mat[0, 2] for f1 in F_L_Fdag]

    NL_chem_G2_mue = [f1[1][1]*chi_mat[1, 0] for f1 in F_L_Fdag]
    NL_chem_G2_mumu = [f1[1][1]*chi_mat[1, 1] for f1 in F_L_Fdag]
    NL_chem_G2_mutau = [f1[1][1]*chi_mat[1, 2] for f1 in F_L_Fdag]

    NL_chem_G2_taue = [f1[2][2]*chi_mat[2, 0] for f1 in F_L_Fdag]
    NL_chem_G2_taumu = [f1[2][2]*chi_mat[2, 1] for f1 in F_L_Fdag]
    NL_chem_G2_tautau = [f1[2][2]*chi_mat[2, 2] for f1 in F_L_Fdag]

    NL_chem_S2_ee = [f2[0][0]*chi_mat[0, 0] for f2 in FstarM_L_MFtrans]
    NL_chem_S2_emu = [f2[0][0]*chi_mat[0, 1] for f2 in FstarM_L_MFtrans]
    NL_chem_S2_etau = [f2[0][0]*chi_mat[0, 2] for f2 in FstarM_L_MFtrans]

    NL_chem_S2_mue = [f2[1][1]*chi_mat[1, 0] for f2 in FstarM_L_MFtrans]
    NL_chem_S2_mumu = [f2[1][1]*chi_mat[1, 1] for f2 in FstarM_L_MFtrans]
    NL_chem_S2_mutau = [f2[1][1]*chi_mat[1, 2] for f2 in FstarM_L_MFtrans]

    NL_chem_S2_taue = [f2[2][2]*chi_mat[2, 0] for f2 in FstarM_L_MFtrans]
    NL_chem_S2_taumu = [f2[2][2]*chi_mat[2, 1] for f2 in FstarM_L_MFtrans]
    NL_chem_S2_tautau = [f2[2][2]*chi_mat[2, 2] for f2 in FstarM_L_MFtrans]

    #Terms non linear in RNb-1 and chemical potential, in the equation for the chemical potential. Will need to multiply by G2 or S2, and M0/Tew * fact2
    NLb_chem_G2_ee = [f1[0][0]*chi_mat[0, 0] for f1 in Fstar_L_Ftrans]
    NLb_chem_G2_emu = [f1[0][0]*chi_mat[0, 1] for f1 in Fstar_L_Ftrans]
    NLb_chem_G2_etau = [f1[0][0]*chi_mat[0, 2] for f1 in Fstar_L_Ftrans]

    NLb_chem_G2_mue = [f1[1][1]*chi_mat[1, 0] for f1 in Fstar_L_Ftrans]
    NLb_chem_G2_mumu = [f1[1][1]*chi_mat[1, 1] for f1 in Fstar_L_Ftrans]
    NLb_chem_G2_mutau = [f1[1][1]*chi_mat[1, 2] for f1 in Fstar_L_Ftrans]

    NLb_chem_G2_taue = [f1[2][2]*chi_mat[2, 0] for f1 in Fstar_L_Ftrans]
    NLb_chem_G2_taumu = [f1[2][2]*chi_mat[2, 1] for f1 in Fstar_L_Ftrans]
    NLb_chem_G2_tautau = [f1[2][2]*chi_mat[2, 2] for f1 in Fstar_L_Ftrans]

    NLb_chem_S2_ee = [f2[0][0]*chi_mat[0, 0] for f2 in FM_L_MFdag]
    NLb_chem_S2_emu = [f2[0][0]*chi_mat[0, 1] for f2 in FM_L_MFdag]
    NLb_chem_S2_etau = [f2[0][0]*chi_mat[0, 2] for f2 in FM_L_MFdag]

    NLb_chem_S2_mue = [f2[1][1]*chi_mat[1, 0] for f2 in FM_L_MFdag]
    NLb_chem_S2_mumu = [f2[1][1]*chi_mat[1, 1] for f2 in FM_L_MFdag]
    NLb_chem_S2_mutau = [f2[1][1]*chi_mat[1, 2] for f2 in FM_L_MFdag]

    NLb_chem_S2_taue = [f2[2][2]*chi_mat[2, 0] for f2 in FM_L_MFdag]
    NLb_chem_S2_taumu = [f2[2][2]*chi_mat[2, 1] for f2 in FM_L_MFdag]
    NLb_chem_S2_tautau = [f2[2][2]*chi_mat[2, 2] for f2 in FM_L_MFdag]

    Bvec_G2_fact2, Bvec_S2_fact2 = list(), list()
    for k in range(0, 18):
        b_G2 = np.zeros(21*21).reshape(21,21)
        b_S2 = np.zeros(21*21).reshape(21,21)
        Bvec_G2_fact2.append(b_G2)
        Bvec_S2_fact2.append(b_S2)
    
    NLbe_G2 = np.array([0]*18 + [0]*3+
                       [0]*18 + [0]*3+
                       [0]*18 + [0]*3+
                       [0]*18 + [0]*3+
                       [0]*18 + [0]*3+
                       [0]*18 + [0]*3+
                       [0]*18 + [0]*3+
                       [0]*18 + [0]*3+
                       [0]*18 + [0]*3+
                       [0]*18 + [0]*3+
                       [0]*18 + [0]*3+
                       [0]*18 + [0]*3+
                       [0]*18 + [0]*3+
                       [0]*18 + [0]*3+
                       [0]*18 + [0]*3+
                       [0]*18 + [0]*3+
                       [0]*18 + [0]*3+
                       [0]*18 + [0]*3+
                       [NL_chem_G2_ee[0]] + [NL_chem_G2_ee[1]] + [NL_chem_G2_ee[2]] + [NL_chem_G2_ee[3]] + [NL_chem_G2_ee[4]] + [NL_chem_G2_ee[5]] + [NL_chem_G2_ee[6]] + [NL_chem_G2_ee[7]] + [NL_chem_G2_ee[8]] +
                       [NLb_chem_G2_ee[0]] + [NLb_chem_G2_ee[1]] + [NLb_chem_G2_ee[2]] + [NLb_chem_G2_ee[3]] + [NLb_chem_G2_ee[4]] + [NLb_chem_G2_ee[5]] + [NLb_chem_G2_ee[6]] + [NLb_chem_G2_ee[7]] + [NLb_chem_G2_ee[8]] +
                       [0]*3 +
                       [NL_chem_G2_emu[0]] + [NL_chem_G2_emu[1]] + [NL_chem_G2_emu[2]] + [NL_chem_G2_emu[3]] + [NL_chem_G2_emu[4]] + [NL_chem_G2_emu[5]] + [NL_chem_G2_emu[6]] + [NL_chem_G2_emu[7]] + [NL_chem_G2_emu[8]] +
                       [NLb_chem_G2_emu[0]] + [NLb_chem_G2_emu[1]] + [NLb_chem_G2_emu[2]] + [NLb_chem_G2_emu[3]] + [NLb_chem_G2_emu[4]] + [NLb_chem_G2_emu[5]] + [NLb_chem_G2_emu[6]] + [NLb_chem_G2_emu[7]] + [NLb_chem_G2_emu[8]] +
                       [0]*3 +
                       [NL_chem_G2_etau[0]] + [NL_chem_G2_etau[1]] + [NL_chem_G2_etau[2]] + [NL_chem_G2_etau[3]] + [NL_chem_G2_etau[4]] + [NL_chem_G2_etau[5]] + [NL_chem_G2_etau[6]] + [NL_chem_G2_etau[7]] + [NL_chem_G2_etau[8]] +
                       [NLb_chem_G2_etau[0]] + [NLb_chem_G2_etau[1]] + [NLb_chem_G2_etau[2]] + [NLb_chem_G2_etau[3]] + [NLb_chem_G2_etau[4]] + [NLb_chem_G2_etau[5]] + [NLb_chem_G2_etau[6]] + [NLb_chem_G2_etau[7]] + [NLb_chem_G2_etau[8]] +
                       [0]*3,
                       dtype=np.complex128).reshape(21,21)
    
    NLbe_S2 = np.array([0]*18 + [0]*3+
                       [0]*18 + [0]*3+
                       [0]*18 + [0]*3+
                       [0]*18 + [0]*3+
                       [0]*18 + [0]*3+
                       [0]*18 + [0]*3+
                       [0]*18 + [0]*3+
                       [0]*18 + [0]*3+
                       [0]*18 + [0]*3+
                       [0]*18 + [0]*3+
                       [0]*18 + [0]*3+
                       [0]*18 + [0]*3+
                       [0]*18 + [0]*3+
                       [0]*18 + [0]*3+
                       [0]*18 + [0]*3+
                       [0]*18 + [0]*3+
                       [0]*18 + [0]*3+
                       [0]*18 + [0]*3+
                       [NL_chem_S2_ee[0]] + [NL_chem_S2_ee[1]] + [NL_chem_S2_ee[2]] + [NL_chem_S2_ee[3]] + [NL_chem_S2_ee[4]] + [NL_chem_S2_ee[5]] + [NL_chem_S2_ee[6]] + [NL_chem_S2_ee[7]] + [NL_chem_S2_ee[8]] +
                       [NLb_chem_S2_ee[0]] + [NLb_chem_S2_ee[1]] + [NLb_chem_S2_ee[2]] + [NLb_chem_S2_ee[3]] + [NLb_chem_S2_ee[4]] + [NLb_chem_S2_ee[5]] + [NLb_chem_S2_ee[6]] + [NLb_chem_S2_ee[7]] + [NLb_chem_S2_ee[8]] +
                       [0]*3 +
                       [NL_chem_S2_emu[0]] + [NL_chem_S2_emu[1]] + [NL_chem_S2_emu[2]] + [NL_chem_S2_emu[3]] + [NL_chem_S2_emu[4]] + [NL_chem_S2_emu[5]] + [NL_chem_S2_emu[6]] + [NL_chem_S2_emu[7]] + [NL_chem_S2_emu[8]] +
                       [NLb_chem_S2_emu[0]] + [NLb_chem_S2_emu[1]] + [NLb_chem_S2_emu[2]] + [NLb_chem_S2_emu[3]] + [NLb_chem_S2_emu[4]] + [NLb_chem_S2_emu[5]] + [NLb_chem_S2_emu[6]] + [NLb_chem_S2_emu[7]] + [NLb_chem_S2_emu[8]] +
                       [0]*3 +
                       [NL_chem_S2_etau[0]] + [NL_chem_S2_etau[1]] + [NL_chem_S2_etau[2]] + [NL_chem_S2_etau[3]] + [NL_chem_S2_etau[4]] + [NL_chem_S2_etau[5]] + [NL_chem_S2_etau[6]] + [NL_chem_S2_etau[7]] + [NL_chem_S2_etau[8]] +
                       [NLb_chem_S2_etau[0]] + [NLb_chem_S2_etau[1]] + [NLb_chem_S2_etau[2]] + [NLb_chem_S2_etau[3]] + [NLb_chem_S2_etau[4]] + [NLb_chem_S2_etau[5]] + [NLb_chem_S2_etau[6]] + [NLb_chem_S2_etau[7]] + [NLb_chem_S2_etau[8]] +
                       [0]*3,
                       dtype=np.complex128).reshape(21,21)

    Bvec_G2_fact2.append((NLbe_G2) + np.transpose(NLbe_G2))
    Bvec_S2_fact2.append((NLbe_S2) + np.transpose(NLbe_S2))

    NLbmu_G2 = np.array([0]*18 + [0]*3+
                        [0]*18 + [0]*3+
                        [0]*18 + [0]*3+
                        [0]*18 + [0]*3+
                        [0]*18 + [0]*3+
                        [0]*18 + [0]*3+
                        [0]*18 + [0]*3+
                        [0]*18 + [0]*3+
                        [0]*18 + [0]*3+
                        [0]*18 + [0]*3+
                        [0]*18 + [0]*3+
                        [0]*18 + [0]*3+
                        [0]*18 + [0]*3+
                        [0]*18 + [0]*3+
                        [0]*18 + [0]*3+
                        [0]*18 + [0]*3+
                        [0]*18 + [0]*3+
                        [0]*18 + [0]*3+
                        [NL_chem_G2_mue[0]] + [NL_chem_G2_mue[1]] + [NL_chem_G2_mue[2]] + [NL_chem_G2_mue[3]] + [NL_chem_G2_mue[4]] + [NL_chem_G2_mue[5]] + [NL_chem_G2_mue[6]] + [NL_chem_G2_mue[7]] + [NL_chem_G2_mue[8]] +
                        [NLb_chem_G2_mue[0]] + [NLb_chem_G2_mue[1]] + [NLb_chem_G2_mue[2]] + [NLb_chem_G2_mue[3]] + [NLb_chem_G2_mue[4]] + [NLb_chem_G2_mue[5]] + [NLb_chem_G2_mue[6]] + [NLb_chem_G2_mue[7]] + [NLb_chem_G2_mue[8]] +
                        [0]*3 + 
                        [NL_chem_G2_mumu[0]] + [NL_chem_G2_mumu[1]] + [NL_chem_G2_mumu[2]] + [NL_chem_G2_mumu[3]] + [NL_chem_G2_mumu[4]] + [NL_chem_G2_mumu[5]] + [NL_chem_G2_mumu[6]] + [NL_chem_G2_mumu[7]] + [NL_chem_G2_mumu[8]] +
                        [NLb_chem_G2_mumu[0]] + [NLb_chem_G2_mumu[1]] + [NLb_chem_G2_mumu[2]] + [NLb_chem_G2_mumu[3]] + [NLb_chem_G2_mumu[4]] + [NLb_chem_G2_mumu[5]] + [NLb_chem_G2_mumu[6]] + [NLb_chem_G2_mumu[7]] + [NLb_chem_G2_mumu[8]] +
                        [0]*3 + 
                        [NL_chem_G2_mutau[0]] + [NL_chem_G2_mutau[1]] + [NL_chem_G2_mutau[2]] + [NL_chem_G2_mutau[3]] + [NL_chem_G2_mutau[4]] + [NL_chem_G2_mutau[5]] + [NL_chem_G2_mutau[6]] + [NL_chem_G2_mutau[7]] + [NL_chem_G2_mutau[8]] + 
                        [NLb_chem_G2_mutau[0]] + [NLb_chem_G2_mutau[1]] + [NLb_chem_G2_mutau[2]] + [NLb_chem_G2_mutau[3]] + [NLb_chem_G2_mutau[4]] + [NLb_chem_G2_mutau[5]] + [NLb_chem_G2_mutau[6]] + [NLb_chem_G2_mutau[7]] + [NLb_chem_G2_mutau[8]] +
                        [0]*3,
                        dtype=np.complex128).reshape(21,21)
    
    NLbmu_S2 = np.array([0]*18 + [0]*3+
                        [0]*18 + [0]*3+
                        [0]*18 + [0]*3+
                        [0]*18 + [0]*3+
                        [0]*18 + [0]*3+
                        [0]*18 + [0]*3+
                        [0]*18 + [0]*3+
                        [0]*18 + [0]*3+
                        [0]*18 + [0]*3+
                        [0]*18 + [0]*3+
                        [0]*18 + [0]*3+
                        [0]*18 + [0]*3+
                        [0]*18 + [0]*3+
                        [0]*18 + [0]*3+
                        [0]*18 + [0]*3+
                        [0]*18 + [0]*3+
                        [0]*18 + [0]*3+
                        [0]*18 + [0]*3+
                        [NL_chem_S2_mue[0]] + [NL_chem_S2_mue[1]] + [NL_chem_S2_mue[2]] + [NL_chem_S2_mue[3]] + [NL_chem_S2_mue[4]] + [NL_chem_S2_mue[5]] + [NL_chem_S2_mue[6]] + [NL_chem_S2_mue[7]] + [NL_chem_S2_mue[8]] +
                        [NLb_chem_S2_mue[0]] + [NLb_chem_S2_mue[1]] + [NLb_chem_S2_mue[2]] + [NLb_chem_S2_mue[3]] + [NLb_chem_S2_mue[4]] + [NLb_chem_S2_mue[5]] + [NLb_chem_S2_mue[6]] + [NLb_chem_S2_mue[7]] + [NLb_chem_S2_mue[8]] +
                        [0]*3 + 
                        [NL_chem_S2_mumu[0]] + [NL_chem_S2_mumu[1]] + [NL_chem_S2_mumu[2]] + [NL_chem_S2_mumu[3]] + [NL_chem_S2_mumu[4]] + [NL_chem_S2_mumu[5]] + [NL_chem_S2_mumu[6]] + [NL_chem_S2_mumu[7]] + [NL_chem_S2_mumu[8]] +
                        [NLb_chem_S2_mumu[0]] + [NLb_chem_S2_mumu[1]] + [NLb_chem_S2_mumu[2]] + [NLb_chem_S2_mumu[3]] + [NLb_chem_S2_mumu[4]] + [NLb_chem_S2_mumu[5]] + [NLb_chem_S2_mumu[6]] + [NLb_chem_S2_mumu[7]] + [NLb_chem_S2_mumu[8]] +
                        [0]*3 + 
                        [NL_chem_S2_mutau[0]] + [NL_chem_S2_mutau[1]] + [NL_chem_S2_mutau[2]] + [NL_chem_S2_mutau[3]] + [NL_chem_S2_mutau[4]] + [NL_chem_S2_mutau[5]] + [NL_chem_S2_mutau[6]] + [NL_chem_S2_mutau[7]] + [NL_chem_S2_mutau[8]] +
                        [NLb_chem_S2_mutau[0]] + [NLb_chem_S2_mutau[1]] + [NLb_chem_S2_mutau[2]] + [NLb_chem_S2_mutau[3]] + [NLb_chem_S2_mutau[4]] + [NLb_chem_S2_mutau[5]] + [NLb_chem_S2_mutau[6]] + [NLb_chem_S2_mutau[7]] + [NLb_chem_S2_mutau[8]] +
                        [0]*3,
                        dtype=np.complex128).reshape(21,21)

    Bvec_G2_fact2.append(NLbmu_G2 + np.transpose(NLbmu_G2))
    Bvec_S2_fact2.append(NLbmu_S2 + np.transpose(NLbmu_S2))

    NLbtau_G2 = np.array([0]*18 + [0]*3+
                         [0]*18 + [0]*3+
                         [0]*18 + [0]*3+
                         [0]*18 + [0]*3+
                         [0]*18 + [0]*3+
                         [0]*18 + [0]*3+
                         [0]*18 + [0]*3+
                         [0]*18 + [0]*3+
                         [0]*18 + [0]*3+
                         [0]*18 + [0]*3+
                         [0]*18 + [0]*3+
                         [0]*18 + [0]*3+
                         [0]*18 + [0]*3+
                         [0]*18 + [0]*3+
                         [0]*18 + [0]*3+
                         [0]*18 + [0]*3+
                         [0]*18 + [0]*3+
                         [0]*18 + [0]*3+
                         [NL_chem_G2_taue[0]] + [NL_chem_G2_taue[1]] + [NL_chem_G2_taue[2]] + [NL_chem_G2_taue[3]] + [NL_chem_G2_taue[4]] + [NL_chem_G2_taue[5]] + [NL_chem_G2_taue[6]] + [NL_chem_G2_taue[7]] + [NL_chem_G2_taue[8]] +
                         [NLb_chem_G2_taue[0]] + [NLb_chem_G2_taue[1]] + [NLb_chem_G2_taue[2]] + [NLb_chem_G2_taue[3]] + [NLb_chem_G2_taue[4]] + [NLb_chem_G2_taue[5]] + [NLb_chem_G2_taue[6]] + [NLb_chem_G2_taue[7]] + [NLb_chem_G2_taue[8]] +
                         [0]*3 + 
                         [NL_chem_G2_taumu[0]] + [NL_chem_G2_taumu[1]] + [NL_chem_G2_taumu[2]] + [NL_chem_G2_taumu[3]] + [NL_chem_G2_taumu[4]] + [NL_chem_G2_taumu[5]] + [NL_chem_G2_taumu[6]] + [NL_chem_G2_taumu[7]] + [NL_chem_G2_taumu[8]] +
                         [NLb_chem_G2_taumu[0]] + [NLb_chem_G2_taumu[1]] + [NLb_chem_G2_taumu[2]] + [NLb_chem_G2_taumu[3]] + [NLb_chem_G2_taumu[4]] + [NLb_chem_G2_taumu[5]] + [NLb_chem_G2_taumu[6]] + [NLb_chem_G2_taumu[7]] + [NLb_chem_G2_taumu[8]] +
                         [0]*3 + 
                         [NL_chem_G2_tautau[0]] + [NL_chem_G2_tautau[1]] + [NL_chem_G2_tautau[2]] + [NL_chem_G2_tautau[3]] + [NL_chem_G2_tautau[4]] + [NL_chem_G2_tautau[5]] + [NL_chem_G2_tautau[6]] + [NL_chem_G2_tautau[7]] + [NL_chem_G2_tautau[8]] +
                         [NLb_chem_G2_tautau[0]] + [NLb_chem_G2_tautau[1]] + [NLb_chem_G2_tautau[2]] + [NLb_chem_G2_tautau[3]] + [NLb_chem_G2_tautau[4]] + [NLb_chem_G2_tautau[5]] + [NLb_chem_G2_tautau[6]] + [NLb_chem_G2_tautau[7]] + [NLb_chem_G2_tautau[8]] +
                         [0]*3,
                         dtype=np.complex128).reshape(21,21)
        
    NLbtau_S2 = np.array([0]*18 + [0]*3+
                         [0]*18 + [0]*3+
                         [0]*18 + [0]*3+
                         [0]*18 + [0]*3+
                         [0]*18 + [0]*3+
                         [0]*18 + [0]*3+
                         [0]*18 + [0]*3+
                         [0]*18 + [0]*3+
                         [0]*18 + [0]*3+
                         [0]*18 + [0]*3+
                         [0]*18 + [0]*3+
                         [0]*18 + [0]*3+
                         [0]*18 + [0]*3+
                         [0]*18 + [0]*3+
                         [0]*18 + [0]*3+
                         [0]*18 + [0]*3+
                         [0]*18 + [0]*3+
                         [0]*18 + [0]*3+
                         [NL_chem_S2_taue[0]] + [NL_chem_S2_taue[1]] + [NL_chem_S2_taue[2]] + [NL_chem_S2_taue[3]] + [NL_chem_S2_taue[4]] + [NL_chem_S2_taue[5]] + [NL_chem_S2_taue[6]] + [NL_chem_S2_taue[7]] + [NL_chem_S2_taue[8]] +
                         [NLb_chem_S2_taue[0]] + [NLb_chem_S2_taue[1]] + [NLb_chem_S2_taue[2]] + [NLb_chem_S2_taue[3]] + [NLb_chem_S2_taue[4]] + [NLb_chem_S2_taue[5]] + [NLb_chem_S2_taue[6]] + [NLb_chem_S2_taue[7]] + [NLb_chem_S2_taue[8]] +
                         [0]*3 + 
                         [NL_chem_S2_taumu[0]] + [NL_chem_S2_taumu[1]] + [NL_chem_S2_taumu[2]] + [NL_chem_S2_taumu[3]] + [NL_chem_S2_taumu[4]] + [NL_chem_S2_taumu[5]] + [NL_chem_S2_taumu[6]] + [NL_chem_S2_taumu[7]] + [NL_chem_S2_taumu[8]] +
                         [NLb_chem_S2_taumu[0]] + [NLb_chem_S2_taumu[1]] + [NLb_chem_S2_taumu[2]] + [NLb_chem_S2_taumu[3]] + [NLb_chem_S2_taumu[4]] + [NLb_chem_S2_taumu[5]] + [NLb_chem_S2_taumu[6]] + [NLb_chem_S2_taumu[7]] + [NLb_chem_S2_taumu[8]] +
                         [0]*3 + 
                         [NL_chem_S2_tautau[0]] + [NL_chem_S2_tautau[1]] + [NL_chem_S2_tautau[2]] + [NL_chem_S2_tautau[3]] + [NL_chem_S2_tautau[4]] + [NL_chem_S2_tautau[5]] + [NL_chem_S2_tautau[6]] + [NL_chem_S2_tautau[7]] + [NL_chem_S2_tautau[8]] +
                         [NLb_chem_S2_tautau[0]] + [NLb_chem_S2_tautau[1]] + [NLb_chem_S2_tautau[2]] + [NLb_chem_S2_tautau[3]] + [NLb_chem_S2_tautau[4]] + [NLb_chem_S2_tautau[5]] + [NLb_chem_S2_tautau[6]] + [NLb_chem_S2_tautau[7]] + [NLb_chem_S2_tautau[8]] +
                         [0]*3,
                         dtype=np.complex128).reshape(21,21)

    Bvec_G2_fact2.append(NLbtau_G2 + np.transpose(NLbtau_G2))
    Bvec_S2_fact2.append(NLbtau_S2 + np.transpose(NLbtau_S2))

    #print('Coefficients computed in', time.time() - start_time, 's')

    return A_G0, A_S0, A_G0_fact2, A_S0_fact2, A_G1, A_S1, A_G1_fact2, A_S1_fact2, A_ham_0, A_ham_LNC, A_ham_LNV, A_L, Bvec_G2, Bvec_S2, Bvec_G2_fact2, Bvec_S2_fact2

#-------------------------------------------------------------------------------------#
#             Quantum Kinetic Equations for ARS Leptogenesis linearised               #
#-------------------------------------------------------------------------------------#
import time

@njit
def set_rates(z, Mav, T,set_rates_kn2):
    # print("wow")

    Tew = 131.7
    x = Tew/T

    h0      = interp_one_over_y0(z)#np.pi**2/(18*zeta(3))
    hLNC    = interp_hp(z)#np.pi**2/(144*zeta(3))#
    #the M dependence is moved to the definition of the interaction Hamiltonian
    hLNV    = interp_hm(z)/(Mav**2)#The dependence over 1/T^2 is kept inside hLNV\

    #Relativistic valeus: 
    G0_rel = 0.013296432287963416#G0_M_fun(M[0,0], x) #9.15e-3
    G1_rel = 0.006737226958120181e-2#7.29e-3
    S0_rel = 0.0432819999972177/T**2 
    S1_rel = 2.5e-2/T**2

    G0      = G0_M_fun(Mav, x)
    
    G1      = G1_M_fun(Mav, x) * 2/(z**2 * set_rates_kn2)

    #the M dependence is moved to the various coefficients of the DMEs
    S0      = S0_M_fun(Mav, x) /(Mav**2) #The dependence over 1/T^2 is kept inside S0
    S1      = S1_M_fun(Mav, x) * 2/(z**2 * set_rates_kn2) /(Mav**2) #The dependence over 1/T^2 is kept inside S1

    #The non-linear terms here do not include the non-relativistic corrections
    G2      = -2.19e-3 
    S2      = -1.651e-2 /T**2

    return h0, hLNC, hLNV, G0, G1, S0, S1, G2, S2

@njit
def lin_RHS(x, y, M1, Tew, M0, Lambda,
            A_G0, A_S0, A_G0_fact2, A_S0_fact2, A_G1, A_S1, A_G1_fact2, A_S1_fact2, A_ham_0, A_ham_LNC, A_ham_LNV, A_L,
            C_col_precomp,
            Bvec_G2, Bvec_S2, Bvec_G2_fact2, Bvec_S2_fact2,kn1,kn2,I21):

    
    #If the splitting is small, Mav \simeq M1
    Mav = M1
    # print(Bvec_G2[1].dtype,Bvec_S2[1].dtype,Bvec_G2_fact2[1].dtype,Bvec_S2_fact2[1].dtype,y.dtype,A.dtype,R.dtype,C_col.dtype )

    # lambda x: Bvec_G2_fact2[x]=Bvec_G2_fact2[x].astype(np.complex128) for x in range(len(Bvec_G2_fact2))


    # Calculate T
    T       = Tew / x

    z       = Mav*x/Tew

    h0, hLNC, hLNV, G0, G1, S0, S1, G2, S2 = set_rates(z, Mav, T,kn2)

    # This is expensive so cache it
    #fdyneq = f_DYNeq(M2, x, Tew, gss)
    #fyneq = f_YNeq(M2, x, Tew, gss)
    temp_z = Mav*x/Tew
    zz  = temp_z.real
    if kn2 == 0:
        fact    = (Mav/Tew)  * 8*zz/(15+8*zz)
        fact2   = 9*1.20205/(2*np.pi**2) * (zz**2 * kn2/2)
    else:
        fact    = (Mav/Tew)   * (kn1/kn2)
        fact2   = (9*1.20205/(2*np.pi**2)) * (zz**2 * kn2/2)


    # Define matrix of the linearised equations
    A = fact*A_L +  (M0/Tew) *(h0 * A_ham_0 * x*x/(Tew*Tew) + hLNC * A_ham_LNC + hLNV * A_ham_LNV +
                               A_G0 * G0 + A_S0 * S0 +
                               A_G0_fact2 * G0*fact2 + A_S0_fact2 * S0*fact2 + A_G1*G1 + A_S1*S1 + A_G1_fact2*G1*fact2 + A_S1_fact2*S1*fact2)
    #The factors of 2 in front of A_G0_fact2 and A_S0_fact2 and the factors of 0.5 in front of A_G1 and A_S1 are inserted by hand
    #to match with the equations given in arxiv:2407.10560


    C_col = C_col_precomp

    #Non-linear terms, add them as 0.5*NL@y in the eqtns
    n = Bvec_G2.shape[0]
    m = y.shape[0]

    NL = np.zeros((n, m), dtype=np.complex128)

    for i in range(n):
        v1 = y @ Bvec_G2[i]
        v2 = y @ Bvec_S2[i]
        v3 = y @ Bvec_G2_fact2[i]
        v4 = y @ Bvec_S2_fact2[i]

        NL[i] = (
            G2 * (M0/Tew) * v1
            + S2 * (M0/Tew) * v2
            + G2 * (M0/Tew) * fact2 * v3
            + S2 * (M0/Tew) * fact2 * v4
        )



    #Regulator
    R = np.linalg.inv(-A*x/ Lambda + I21)
    A_regulated = A @ R

    #Equations are written for RN-1 and RNb-1
    eqtns = (A_regulated@y + R@C_col*fact)# + 0.5*NL@y)

    return np.real(eqtns)

@njit
def Jac(x, M1, Tew, M0, Lambda,
            A_G0, A_S0, A_G0_fact2, A_S0_fact2, A_G1, A_S1, A_G1_fact2, A_S1_fact2, A_ham_0, A_ham_LNC, A_ham_LNV, A_L,
            C_col_precomp,
            kn1,kn2,I21):

    #If the splitting is small, Mav \simeq M1
    Mav = M1

    # Calculate T
    T       = Tew / x

    z       = Mav*x/Tew

    h0, hLNC, hLNV, G0, G1, S0, S1, G2, S2 = set_rates(z, Mav, T,kn2)


    # This is expensive so cache it
    #fdyneq = f_DYNeq(M2, x, Tew, gss)
    #fyneq = f_YNeq(M2, x, Tew, gss)
    temp_z = Mav*x/Tew
    zz  = temp_z.real
    if kn2 == 0:
        fact    = (Mav/Tew)  * 8*zz/(15+8*zz)
        fact2   = 9*1.20205/(2*np.pi**2) * (zz**2 * kn2/2)
    else:
        fact    = (Mav/Tew)   * (kn1/kn2)
        fact2   = (9*1.20205/(2*np.pi**2)) * (zz**2 * kn2/2)


    # Define matrix of the linearised equations
    A = fact*A_L +  (M0/Tew) *(h0 * A_ham_0 * x*x/(Tew*Tew) + hLNC * A_ham_LNC + hLNV * A_ham_LNV +
                               A_G0 * G0 + A_S0 * S0 +
                               A_G0_fact2 * G0*fact2 + A_S0_fact2 * S0*fact2 + A_G1*G1 + A_S1*S1 + A_G1_fact2*G1*fact2 + A_S1_fact2*S1*fact2)
    #The factors of 2 in front of A_G0_fact2 and A_S0_fact2 and the factors of 0.5 in front of A_G1 and A_S1 are inserted by hand
    #to match with the equations given in arxiv:2407.10560

    C_col = C_col_precomp

    #Regulator
    R = np.linalg.inv(-A*x/ Lambda + I21)
    A_regulated = A @ R


    return  np.real(A_regulated), np.real(R@C_col*fact)


class EtaB_ARS_3RHN(ulysses.ULSBase):

    def shortname(self): return "BEARS_3RHN"

    def flavourindices(self): return [7,8,9]

    def flavourlabels(self): return [r"$\mu_{\Delta_e}$",  r"$\mu_{\Delta_\mu}$", r"$\mu_{\Delta_\tau}$"]

    def extendedindices(self): return [1,2,3,4,5,6]

    def extendedlabels(self): return [r"$N_{N_1}$", r"$N_{N_2}$", r"$N_{N_3}$",r"$|N_{N_1}-N_{\overline{N}_1}|$", r"$|N_{N_2}-N_{\overline{N}_2}|$", r"$|N_{N_3}-N_{\overline{N}_3}|$"]    
    

    def RHS(self, x, y, M1, Tew, M0, Lambda,
            A_G0, A_S0, A_G0_fact2, A_S0_fact2, A_G1, A_S1, A_G1_fact2, A_S1_fact2, A_ham_0, A_ham_LNC, A_ham_LNV, A_L,
            C_col_precomp,
            Bvec_G2, Bvec_S2, Bvec_G2_fact2, Bvec_S2_fact2,kn1,kn2,I21):

        return lin_RHS(x, y, M1, Tew, M0, Lambda,
                       A_G0, A_S0, A_G0_fact2, A_S0_fact2, A_G1, A_S1, A_G1_fact2, A_S1_fact2, A_ham_0, A_ham_LNC, A_ham_LNV, A_L,
                       C_col_precomp,
                       Bvec_G2, Bvec_S2, Bvec_G2_fact2, Bvec_S2_fact2,kn1,kn2,I21)
    

    def Jac_for_RHS(self, x,kn1,kn2,I21,A_G0, A_S0, A_G0_fact2, A_S0_fact2, A_G1, A_S1, A_G1_fact2, A_S1_fact2, A_ham_0, A_ham_LNC, A_ham_LNV, A_L,
                    C_col_precomp, Lambda = 1e3):
        # Global constants (masses in GeV)
        Tew         = 131.7
        gss         = 106.75
        M0          = 7.112582895088419e+17

        return Jac(x, self.M1, Tew, M0, Lambda,
                       A_G0, A_S0, A_G0_fact2, A_S0_fact2, A_G1, A_S1, A_G1_fact2, A_S1_fact2, A_ham_0, A_ham_LNC, A_ham_LNV, A_L,
                       C_col_precomp,
                       kn1,kn2,I21)
    
    @property
    def EtaB(self):

        # Global constants (masses in GeV)
        Tew         = 131.7
        gss         = 106.75
        M0          = 7.112582895088419e+17
        zeta3       = zeta(3.)

        # Physical constants
        mp          = 1.672621898e-24  # Proton mass in g
        ngamma      = 410.7  # Present photon number density in cm^-3
        rhoc        = 1.87840e-29  # Critical density of the Universe in h^2 g cm^-3
        gstaro      = 43 / 11  # Entropic effective degrees of freedom at present
        ToYb        = 45 * zeta3 / (gstaro * np.pi**4)
        ToOmegab    = mp * ngamma / rhoc
        self.Lambda = 1e3 if self.Lambda is None else self.Lambda
        self.evolname = r"$T_{\rm{ew}}/T$"

        # Model parameters
        Fmat        = self.h
        dMval_21    = self.M2 - self.M1
        dMval_31    = self.M3 - self.M1

        #initial conditions for the Gell-Mann coefficients for RN-1 and RNb-1 and the chemical potentials
        y0          = np.array([self.initial_abundance-1.] + [0.] * 8 + [self.initial_abundance-1.] + [0.]*8 + [0, 0, 0], dtype=np.complex128)

        # Mass matrices and differences
        M_mat       = np.diag([self.M1, self.M1 + dMval_21, self.M1 + dMval_31])
        deltaM2_21  = 2 * self.M1 * dMval_21 + dMval_21**2
        deltaM2_31  = 2 * self.M1 * dMval_31 + dMval_31**2

        # When N3 is decoupled (all Yukawa couplings to N3 are zero):
        # (1) The N3 Hamiltonian term deltaM2_31 ≈ M3^2 creates extreme numerical
        #     stiffness; zero it so the regulator matrix stays well-conditioned.
        # (2) The A_L and C_col terms must be projected onto the 2-neutrino subspace
        #     so the expanding-universe source contributes exactly zero to d(rN33)/dx.
        n3_decoupled = np.linalg.norm(Fmat[:, 2]) < 1e-17
        print 
        if n3_decoupled:
            deltaM2_31 = 0.0
            Fmat[:, 2] = 0.0

        Dm2_mat     = np.diag([0, deltaM2_21, deltaM2_31])

        # Chi matrix
        chi_mat     = 1./711. * np.array([[257., 20., 20.], [20., 257., 20.], [20., 20., 257.]], dtype=np.complex128)

        #Compute here just once the coefficients of the linearised equations
        A_G0, A_S0, A_G0_fact2, A_S0_fact2, A_G1, A_S1, A_G1_fact2, A_S1_fact2, A_ham_0, A_ham_LNC, A_ham_LNV, A_L, Bvec_G2, Bvec_S2, Bvec_G2_fact2, Bvec_S2_fact2 = compute_linearised_coefficients(Fmat, M_mat, chi_mat, Dm2_mat, n3_decoupled=n3_decoupled)

        # Precompute the C_col source vector.  When N3 is decoupled use the
        # projected identity diag(1,1,0) so the source is consistent with A_L.
        if n3_decoupled:
            C_col_precomp = np.array(GellMann_coefficients_decoupled(np.identity(3)) + GellMann_coefficients_decoupled(np.identity(3)) + [0.,0.,0.], dtype=np.complex128)
        else:
            C_col_precomp = np.array(GellMann_coefficients(np.identity(3)) + GellMann_coefficients(np.identity(3)) + [0.,0.,0.], dtype=np.complex128)

        # Vectors for calculations
        Lvec        = np.zeros((24, 3), dtype=np.complex128)
        Rvec        = np.zeros((24, 3), dtype=np.complex128)
        acr         = np.zeros((24, 3), dtype=np.complex128)

        Bvec_G2 = np.array(Bvec_G2, dtype=np.complex128)
        Bvec_S2 = np.array(Bvec_S2, dtype=np.complex128)
        Bvec_G2_fact2 = np.array(Bvec_G2_fact2, dtype=np.complex128)
        Bvec_S2_fact2 = np.array(Bvec_S2_fact2, dtype=np.complex128)
        I21 = np.identity(21)

        ys          = solve_ivp(lambda x, q: self.RHS(x, q, self.M1, Tew, M0, self.Lambda,
                                                       A_G0, A_S0, A_G0_fact2, A_S0_fact2, A_G1, A_S1, A_G1_fact2, A_S1_fact2, A_ham_0, A_ham_LNC, A_ham_LNV, A_L,
                                                       C_col_precomp,
                                                       Bvec_G2, Bvec_S2, Bvec_G2_fact2, Bvec_S2_fact2,kn(1,self.M1*x/Tew),kn(2,self.M1*x/Tew),I21),
                                                       [1e-6, 1], y0, method='BDF', jac = lambda x, q: self.Jac_for_RHS(x,kn(1,self.M1*x/Tew),kn(2,self.M1*x/Tew),I21, A_G0, A_S0, A_G0_fact2, A_S0_fact2, A_G1, A_S1, A_G1_fact2, A_S1_fact2, A_ham_0, A_ham_LNC, A_ham_LNV, A_L,
                                                                                                         C_col_precomp, Lambda = self.Lambda)[0], atol=1e-13, rtol = 1e-6, max_step = 0.1)
        
        #ys          = odeint(lambda x, t: self.RHS(t, x, self.M1, Tew, M0, self.Lambda,
        #                                               A_G0, A_S0, A_G0_fact2, A_S0_fact2, A_G1, A_S1, A_G1_fact2, A_S1_fact2, A_ham_0, A_ham_LNC, A_ham_LNV, A_L, Bvec_G2, Bvec_S2, Bvec_G2_fact2, Bvec_S2_fact2),
        #                                               y0, np.logspace(-6, 0, 100), #, Dfun = lambda t, x: self.Jac_for_RHS(t, Lambda = self.Lambda)[0],
        #                                               atol=1e-13, rtol = 1e-6, mxstep = 10**6, full_output=True)[0]
        
        t, muD1, muD2, muD3             = [ys.t, ys.y[18], ys.y[19], ys.y[20]]
        self.zs = t
        #t = np.logspace(-6, 0, 100)
        #muD1, muD2, muD3 = [ys[:, 18], ys[:, 19], ys[:, 20]]
        
        # Calculate YB solutions for e, mu, and tau components
        YB_sol_e, YB_sol_mu, YB_sol_tau = map(lambda x: np.real(x[-1]), [muD1, muD2, muD3])
        YB_sol                          = YB_sol_e + YB_sol_mu + YB_sol_tau

        # Conversion factor
        cf          = (28/79) * np.pi**2 / (27 * 6 * zeta(3))

        # Convert YB solutions to YB values
        YBe, YBmu, YBtau, YB            = map(lambda sol: f_convertYBLtoYB(f_convertmutoY(sol, gss)), [YB_sol_e, YB_sol_mu, YB_sol_tau, YB_sol])

        # Calculate etaB values
        etaBe, etaBmu, etaBtau, etaB    = map(lambda yb: yb / ToYb, [YBe, YBmu, YBtau, YB])

        # Calculate Neq values
        Neq         = [(3/16) * (self.M1 * x / Tew)**2 * kn(2, self.M1 * x / Tew) for x in t]

        #Calculate heavy neutrino abundances
        L = gellmann_matrices()
        rN1_minus_id = np.sum([ys.y[j]*l[0,0] for j, l in enumerate(L)], axis = 0)
        rN2_minus_id = np.sum([ys.y[j]*l[1,1] for j, l in enumerate(L)], axis = 0)
        rN3_minus_id = np.sum([ys.y[j]*l[2,2] for j, l in enumerate(L)], axis = 0)

        rN12         = np.sum([ys.y[j]*l[0,1] for j, l in enumerate(L)], axis = 0)
        rN13         = np.sum([ys.y[j]*l[0,2] for j, l in enumerate(L)], axis = 0)
        rN23         = np.sum([ys.y[j]*l[1,2] for j, l in enumerate(L)], axis = 0)

        rNbar1_minus_id = np.sum([ys.y[j+9]*l[0,0] for j, l in enumerate(L)], axis = 0)
        rNbar2_minus_id = np.sum([ys.y[j+9]*l[1,1] for j, l in enumerate(L)], axis = 0)
        rNbar3_minus_id = np.sum([ys.y[j+9]*l[2,2] for j, l in enumerate(L)], axis = 0)

        rNbar12         = np.sum([ys.y[j+9]*l[0,1] for j, l in enumerate(L)], axis = 0)
        rNbar13         = np.sum([ys.y[j+9]*l[0,2] for j, l in enumerate(L)], axis = 0)
        rNbar23         = np.sum([ys.y[j+9]*l[1,2] for j, l in enumerate(L)], axis = 0)

        #rN1_minus_id = np.sum([ys[:, j]*l[0,0] for j, l in enumerate(L)], axis = 0)
        #rN2_minus_id = np.sum([ys[:, j]*l[1,1] for j, l in enumerate(L)], axis = 0)
        #rN3_minus_id = np.sum([ys[:, j]*l[2,2] for j, l in enumerate(L)], axis = 0)

        #Plot haestetics     
        lw              = 0.9
        x1, x2          = -6, 0
        y1, y2          = -13, 1
        xticks, xminor  = ticks_log(x1, x2)
        yticks, yminor  = ticks_log(y1, y2)
        prop_cycle      = plt.rcParams['axes.prop_cycle']
        colors          = prop_cycle.by_key()['color']

        if self.plot == True:
            fig, ax = plt.subplots()
            ax.axhline(6.1e-10, color='grey', linewidth=0.6)
            if self._zcut < 1:
                ax.axvline(self._zcut, color='grey', linewidth=0.6, linestyle='-.')

            # Plotting various components
            ax.plot(t, np.abs((rN1_minus_id + 1) * Neq), label=r"$N_{N_1}$", linewidth=lw, color = colors[4], ls = '-')
            ax.plot(t, np.abs((rN2_minus_id + 1) * Neq), label=r"$N_{N_2}$", linewidth=lw, color = colors[5], ls = '-.')
            ax.plot(t, np.abs((rN3_minus_id + 1) * Neq), label=r"$N_{N_3}$", linewidth=lw, color = colors[6], ls = ':')
            ax.plot(t, np.real(rN1_minus_id - rNbar1_minus_id) * Neq, label=r"$|N_{N_1}-N_{\overline{N}_1}|$", linewidth=lw, color = 'k', ls = '-')
            ax.plot(t, np.real(rN2_minus_id - rNbar2_minus_id) * Neq, label=r"$|N_{N_2}-N_{\overline{N}_2}|$", linewidth=lw, color = 'k', ls = '-.')
            ax.plot(t, np.real(rN3_minus_id - rNbar3_minus_id) * Neq, label=r"$|N_{N_3}-N_{\overline{N}_3}|$", linewidth=lw, color = 'k', ls = ':')
            #ax.plot(t, np.real(rN13) * Neq, label=r"Re($N_{N_1N_3}$)", linewidth=lw, color = 'gold', ls = '-')
            #ax.plot(t, np.real(rN23) * Neq, label=r"Re($N_{N_2N_3}$)", linewidth=lw, color = 'gold', ls = '-.')
            #ax.plot(t, np.real(rN13-rN23) * Neq, label=r"Re($N_{N_1N_3}-N_{N_2 N_3}$)", linewidth=lw, color = 'gold', ls = '-')

            ax.plot(t, np.real(muD1), label=r"$\mu_{\Delta_e}$", linewidth=lw, color = colors[0])
            ax.plot(t, np.real(muD2), label=r"$\mu_{\Delta_\mu}$", linewidth=lw, color = colors[1])
            ax.plot(t, np.real(muD3), label=r"$\mu_{\Delta_\tau}$", linewidth=lw, color = colors[2])
            ax.plot(t, cf * np.real(muD1 + muD2 + muD3), label=r"$\eta_{B}$", linewidth=lw, color = colors[3])
            ax.plot(t, Neq, linewidth=lw, color='k', linestyle=':')

            # Plotting negative components with dashed lines
            ax.plot(t, -np.real(muD1), linestyle='--', color=colors[0], linewidth=lw)
            ax.plot(t, -np.real(muD2), linestyle='--', color=colors[1], linewidth=lw)
            ax.plot(t, -np.real(muD3), linestyle='--', color=colors[2], linewidth=lw)
            ax.plot(t, -cf * np.real(muD1 + muD2 + muD3), linestyle='--', linewidth=lw, color=colors[3])
            ax.plot(t, -np.real(rN1_minus_id - rNbar1_minus_id) * Neq, linewidth=lw, color = 'k', linestyle='-')
            ax.plot(t, -np.real(rN2_minus_id - rNbar2_minus_id) * Neq, linewidth=lw, color = 'k', linestyle='-.')
            ax.plot(t, -np.real(rN3_minus_id - rNbar3_minus_id) * Neq, linewidth=lw, color = 'k', linestyle=':')

            # Setting scales, ticks, and labels
            ax.set_xscale('log')
            ax.set_yscale('log')
            ax.set_xticks(xticks, minor=xminor)
            ax.set_yticks(yticks, minor=yminor)
            ax.set_xlim([10**x1, 10**x2])
            ax.set_ylim([10**y1, 10**y2])
            ax.axvline(1, lw=0.8, ls=':', color='k')
            ax.xaxis.set_ticks_position('both')
            ax.yaxis.set_ticks_position('both')
            ax.tick_params(direction='in', which='both', labelsize=13, width=0.6)
            ax.set_xlabel(r"$T_{\rm{ew}}/T$", fontsize=13)
            ax.set_ylabel(r"$|\eta_B|, \, |\mu_{\Delta_\alpha}|,\, |N_{N_j}-N_{\overline{N}_j}|,\,N_{N_j}$", fontsize=13)
            ax.legend(loc='lower left', shadow=True, ncol=1, prop={'size': 12})

            ax.text(1.5e-6, 8.5e-1, r'$N_N^{\rm{eq}}$', fontsize=11, color='k')

            ax.grid(which='both', linestyle='-', linewidth='0.6', color='grey', alpha=0.2)

            if self.save_plot == True:
                
                fig.savefig(self.path, dpi=300, bbox_inches='tight')
                print(f"Plot saved to {self.path}")

            
            # Show plot
            plt.show()
            plt.close()

        evol_joint = np.array([t,np.abs((rN1_minus_id + 1) * Neq),np.abs((rN2_minus_id + 1) * Neq), np.abs((rN3_minus_id + 1) * Neq), np.real(rN1_minus_id - rNbar1_minus_id) * Neq, 
                               np.real(rN2_minus_id - rNbar2_minus_id) * Neq, np.real(rN3_minus_id - rNbar3_minus_id) * Neq, np.real(muD1), np.real(muD2), np.real(muD3)])
        
        evol_joint = np.transpose(evol_joint)

        self.setEvolDataARS(evol_joint)
                
        return etaB


#---------------------------------------------#
#       Integrated Rates  from Juraj          #
#---------------------------------------------#

import os
import pandas as pd
from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator
import re


# Function to extract mass values from file names
def extract_mass(file_name):
    # Find all sequences of digits in the file name
    num = float(re.findall(r'\d+\.?\d*', file_name)[0])
    # Adjust numbers based on the presence of "TeV" or "MeV"
    if "TeV" in file_name:
        num = num * 1e3
    elif "MeV" in file_name:
        num = num * 1e-3
    return num


data_dir = os.path.dirname(ulysses.__file__)
    
PATHrates_nonrel = os.path.join(data_dir, './data/ARS_rates_nonrel/')

# List all files in the directory
file_names = [f for f in os.listdir(PATHrates_nonrel) if os.path.isfile(os.path.join(PATHrates_nonrel, f))]
# List to store DataFrames
df_list = []

# Process each file
for file_name in file_names:
    # Extract mass value from the file name
    M = extract_mass(file_name)
    
    # Read the data from the file into a DataFrame
    df = pd.read_csv(PATHrates_nonrel + file_name, sep=',', header=None)
    df.columns = ['T', 'gamma1_p', 'gamma1_m', 'gamma_p', 'gamma_m', 'hp', 'hm', 'av_inv_k0', 'dYdzeta']
    
    # Add new columns and perform calculations
    df['M'] = M

    df['z'] = df['M'] / df['T']
    df['x'] = 131.7 / df['T']
    # print(df['x'])
    for c in ['gamma1_p', 'gamma1_m', 'gamma_p', 'gamma_m', 'hp', 'hm']:
        df[c] = df[c] / df['T']
    df['av_inv_y0'] = df['av_inv_k0'] * df['T']
    
    # Reorder columns
    df = df[['M', 'x', 'T', 'z', 'gamma1_p', 'gamma1_m', 'gamma_p', 'gamma_m', 'hp', 'hm', 'av_inv_y0', 'dYdzeta']]
        
    # Append the DataFrame to the list
    df_list.append(df)

# Concatenate all DataFrames into a single DataFrame
combined_df = pd.concat(df_list, ignore_index=True)
combined_df = combined_df.sort_values(by=['M', 'x']).reset_index(drop=True)

#Interpolate Juraj's rates
# Ms = np.log10(combined_df['M'].values)
# Xs = np.log10(combined_df['x'].values)
# G0tab = combined_df['gamma_p'].values
# G1tab = combined_df['gamma1_p'].values
# S0tab = combined_df['gamma_m'].values
# S1tab = combined_df['gamma1_m'].values
# # print(len(Ms))


#Faster interpolation

@njit
def _interp_core(logM, logx,
                 logM_grid,
                 x_flat, v_flat,
                 offsets, lengths):


    nM = len(logM_grid)

    # ------------------------
    # MASS HANDLING
    # ------------------------

    if logM <= logM_grid[0]:
        i = 0
        use_nearest_M = True
    elif logM >= logM_grid[nM-1]:
        i = nM - 1
        use_nearest_M = True
    else:
        i = np.searchsorted(logM_grid, logM) - 1
        use_nearest_M = False

    # ------------------------
    # If outside mass range → full nearest
    # ------------------------

    if use_nearest_M:

        off = offsets[i]
        length = lengths[i]

        xi = x_flat[off:off+length]
        vi = v_flat[off:off+length]

        j = np.searchsorted(xi, logx)

        if j == 0:
            return vi[0]
        elif j >= length:
            return vi[length-1]
        else:
            if abs(logx - xi[j-1]) < abs(logx - xi[j]):
                return vi[j-1]
            else:
                return vi[j]

    # ------------------------
    # Otherwise do bilinear
    # ------------------------

    off_i = offsets[i]
    off_ip1 = offsets[i+1]

    len_i = lengths[i]
    len_ip1 = lengths[i+1]

    xi = x_flat[off_i:off_i+len_i]
    vi = v_flat[off_i:off_i+len_i]

    xi2 = x_flat[off_ip1:off_ip1+len_ip1]
    vi2 = v_flat[off_ip1:off_ip1+len_ip1]

    # ---- check x bounds for both masses
    if (logx < xi[0] or logx > xi[len_i-1] or
        logx < xi2[0] or logx > xi2[len_ip1-1]):

        # nearest in full 2D sense (structured version)

        # nearest mass
        if abs(logM - logM_grid[i]) < abs(logM - logM_grid[i+1]):
            idxM = i
        else:
            idxM = i + 1

        off = offsets[idxM]
        length = lengths[idxM]

        xi = x_flat[off:off+length]
        vi = v_flat[off:off+length]

        j = np.searchsorted(xi, logx)

        if j == 0:
            return vi[0]
        elif j >= length:
            return vi[length-1]
        else:
            if abs(logx - xi[j-1]) < abs(logx - xi[j]):
                return vi[j-1]
            else:
                return vi[j]

    # ---- regular bilinear interpolation

    j = np.searchsorted(xi, logx) - 1
    if j < 0:
        j = 0
    elif j >= len_i - 1:
        j = len_i - 2

    tx = (logx - xi[j]) / (xi[j+1] - xi[j])
    val_i = vi[j] * (1 - tx) + vi[j+1] * tx

    j2 = np.searchsorted(xi2, logx) - 1
    if j2 < 0:
        j2 = 0
    elif j2 >= len_ip1 - 1:
        j2 = len_ip1 - 2

    tx2 = (logx - xi2[j2]) / (xi2[j2+1] - xi2[j2])
    val_ip1 = vi2[j2] * (1 - tx2) + vi2[j2+1] * tx2

    tM = (logM - logM_grid[i]) / (logM_grid[i+1] - logM_grid[i])

    return val_i * (1 - tM) + val_ip1 * tM



def make_interpolator(unique_M, x_arrays, value_arrays):

    logM_grid = np.log10(np.asarray(unique_M))

    lengths = np.array([len(x) for x in x_arrays], dtype=np.int64)
    offsets = np.zeros(len(lengths), dtype=np.int64)
    offsets[1:] = np.cumsum(lengths)[:-1]

    x_flat = np.concatenate([np.log10(x) for x in x_arrays])
    v_flat = np.concatenate(value_arrays)

    return logM_grid, x_flat, v_flat, offsets, lengths



Ms = np.sort(combined_df['M'].unique())

Xs = []
G0tab = []
G1tab = []
S0tab = []
S1tab = []

for M in Ms:
    sub = combined_df[combined_df['M'] == M]
    Xs.append(sub['x'].values)
    G0tab.append(sub['gamma_p'].values)
    G1tab.append(sub['gamma1_p'].values)
    S0tab.append(sub['gamma_m'].values)
    S1tab.append(sub['gamma1_m'].values)



# Create the interpolation function


# points = np.column_stack((Ms, Xs))
# G0_interpolator = LinearNDInterpolator(points, G0tab)
# G0_extrapolator = NearestNDInterpolator(points, G0tab)
# G1_interpolator = LinearNDInterpolator(points, G1tab)
# G1_extrapolator = NearestNDInterpolator(points, G1tab)
# S0_interpolator = LinearNDInterpolator(points, S0tab)
# S0_extrapolator = NearestNDInterpolator(points, S0tab)
# S1_interpolator = LinearNDInterpolator(points, S1tab)
# S1_extrapolator = NearestNDInterpolator(points, S1tab)

G0_logM, G0_xflat, G0_vflat, G0_offsets, G0_lengths = make_interpolator(Ms, Xs, G0tab)
G1_logM, G1_xflat, G1_vflat, G1_offsets, G1_lengths = make_interpolator(Ms, Xs, G1tab)
S0_logM, S0_xflat, S0_vflat, S0_offsets, S0_lengths = make_interpolator(Ms, Xs, S0tab)
S1_logM, S1_xflat, S1_vflat, S1_offsets, S1_lengths = make_interpolator(Ms, Xs, S1tab)



# Define the function to return interpolated G0
@njit
def G0_M_fun(M, x):
    log_M = np.log10(M)
    log_x = np.log10(x)

    interpolated_value = _interp_core(log_M,log_x,G0_logM, G0_xflat, G0_vflat, G0_offsets, G0_lengths)

    return interpolated_value

@njit
def G1_M_fun(M, x):
    log_M = np.log10(M)
    log_x = np.log10(x)
    interpolated_value = _interp_core(log_M,log_x,G1_logM, G1_xflat, G1_vflat, G1_offsets, G1_lengths)

    return interpolated_value

@njit
def S0_M_fun(M, x):
    log_M = np.log10(M)
    log_x = np.log10(x)
    interpolated_value = _interp_core(log_M,log_x, S0_logM, S0_xflat, S0_vflat, S0_offsets, S0_lengths)

    return interpolated_value

@njit
def S1_M_fun(M, x):
    log_M = np.log10(M)
    log_x = np.log10(x)
    interpolated_value = _interp_core(log_M,log_x,S1_logM, S1_xflat, S1_vflat, S1_offsets, S1_lengths)

    return interpolated_value


    """
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
        

        M_f  = os.path.join(data_dir, "./data/Log_M.txt")
        M2_f  = os.path.join(data_dir, "./data/Log_M.txt")
        x_f  = os.path.join(data_dir, "./data/Log_x.txt")
        G0M_f = os.path.join(data_dir, "./data/MyG0log_massdep.txt")
        G1M_f = os.path.join(data_dir, "./data/MyG1log_massdep.txt")
        G2M_f = os.path.join(data_dir, "./data/MyG2log_massdep.txt")
        S0M_f = os.path.join(data_dir, "./data/MyS0log_massdep.txt")
        S1M_f = os.path.join(data_dir, "./data/MyS1log_massdep.txt")
        S2M_f = os.path.join(data_dir, "./data/MyS2log_massdep.txt")
        
        #S0M_f = os.path.join(data_dir, "./data/s0log_massdep.txt")
        #S1M_f = os.path.join(data_dir, "./data/s1log_massdep.txt")
        #S2M_f = os.path.join(data_dir, "./data/s2log_massdep.txt")

        Mtab   = np.loadtxt(M_f, skiprows=0)
        M2tab  = np.loadtxt(M2_f, skiprows=0)
        xtab   = np.loadtxt(x_f, skiprows=0)
        G0Mtab = np.loadtxt(G0M_f, skiprows=0).reshape((len(Mtab), len(xtab)))
        G1Mtab = np.loadtxt(G1M_f, skiprows=0).reshape((len(Mtab), len(xtab)))
        G2Mtab = np.loadtxt(G2M_f, skiprows=0)
        S0Mtab = np.loadtxt(S0M_f, skiprows=0).reshape((len(Mtab), len(xtab)))
        S1Mtab = np.loadtxt(S1M_f, skiprows=0).reshape((len(Mtab), len(xtab)))
        S2Mtab = np.loadtxt(S2M_f, skiprows=0)
        
        self.G0MInt_ = RectBivariateSpline(Mtab, xtab, G0Mtab) # 2-D Interpolation
        self.G1MInt_ = RectBivariateSpline(Mtab, xtab, G1Mtab) # 2-D Interpolation
        self.G2MInt_ = RectBivariateSpline(M2tab, xtab, G2Mtab) # 2-D Interpolation

        self.S0MInt_ = RectBivariateSpline(Mtab, xtab, S0Mtab) # 2-D Interpolation
        self.S1MInt_ = RectBivariateSpline(Mtab, xtab, S1Mtab) # 2-D Interpolation
        self.S2MInt_ = RectBivariateSpline(M2tab, xtab, S2Mtab) # 2-D Interpolation

    
    def G0_fun(self,x): return interpolate.splev(math.log(x), self.G0Int_, der=0)

    def G1_fun(self,x): return interpolate.splev(math.log(x), self.G1Int_, der=0)

    def G2_fun(self,x): return interpolate.splev(math.log(x), self.G2Int_, der=0)

    

    def G0_M_fun(self, M, x):

        if M < 0.1:
            M = 0.1
        elif M > 5e5:
            M = 5e5

        return self.G0MInt_(np.log(M), np.log(x))[0,0]

    def G1_M_fun(self,M, x):

        if M < 0.1:
            M = 0.1
        elif M > 5e5:
            M = 5e5

        return self.G1MInt_(np.log(M), np.log(x))[0,0]

    def G2_M_fun(self,M, x):

        if M < 0.1:
            M = 0.1
        elif M > 100:
            M = 100

        return self.G2MInt_(np.log(M), np.log(x))[0,0]


    def S0_M_fun(self,M, x):

        if M < 0.1:
            M = 0.1
        elif M > 5e5:
            M = 5e5

        return self.S0MInt_(np.log(M), np.log(x))[0,0]

    def S1_M_fun(self,M, x):

        if M < 0.1:
            M = 0.1
        elif M > 5e5:
            M = 5e5

        return self.S1MInt_(np.log(M), np.log(x))[0,0]

    def S2_M_fun(self,M, x):

        if M < 0.1:
            M = 0.1
        elif M > 100.:
            M = 100.

        return self.S2MInt_(np.log(M), np.log(x))[0,0]

    def RHS(self, x, y, Fmat, M1, deltaM_21, deltaM_31, Tew, gss, M0, M_mat, Dm2_mat, chi_mat, Lvec, Rvec, acr, Lambda, A_G0, A_S0, A_G0_fact2, A_S0_fact2, A_G1, A_S1, A_G1_fact2, A_S1_fact2, A_ham_0, A_ham_LNC, A_ham_LNV, A_L):

        funcs = [self.G0_M_fun, self.G1_M_fun, self.G2_M_fun, self.S0_M_fun, self.S1_M_fun, self.S2_M_fun]

        return lin_RHS(x, y, Fmat, M1, deltaM_21, deltaM_31, Tew, gss, M0, M_mat, Dm2_mat, chi_mat, Lvec, Rvec, acr, funcs, Lambda, A_G0, A_S0, A_G0_fact2, A_S0_fact2, A_G1, A_S1, A_G1_fact2, A_S1_fact2, A_ham_0, A_ham_LNC, A_ham_LNV, A_L, use_interpolation=True)
    """
