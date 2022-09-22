###################################################################################################
#                                                                                                 #
#                               Schwarzschild and Kerr BHs Library                                #
#                                                                                                 #
#         Authors: Andrew Cheek, Lucien Heurtier, Yuber F. Perez-Gonzalez, Jessica Turner         #
#                                   Based on: arXiv:2107.00013                                    #
#                                                                                                 #
###################################################################################################
import ulysses
import numpy as np
import pandas as pd
from scipy import interpolate
import scipy.integrate as integrate
from scipy.integrate import quad, ode, solve_ivp, odeint
from scipy.optimize import root
from scipy.special import zeta, kn, spherical_jn, jv
from scipy.interpolate import RectBivariateSpline
import mpmath
from mpmath import polylog

from numpy import sqrt, log, exp, log10, pi, logspace, linspace, seterr, min, max, append
from numpy import loadtxt, zeros, floor, ceil, unique, sort, cbrt, concatenate, delete, real

from collections import OrderedDict
olderr = np.seterr(all='ignore')

# Particle masses, in GeV

mW   = 80.379
mZ   = 91.1876
mH   = 125.18
me   = 0.5109989461e-3
mmu  = 105.6583745e-3
mtau = 1.77686
mu   = 2.2e-3
md   = 4.6e-3
ms   = 95e-3
mc   = 1.275
mb   = 4.18
mt   = 173.1
mg   = 0.6      # Ficticious gluon mass ---> indicates the QCD phase transition, following PRD41(1990)3052

# Degrees of freedom of the SM ---> Before the EW phase transition

gW  = 4.        # W
gZ  = 2.        # Z
gH  = 4.        # Higgs
gp  = 2.        # photon
gg  = 2.        # graviton
ggl = 16.       # gluons
gl  = 2.*2.     # leptons
gq  = 2.*2.*3.  # quarks
gnu = 2.        # LH neutrino

gf = 3.*gnu + 3.*gl + 6.*gq   # Total number of SM fermion dofs 
gs = gH                       # Total number of SM scalar dofs
gv = gW + gZ + gp + gg + ggl  # Total number of SM vector dofs

# Constants

c     = 299792.458       # in km/s
gamma = np.sqrt(3.)**-3. # Collapse factor
GCF   = 6.70883e-39      # Gravitational constant in GeV^-2
mPL   = GCF**-0.5        # Planck mass in GeV
v     = 174              # Higgs vev
csp   = 0.35443          # sphaleron conversion factor
GF    = 1.1663787e-5     # Fermi constant in GeV^-2

# Conversion factors

GeV_in_g     = 1.782661907e-24  # 1 GeV in g
Mpc_in_cm    = 3.085677581e24   # 1 Mpc in cm

cm_in_invkeV = 5.067730938543699e7       # 1 cm in keV^-1
year_in_s    = 3.168808781402895e-8      # 1 year in s
GeV_in_invs  = cm_in_invkeV * c * 1.e11  # 1 GeV in s^-1

MPL   = mPL * GeV_in_g        # Planck mass in g
kappa = mPL**4 * GeV_in_g**3  # Evaporation constant in g^3 * GeV -- from PRD41(1990)3052

# BH Temperature in GeV

def TBH(M, astar):

    M_GeV = M/GeV_in_g
    
    return (1./(4.*np.pi*GCF*M_GeV))*(np.sqrt(abs(1. - astar**2))/(1. + np.sqrt(abs(1. - astar**2)))) # M in g

#-------------------------------------------------------------------------------------------------------------------------------------#
#                                                  Momentum Integrated Rate for Kerr BHs                                              #
#-------------------------------------------------------------------------------------------------------------------------------------#


def Gamma_S(M, ast, m):# Scalar, in GeV

    GM = GCF * (M/GeV_in_g) # in GeV^-1

    g1, g2, g3, g4, g5, g6, g7, g8 = [0.3916590458886117, 0.15958502152092527, 1.9234512633493943, 0.00017703648060188705,
                                      11.260851723128718, 26.429307340860653, 28.982647658008826, 11.91716289208516]
    
    TKBH = TBH(M, ast)
    
    hs = 10.**(g1 - g2*ast + g3*ast**2 + (g4*ast**2)/(-1.025 + ast)**2 - g5*ast**3 + g6*ast**4 - g7*ast**5 + g8*ast**6)

    if m > 0.:
        
        if ast <= 1.e-5:

            B, C, nu = [7.50218, -2.94379, 0.420819]

            z = m/TKBH

        else:

            a0, a1, a2, a3, a4, a5 = [0.90067, -0.287571, 2.06242, -6.031, 4.3491, -0.000203751]
            b0, b1, b2, b3, b4, b5 = [7.68412, -1.1945, 3.42557, -19.2999, 11.6408, -0.00076169]
            c0, c1, c2, c3, c4, c5 = [-0.438955, -0.570661, 2.3257, -0.981604, -0.97489, 0.000357955]

            B  = 10.**(a0 + a1*ast + a2*ast**2 + a3*ast**3 + a4*ast**4 + (a5*ast**2)/(ast-1.025)**2)
            C  = b0 + b1*ast + b2*ast**2 + b3*ast**3 + b4*ast**4 + (b5*ast**2)/(ast-1.025)**2
            nu = 10.**(c0 + c1*ast + c2*ast**2 + c3*ast**3 + c4*ast**4 + (c5*ast**2)/(ast-1.025)**2)

            z = GM * m

        In = hs * (1. - (1. + np.exp(-B * np.log10(abs(z)) - C))**(-nu)) # DM emission rate including greybody factors
        
    else:
        
        In = hs
    
    return  (27/(1024. * np.pi**4 * GM)) * In


def Gamma_F(M, ast, m):# Fermion

    GM = GCF * (M/GeV_in_g) # in GeV^-1

    TKBH = TBH(M, ast)
    
    hf = 10.**(-0.0469657 - 0.0112247*ast + 0.790198*ast**2 + (0.000174539*ast**2)/(-1.025 + ast)**2 -
               0.80843*ast**3 + 0.535612*ast**4)

    if m > 0.:
        
        if ast <= 1.e-5:

            B, C, nu = [12.3573, -8.74364, 0.304554]

            z = m/TKBH

        else:

            a0, a1, a2, a3, a4, a5 = [1.02775, 0.251744, -1.91938, 3.71237, -2.57412, 0.000208832]
            b0, b1, b2, b3, b4, b5 = [8.64208, 0.646042, -11.8172, 18.2938, -13.6375, -0.00100791]
            c0, c1, c2, c3, c4, c5 = [-0.494509, -0.169789, 1.60784, -2.20497, 1.33059, -0.000398046]

            B  = 10.**(a0 + a1*ast + a2*ast**2 + a3*ast**3 + a4*ast**4 + (a5*ast**2)/(ast-1.025)**2)
            C  = b0 + b1*ast + b2*ast**2 + b3*ast**3 + b4*ast**4 + (b5*ast**2)/(ast-1.025)**2
            nu = 10.**(c0 + c1*ast + c2*ast**2 + c3*ast**3 + c4*ast**4 + (c5*ast**2)/(ast-1.025)**2)

            z = GM * m

        In = hf * (1. - (1. + np.exp(-B * np.log10(abs(z)) - C))**(-nu)) # DM emission rate including greybody factors
                        
    else:
        
        In = hf
    
    return  2. * (27/(1024. * np.pi**4 * GM)) * In


def Gamma_V(M, ast, m):# Vector

    GM = GCF * (M/GeV_in_g) # in GeV^-1

    TKBH = TBH(M, ast)
    
    hv = 10.**(-0.556741 - 0.194255*ast + 3.85504*ast**2 + (0.000192461*ast**2)/(-1.025 + ast)**2 - 4.67469*ast**3 + 2.56983*ast**4)

    if m > 0.:
        
        if ast <= 1.e-5:

            B, C, nu = [13.465, -9.8134, 0.304932]

            z = m/TKBH

        else:

            a0, a1, a2, a3, a4, a5 = [1.12764, -0.00975927, 0.172785, -0.152433, -0.362799, -0.000059454]
            b0, b1, b2, b3, b4, b5 = [8.99996, -1.07481, 3.30235, -12.0108, 3.96199, -0.000601034]
            c0, c1, c2, c3, c4, c5 = [-0.529893, -0.0088484, 0.0297965, -0.677536, 1.14249, 0.0000377463]

            B  = 10.**(a0 + a1*ast + a2*ast**2 + a3*ast**3 + a4*ast**4 + (a5*ast**2)/(ast-1.025)**2)
            C  = b0 + b1*ast + b2*ast**2 + b3*ast**3 + b4*ast**4 + (b5*ast**2)/(ast-1.025)**2
            nu = 10.**(c0 + c1*ast + c2*ast**2 + c3*ast**3 + c4*ast**4 + (c5*ast**2)/(ast-1.025)**2)

            z = GM * m

        In = hv * (1. - (1. + np.exp(-B * np.log10(abs(z)) - C))**(-nu)) # DM emission rate including greybody factors
                                
    else:
        
        In = hv
    
    return  3. * (27/(1024. * np.pi**4 * GM)) * In

def Gamma_G(M, ast, m):# Spin 2

    GM = GCF * (M/GeV_in_g) # in GeV^-1

    TKBH = TBH(M, ast)
    
    hg = 10.**(-1.59185 - 0.860299*ast + 11.1725*ast**2 + (0.00020413*ast**2)/(-1.025 + ast)**2 - 13.5798*ast**3 + 6.74455*ast**4)

    if m > 0.:
        
        if ast <= 1.e-5:

            B, C, nu = [22.325, -21.2326, 0.12076]

            z = m/TKBH

        else:

            a0, a1, a2, a3, a4, a5 = [1.27648, 0.234641, -0.967098, 1.59757, -1.05198, -0.000101089]
            b0, b1, b2, b3, b4, b5 = [9.06958, 2.94595, -25.4011, 39.1391, -24.1023, -0.00125056]
            c0, c1, c2, c3, c4, c5 = [-0.646339, -0.0594674, 0.0142806, -0.0170435, 0.253691, 0.0000985018]

            B  = 10.**(a0 + a1*ast + a2*ast**2 + a3*ast**3 + a4*ast**4 + (a5*ast**2)/(ast-1.025)**2)
            C  = b0 + b1*ast + b2*ast**2 + b3*ast**3 + b4*ast**4 + (b5*ast**2)/(ast-1.025)**2
            nu = 10.**(c0 + c1*ast + c2*ast**2 + c3*ast**3 + c4*ast**4 + (c5*ast**2)/(ast-1.025)**2)

            z = GM * m

        In = hg * (1. - (1. + np.exp(-B * np.log10(abs(z)) - C))**(-nu)) # DM emission rate including greybody factors
                                        
    else:
        
        In = hg
    
    return  5. * (27/(1024. * np.pi**4 * GM)) * In

def Gamma_GO(M, ast, m):# Geometric optics limit

    GM = GCF * (M/GeV_in_g) # in GeV^-1

    TKBH = TBH(M, ast)

    zBH = m/TKBH

    In = - zBH * polylog(2, -np.exp(-zBH)) - polylog(3, -np.exp(-zBH))# DM emission rate including greybody factors
    
    return  2. * (27/(1024. * np.pi**4 * GM)) * In

#-------------------------------------------------------------------------------------------------------------------------------------#
#                                         Total f functions ---> related to the mass rate, dM/dt                                      #
#                                                   Counting SM dofs + Dark Radiation                                                 #
#-------------------------------------------------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------#
#                     f and g fitted functions                      #
#            taken from PRD14(1976)3260 + 9801044[gr-qc]            #
#-------------------------------------------------------------------#

def fs(astar): return (-4.11848 - 0.418271*astar + 2.58436*astar**2 + (0.0000849804*astar**2)/(astar - 1.025)**2
                       -5.76425*astar**3 + 4.01629*astar**4)

def ff(astar): return (-4.38503 - 0.0168371*astar + 1.1853*astar**2 + (0.000261808*astar**2)/(astar - 1.025)**2
                       - 1.21264*astar**3 + 0.803417*astar**4)
    
def fv(astar): return (-4.7754443 + 0.0629953*astar + 3.15186*astar**2 + (0.000199368*astar**2)/(astar - 1.025)**2
                       - 3.52188*astar**3 + 1.99382*astar**4)

def fg(astar): return (-5.71339 + 0.550866*astar + 7.5178*astar**2 + (0.000330209*astar**2)/(astar - 1.025)**2
                       - 9.50996*astar**3 + 5.34276*astar**4)

#--------------------------------------------------------------------------------#
#              Our interpolated forms including the particle's mass              #
#--------------------------------------------------------------------------------#

# Scalar

def phi_s(M, ast, m):

    GM = GCF * (M/GeV_in_g) # in GeV^-1
    
    TKBH = TBH(M, ast)

    f0 = 10.**fs(ast)

    if m > 0.:
        
        if ast <= 1.e-5:

            B, C, nu = [7.79984, -3.80742, 0.48848]

            z = m/TKBH

        else:

            a0, a1, a2, a3, a4, a5 = [0.862565, 1.06174, -6.40438, 10.3813, -5.12991, 0.000108425]
            b0, b1, b2, b3, b4, b5 = [7.02688, 2.99615, -25.1091, 31.049, -14.5699, -0.00145935]
            c0, c1, c2, c3, c4, c5 = [-0.280799, -1.87129, 11.5739, -20.6905, 11.1745, -0.000245393]

            B  = 10.**(a0 + a1*ast + a2*ast**2 + a3*ast**3 + a4*ast**4 + (a5*ast**2)/(ast-1.025)**2)
            C  = b0 + b1*ast + b2*ast**2 + b3*ast**3 + b4*ast**4 + (b5*ast**2)/(ast-1.025)**2
            nu = 10.**(c0 + c1*ast + c2*ast**2 + c3*ast**3 + c4*ast**4 + (c5*ast**2)/(ast-1.025)**2)

            z = GM * m

        In = f0 * (1. - (1. + np.exp(-B * np.log10(abs(z)) - C))**(-nu))   
        
    else:
        
        In = f0

    return In

# Fermion

def phi_f(M, ast, m):

    GM = GCF * (M/GeV_in_g) # in GeV^-1
    
    TKBH = TBH(M, ast)

    f12  = 10.**ff(ast)

    if m > 0.:
        
        if ast <= 1.e-5:
    
            B, C, nu = [13.0496, -9.91178, 0.3292]
    
            z = m/TKBH

        else:

            a0, a1, a2, a3, a4, a5 = [1.05952, -0.549081, 3.00671, -6.6903, 3.97656, 0.000107068]
            b0, b1, b2, b3, b4, b5 = [8.19678, 2.07543, -19.0044, 22.0031, -12.3065, -0.00101776]
            c0, c1, c2, c3, c4, c5 = [-0.474275, 1.15307, -6.77883, 15.0118, -9.39296, -1.40672e-6]

            B  = 10.**(a0 + a1*ast + a2*ast**2 + a3*ast**3 + a4*ast**4 + (a5*ast**2)/(ast-1.025)**2)
            C  = b0 + b1*ast + b2*ast**2 + b3*ast**3 + b4*ast**4 + (b5*ast**2)/(ast-1.025)**2
            nu = 10.**(c0 + c1*ast + c2*ast**2 + c3*ast**3 + c4*ast**4 + (c5*ast**2)/(ast-1.025)**2)

            z = GM * m

        In = f12 * (1. - (1. + np.exp(-B * np.log10(abs(z)) - C))**(-nu))
        
    else:
        
        In = f12

    return In

# Vector

def phi_v(M, ast, m):

    GM = GCF * (M/GeV_in_g) # in GeV^-1
    
    TKBH = TBH(M, ast)

    f1 = 10.**fv(ast)

    if m > 0.:
        
        if ast <= 1.e-5:

            B, C, nu = [14.0361, -10.7138, 0.307206]

            z = m/TKBH

        else:

            a0, a1, a2, a3, a4, a5 = [1.14383, -0.00734714, -0.0599119, 0.136389, -0.575339, 0.000107856]
            b0, b1, b2, b3, b4, b5 = [8.88373, -1.45578, 2.35702, -14.6024, 7.09857, -0.00110925]
            c0, c1, c2, c3, c4, c5 = [-0.517544, 0.0317647, 0.0742752, -0.408531, 0.943227, -0.000266956]

            B  = 10.**(a0 + a1*ast + a2*ast**2 + a3*ast**3 + a4*ast**4 + (a5*ast**2)/(ast-1.025)**2)
            C  = b0 + b1*ast + b2*ast**2 + b3*ast**3 + b4*ast**4 + (b5*ast**2)/(ast-1.025)**2
            nu = 10.**(c0 + c1*ast + c2*ast**2 + c3*ast**3 + c4*ast**4 + (c5*ast**2)/(ast-1.025)**2)

            z = GM * m

        In = f1 * (1. - (1. + np.exp(-B * np.log10(abs(z)) - C))**(-nu))
        
    else:
        
        In = f1
    
    return In

# Tensor - spin2

def phi_g(M, ast, m):

    GM = GCF * (M/GeV_in_g) # in GeV^-1
    
    TKBH = TBH(M, ast)

    f2 = 10.**fg(ast)

    if m > 0.:
        
        if ast <= 1.e-5:

            B, C, nu = [21.50941, -20.5135, 0.173423]

            z = m/TKBH

        else:

            a0, a1, a2, a3, a4, a5 = [1.29608, 0.0664447, -0.0755711, 0.0676736, -0.302979, -0.000065566]
            b0, b1, b2, b3, b4, b5 = [9.06618, 1.29693, -16.2621, 21.5632, -14.46, -0.000515066]
            c0, c1, c2, c3, c4, c5 = [-0.644022, -0.250355, 1.23193, -2.46744, 1.85603, 0.0000199225]

            B  = 10.**(a0 + a1*ast + a2*ast**2 + a3*ast**3 + a4*ast**4 + (a5*ast**2)/(ast-1.025)**2)
            C  = b0 + b1*ast + b2*ast**2 + b3*ast**3 + b4*ast**4 + (b5*ast**2)/(ast-1.025)**2
            nu = 10.**(c0 + c1*ast + c2*ast**2 + c3*ast**3 + c4*ast**4 + (c5*ast**2)/(ast-1.025)**2)

            z = GM * m

        In = f2 * (1. - (1. + np.exp(-B * np.log10(abs(z)) - C))**(-nu))
        
    else:
        
        In = f2
    
    return In

#------------------------------------------#
#              SM Contribution             #
#------------------------------------------#

def fSM(M, ast):
    
    # Contribution from each particle --> We do not include the Graviton contribution here

    fgr =  0.
    fp  =  gp * phi_v(M, ast, 0.)  # Photon
    fgl = ggl * phi_v(M, ast, 0.6) # Gluon
    fW  =  gW * phi_v(M, ast, mW)  # W
    fZ  =  gZ * phi_v(M, ast, mZ)  # Z
    fH  =  gH * phi_s(M, ast, mH)  # Higgs

    fnu = 3. * gnu * phi_f(M, ast, 0.) # Active neutrinos
    
    fl  = gl * (phi_f(M, ast, me) + phi_f(M, ast, mmu) + phi_f(M, ast, mtau))  # Charged leptons
    
    fq  = gq * (phi_f(M, ast, mu) + phi_f(M, ast, md) + phi_f(M, ast, ms) +
                phi_f(M, ast, mc) + phi_f(M, ast, mb) + phi_f(M, ast, mt))    # Quarks

    
    return fgr + fp + fgl + fW + fZ + fH + fnu + fl + fq

# RH neutrino contribution

def fRH(M, ast, mrh): return gnu * phi_f(M, ast, mrh)

# DM contribution

def fDM(M, ast, mdm): return gnu * phi_f(M, ast, mdm)
def fX(M, ast, mX):   return 3. * phi_v(M, ast, mX)

# Dark Radiation

def fDR(M, ast, s): # Depending on the particle's spin

    if s == 0.:
        f = phi_s(M, ast, 0.)
    elif s == 0.5:
        f = 2.*phi_f(M, ast, 0.)
    elif s == 1.:
        f = 2.*phi_v(M, ast, 0.)
    elif s == 2.:
        f = 2.*phi_g(M, ast, 0.)
        
    return f

#-------------------------------------------------------------------------------------------------------------------------------------#
#                                Total g functions ---> related to the angular momentum rate, da_*/dt                                 #
#                                                   Counting SM dofs + Dark Radiation                                                 #
#-------------------------------------------------------------------------------------------------------------------------------------#

def gs(astar): return (-4.04521 - 0.251753*astar + 2.31411*astar**2 + (0.0000672319*astar**2)/(astar - 1.025)**2
                       - 3.47359*astar**3 + 2.20082*astar**4)
    
def gf(astar): return (-3.51098 - 0.0545547*astar + 0.811723*astar**2 + (0.000157456*astar**2)/(astar - 1.025)**2
                       - 1.62005*astar**3 + 1.23504*astar**4)

def gv(astar): return (-3.63003 + 0.262898*astar + 0.0846463*astar**2 + (0.000180678*astar**2)/(astar - 1.025)**2
                       + 0.614102*astar**3 + 0.00796626*astar**4)

def gG(astar): return (-4.26363 - 0.0849301*astar + 5.63412*astar**2 + (0.000279216*astar**2)/(astar - 1.025)**2
                       - 6.57782*astar**3 + 3.89969*astar**4)

#--------------------------------------------------------------------------------#
#              Our interpolated forms including the particle's mass              #
#--------------------------------------------------------------------------------#

# Scalar

def gam_s(M, ast, m):

    GM = GCF * (M/GeV_in_g) # in GeV^-1
    
    g0 = 10.**gs(ast)
    
    if m > 0.:
              
        z = GM * m # Dimensionless parameter -- gravitational coupling GMm

        a0, a1, a2, a3, a4, a5 = [1.14795, -0.188218, 0.957975, -2.36396, 1.16129, 0.000139667]
        b0, b1, b2, b3, b4, b5 = [8.02772, 0.807776, -10.1362, 3.33735, -0.490688, -0.0014108]
        c0, c1, c2, c3, c4, c5 = [-0.552166, 0.471878, -2.55669, 5.12681, -2.65038, -0.000218156]
    
        B  = 10.**(a0 + a1*ast + a2*ast**2 + a3*ast**3 + a4*ast**4 + (a5*ast**2)/(ast-1.025)**2)
        C  = b0 + b1*ast + b2*ast**2 + b3*ast**3 + b4*ast**4 + (b5*ast**2)/(ast-1.025)**2
        nu = 10.**(c0 + c1*ast + c2*ast**2 + c3*ast**3 + c4*ast**4 + (c5*ast**2)/(ast-1.025)**2)

        In = g0 * (1. - (1. + np.exp(-B * np.log10(abs(z)) - C))**(-nu))
        
    else:
        
        In = g0
    
    return In

# Fermion

def gam_f(M, ast, m):

    GM = GCF * (M/GeV_in_g) # in GeV^-1

    g12 = 10.**gf(ast)

    if m > 0.:
              
        z = GM * m # Dimensionless parameter -- gravitational coupling GMm

        a0, a1, a2, a3, a4, a5 = [1.00516, -0.30863, 1.35856, -3.82659, 2.60927, 0.0000676727]
        b0, b1, b2, b3, b4, b5 = [7.50414, 1.55438, -18.3266, 23.8847, -13.7473, -0.000844954]
        c0, c1, c2, c3, c4, c5 = [-0.432342, 0.319384, -1.43802, 5.66093, -4.65577, 0.0000408488]
    
        B  = 10.**(a0 + a1*ast + a2*ast**2 + a3*ast**3 + a4*ast**4 + (a5*ast**2)/(ast-1.025)**2)
        C  = b0 + b1*ast + b2*ast**2 + b3*ast**3 + b4*ast**4 + (b5*ast**2)/(ast-1.025)**2
        nu = 10.**(c0 + c1*ast + c2*ast**2 + c3*ast**3 + c4*ast**4 + (c5*ast**2)/(ast-1.025)**2)

        In = g12 * (1. - (1. + np.exp(-B * np.log10(abs(z)) - C))**(-nu))        
                
    else:
        
        In = g12

    return In

# Vector

def gam_v(M, ast, m):

    GM = GCF * (M/GeV_in_g) # in GeV^-1

    g1 = 10.**gv(ast)

    if m > 0.:
              
        z = GM * m # Dimensionless parameter -- gravitational coupling GMm

        a0, a1, a2, a3, a4, a5 = [1.12718, 0.000404064, -0.0201877, -0.26672, -0.185709, 0.0000739435]
        b0, b1, b2, b3, b4, b5 = [8.61971, -0.00904373, -3.23898, -7.78225, 4.85716, -0.00109291]
        c0, c1, c2, c3, c4, c5 = [-0.521117, 0.0923845, -0.471924, 0.81145, 0.163192, -0.00021713]
    
        B  = 10.**(a0 + a1*ast + a2*ast**2 + a3*ast**3 + a4*ast**4 + (a5*ast**2)/(ast-1.025)**2)
        C  = b0 + b1*ast + b2*ast**2 + b3*ast**3 + b4*ast**4 + (b5*ast**2)/(ast-1.025)**2
        nu = 10.**(c0 + c1*ast + c2*ast**2 + c3*ast**3 + c4*ast**4 + (c5*ast**2)/(ast-1.025)**2)

        In = g1 * (1. - (1. + np.exp(-B * np.log10(abs(z)) - C))**(-nu))
                        
    else:
        
        In = g1

    return In

# Tensor - spin2

def gam_g(M, ast, m):

    GM = GCF * (M/GeV_in_g) # in GeV^-1

    g2 = 10.**gG(ast)

    if m > 0.:
              
        z = GM * m # Dimensionless parameter -- gravitational coupling GMm

        a0, a1, a2, a3, a4, a5 = [1.27466, 0.481525, -2.48397, 4.3418, -2.68619, 0.000005166558]
        b0, b1, b2, b3, b4, b5 = [8.91067, 2.7764, -21.9831, 28.0939, -16.192, -0.000732776]
        c0, c1, c2, c3, c4, c5 = [-0.638896, -0.390212, 1.97792, -3.60446, 2.40479, -0.00000106767]
    
        B  = 10.**(a0 + a1*ast + a2*ast**2 + a3*ast**3 + a4*ast**4 + (a5*ast**2)/(ast-1.025)**2)
        C  = b0 + b1*ast + b2*ast**2 + b3*ast**3 + b4*ast**4 + (b5*ast**2)/(ast-1.025)**2
        nu = 10.**(c0 + c1*ast + c2*ast**2 + c3*ast**3 + c4*ast**4 + (c5*ast**2)/(ast-1.025)**2)

        In = g2 * (1. - (1. + np.exp(-B * np.log10(abs(z)) - C))**(-nu))
                                
    else:
        In = g2
        #print(In)
        
    return In

#------------------------------------------#
#              SM Contribution             #
#------------------------------------------#

def gSM(M, ast):

    # Contribution from each particle --> We do not include the Graviton contribution here

    fgr =  0.                      # Graviton
    fp  =  gp * gam_v(M, ast, 0.)  # Photon
    fgl = ggl * gam_v(M, ast, 0.6) # Gluon
    fW  =  gW * gam_v(M, ast, mW)  # W
    fZ  =  gZ * gam_v(M, ast, mZ)  # Z
    fH  =  gH * gam_s(M, ast, mH)  # Higgs

    fnu = 3. * gnu * gam_f(M, ast, 0.)                           # Active neutrinos
    
    fl  = gl * (gam_f(M, ast, me) + gam_f(M, ast, mmu) + gam_f(M, ast, mtau))  # Charged leptons
    
    fq  = gq * (gam_f(M, ast, mu) + gam_f(M, ast, md) + gam_f(M, ast, ms) +
                gam_f(M, ast, mc) + gam_f(M, ast, mb) + gam_f(M, ast, mt))    # Quarks

    return fgr + fp + fgl + fW + fZ + fH + fnu + fl + fq

def gRH(M, ast, mRH): return gnu * gam_f(M, ast, mRH)

# Mediator contribution

def gDM(M, ast, mdm): return gnu * gam_f(M, ast, mdm)

def gX(M, ast, mX): return 3. * gam_v(M, ast, mX)

# Dark Radiation

def gDR(M, ast, s):

    if s == 0.:
        f = gam_s(M, ast, 0.)
    elif s == 0.5:
        f = 2.*gam_f(M, ast, 0.)
    elif s == 1.:
        f = 2.*gam_v(M, ast, 0.)
    elif s == 2.:
        f = 2.*gam_g(M, ast, 0.)
        
    return f

#-----------------------------------------------------------------------------#
#                               PBHs lifetime                                 #
#-----------------------------------------------------------------------------#

def ItauSM(tl, v): # Standard Model + Gravitons
    
    M   = v[0]
    ast = v[1]

    FSM = fSM(M, ast) + gg * phi_g(M, ast, 0.)  # 
    GSM = gSM(M, ast) + gg * gam_g(M, ast, 0.)  #

    dMdtl   = - np.log(10.) * 10.**tl * kappa * FSM * M**-2
    dastdtl = - np.log(10.) * 10.**tl * ast * kappa * M**-3 * (GSM - 2.*FSM)

    return [dMdtl, dastdtl]

def ItauFO(tl, v, mDM): # Freeze Out case
    
    M   = v[0]
    ast = v[1]

    FSM = fSM(M, ast)
    FDM = fDM(M, ast, mDM) # DM evaporation contribution
    FT  = FSM + FDM        # Total Evaporation contribution

    GSM = gSM(M, ast)
    GDM = gDM(M, ast, mDM) # DM evaporation contribution
    GT  = GSM + GDM        # Total Evaporation contribution

    dMdtl   = - np.log(10.) * 10.**tl * kappa * FT * M**-2
    dastdtl = - np.log(10.) * 10.**tl * ast * kappa * M**-3 * (GT - 2.*FT)

    return [dMdtl, dastdtl]

def ItauDR(tl, v, s): # Dark Radiation Case
    
    M   = v[0]
    ast = v[1]

    FSM = fSM(M, ast)
    FDR = fDR(M, ast, s) # DM evaporation contribution
    FT  = FSM + FDR      # Total Evaporation contribution

    GSM = gSM(M, ast)
    GDR = gDR(M, ast, s) # DM evaporation contribution
    GT  = GSM + GDR      # Total Evaporation contribution

    dMdtl   = - np.log(10.) * 10.**tl * kappa * FT * M**-2
    dastdtl = - np.log(10.) * 10.**tl * ast * kappa * M**-3 * (GT - 2.*FT)

    return [dMdtl, dastdtl]

def ItauFI(tl, v, mDM, mX): # Freeze In case (Including mediator)
    
    M   = v[0]
    ast = v[1]

    FSM = fSM(M, ast)
    FDM = fDM(M, ast, mDM) # DM evaporation contribution
    FX  = fX(M, ast, mX)   # Mediator contribution
    FT  = FSM + FDM + FX   # Total Evaporation contribution

    GSM = gSM(M, ast)
    GDM = gDM(M, ast, mDM) # DM evaporation contribution
    GX  = gX(M, ast, mX)  # Mediator contribution
    GT  = GSM + GDM + GX   # Total Evaporation contribution

    dMdtl   = - np.log(10.) * 10.**tl * kappa * FT * M**-2
    dastdtl = - np.log(10.) * 10.**tl * ast * kappa * M**-3 * (GT - 2.*FT)

    return [dMdtl, dastdtl]

def ItauRH(tl, v, M1, M2, M3): # Including 3 RH neutrinos
    
    M   = v[0]
    ast = v[1]

    FSM  = fSM(M, ast)
    FRH1 = fRH(M, ast, M1)           # 1 RH neutrino evaporation contribution
    FRH2 = fRH(M, ast, M2)           # 2 RH neutrino evaporation contribution
    FRH3 = fRH(M, ast, M3)           # 3 RH neutrino evaporation contribution
    FT   = FSM + FRH1 + FRH2 + FRH3  # Total Evaporation contribution

    GSM  = gSM(M, ast)
    GRH1 = gRH(M, ast, M1)           # 1 RH neutrino evaporation contribution
    GRH2 = gRH(M, ast, M2)           # 2 RH neutrino evaporation contribution
    GRH3 = gRH(M, ast, M3)           # 3 RH neutrino evaporation contribution
    GT   = GSM + GRH1 + GRH2 + GRH3  # Total Evaporation contribution

    dMdtl   = - np.log(10.) * 10.**tl * kappa * FT * M**-2
    dastdtl = - np.log(10.) * 10.**tl * ast * kappa * M**-3 * (GT - 2.*FT)

    return [dMdtl, dastdtl]

# Determining the scale fator where PBHs evaporate

def afin(aexp, rPBHi, rRadi, t, ail):

    a = [10.**(aexp[0])]

    ain = 10.**ail # Initial scale factor
    
    A = -ain * rPBHi * np.sqrt(GCF * (ain * rPBHi + rRadi))
    B = a[0] * rPBHi * np.sqrt(GCF * (a[0] * rPBHi + rRadi))
    C = 2. * rRadi * (np.sqrt(GCF*(ain * rPBHi + rRadi)) - np.sqrt(GCF*(a[0]*rPBHi + rRadi)))
    D = GCF * np.sqrt(6.*np.pi) * rPBHi**2
    
    return [A + B + C - D*t]

#-------------------------------------------------------------------------------------------------#
#                                   g*(T) and g*S(T) interpolation                                #
#-------------------------------------------------------------------------------------------------#

import os
data_dir = os.path.dirname(ulysses.__file__)

fg  = os.path.join(data_dir, "etab1BE1Fscalefactor_gstar.txt")
fgS = os.path.join(data_dir, "etab1BE1Fscalefactor_gstarS.txt")

gTab = np.loadtxt(fg, skiprows=0)

Ttab = gTab[:,0]
gtab = gTab[:,1]
tck  = interpolate.splrep(Ttab, gtab, s=0)

def gstar(T): return interpolate.splev(T, tck, der=0)

def dgstardT(T): return interpolate.splev(T, tck, der = 1)

gSTab = np.loadtxt(fgS, skiprows=0)

TStab = gSTab[:,0]
gstab = gSTab[:,1]
tckS  = interpolate.splrep(TStab, gstab, s=0)

def gstarS(T): return interpolate.splev(T, tckS, der = 0)

def dgstarSdT(T): return interpolate.splev(T, tckS, der = 1)
