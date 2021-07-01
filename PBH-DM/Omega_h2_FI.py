###################################################################################################
#                                                                                                 #
#                         Primordial Black Hole + Dark Matter Generation.                         #
#                                     Evaporation + Freeze-In                                     #
#                                                                                                 #
#         Authors: Andrew Cheek, Lucien Heurtier, Yuber F. Perez-Gonzalez, Jessica Turner         #
#                    Based on: arXiv:2107.xxxxx (P1) and  arXiv:2107.xxxxx (P2)                   #
#                                                                                                 #
###################################################################################################

import ulysses
import numpy as np
import math
from odeintw import odeintw
import pandas as pd
from scipy import interpolate
import matplotlib.pyplot as plt
import scipy.integrate as integrate
from scipy.integrate import quad, ode, solve_ivp, odeint
from scipy.optimize import root
from scipy.special import zeta, kn
from scipy.interpolate import interp1d, RectBivariateSpline

from math import sqrt, log, exp, log10, pi, atan

import BHProp as bh #Schwarzschild and Kerr BHs library

from collections import OrderedDict
olderr = np.seterr(all='ignore')


# -------------------- Main Parameters ---------------------------
#
#
#     'model' : [int] parameter to specify the model used. The model implemented here is using the value '1'.
#
#     In model 1 : 
#
#          - 'MDM' : mass of the fermionic dark matter
#
#          - 'mf'  : mass of a light fermion in equilibrium with the SM
#
#          - 'Mmed' : mass of the vector mediator X
#
#          - 'gV'   : coupling of the mediator to SM particles
#
#          - 'gD'   : coupling of the mediator to DM particles
#
#-----------------------------------------------------------------

#--------------------------   Credits  -----------------------------#
#
#      If using this code, please cite:
#
#      - Arxiv:    and    Arxiv:
#
#      - JCAP 12 (2017) 013 • e-Print: 1706.03118 (WDM constraints)
#
#-------------------------------------------------------------------#

#------------------------------------------------------------#
#    Stopping function to reach MBH=Mp during evaporation    #
#------------------------------------------------------------#

def PlanckMass(t, v, Mi):

    if (0.01*Mi > bh.mPL * bh.GeV_in_g): Mst = 0.01*Mi
    else: Mst = bh.mPL * bh.GeV_in_g
    
    return v[0] - Mst # Function to stop the solver if the BH is equal or smaller than the Planck mass

#----------------------------------------------------#
#                Mediator decay width                #
#----------------------------------------------------#

def Gamma_med(paramsDM):
    from math import sqrt,log, exp
    MDM , mf , Mmed , gV , gD , g_DM, model = paramsDM
    
    if(model == 1):
        factor=Mmed/(12*pi)
        if(Mmed>2.*MDM):
            term_DM=gD**2*(1+2*(MDM**2)/(Mmed**2))*sqrt(1-4*(MDM**2)/(Mmed**2))
        else:
            term_DM=0
        if(Mmed>2.*mf):
            term_SM=gV**2*(1+2*(mf**2)/(Mmed**2))*sqrt(1-4*(mf**2)/(Mmed**2))
        else:
            term_SM=0
            
    return factor*(term_DM+term_SM)

#----------------------------------------------------#
#      Annihilation channel amplitude^2              #
#----------------------------------------------------#

def int_Omega_M2(s,paramsDM): # Function including the three processes, DM DM -> SM SM, DM DM, X X

    MDM, mf, Mmed, gV, gD, g_DM, model = paramsDM

    # mediator decay width
    GamX = Gamma_med(paramsDM)
    
    
    #------------------------------------
    #   Different processes
    #------------------------------------
    
    ## DM + DM -> SM + SM  (s-channel)

    factorSM = 16.*np.pi*gV**2*gD**2/3.

    IntSM = factorSM*(s+2*(MDM**2))*(s+2*(mf**2))/((s-Mmed**2)**2+(Mmed**2)*GamX**2)

    ## DM + DM -> DM + DM  (s and t-channel)
    
    X1DM = atan(Mmed/GamX) - atan((-4.*MDM**2+Mmed**2+s)/(Mmed*GamX))
    X2DM = -8.* MDM**2*(Mmed**2+s)+16*MDM**4+Mmed**2*(2.*s+GamX**2)+Mmed**4+s**2
    X3DM = Mmed**2*(Mmed**2+GamX**2)
    
    factorDM = 8*pi*gD**4 
    
    denomDM = (3.* Mmed*GamX *(s-4.*MDM**2)*((s-Mmed**2)**2+Mmed**2*GamX**2))
    
    num1DM = 3.*X1DM*Mmed**2*GamX**2*(8.*s*MDM**2-8*Mmed**2*s-6*Mmed**4+s**2)
               
    num2DM = Mmed*GamX*(12.*MDM**2*(2*s*Mmed**2+Mmed**4-3*s**2) + 12*MDM**4*(Mmed**2*(-log(X2DM/X3DM)+4.)+s*(log(X2DM/X3DM)-6.))
                        -32.*MDM**6+s*(-9.*s*Mmed**2+Mmed**4*(6.*log(X2DM/X3DM)-3.) + 2.*s**2*(-3.*log(X2DM/X3DM)+7.)))          
    
    num3DM = (-3.*X1DM*(s-Mmed**2)**2*(2*(s-2*MDM**2)**2+2.*s*Mmed**2+Mmed**4)
              -3.*Mmed**3*GamX**3*(12.*MDM**2+4.*Mmed**2*log(X2DM/X3DM)+s*(2*log(X2DM/X3DM)-3))+9*X1DM*Mmed**4*GamX**4)

    IntDM = factorDM*(num1DM+num2DM+num3DM)/denomDM

    ## DM + DM -> X + X (t-channel)

    if s <= 1.e15*Mmed:
        X1XX = s - 2.*Mmed**2 - np.sqrt((s-4.*MDM**2)*(s-4.*Mmed**2)) + 1.e-10
    else:
        X1XX = 2.*MDM**2 
    X2XX = s - 2.*Mmed**2 + np.sqrt((s-4.*MDM**2)*(s-4.*Mmed**2)) 
    X3XX = np.sqrt((s-4.*MDM**2)*(s-4.*Mmed**2))
    
    factorXX = 16*pi*gD**4
    
    denomXX = X3XX*(s - 2.*Mmed**2)*(Mmed**4 + MDM**2*(s - 4.*Mmed**2))
    
    num1XX = -(-4.*Mmed**6 + 2.*Mmed**4 * s + s * MDM**2 * (-2.*Mmed**2 + s) + MDM**4*(-8.*Mmed**2 + 4.*s))*X3XX            
    num2XX = -(8.*MDM**6 * (4.*Mmed**2 - s) + Mmed**4*(4.*Mmed**4 + s**2) + 4.*MDM**4 *(6.*Mmed**4 - 6.*Mmed**2*s + s**2)
                + MDM**2 *(-24.*Mmed**6 + 8.*Mmed**4*s - 4.*Mmed**2 * s**2 + s**3))*log(X1XX/X2XX)

    if s > 4.*Mmed**2:  # Including a Heaviside theta th(s - 4mX^2)
        IntXX = factorXX*(num1XX + num2XX)/denomXX
    else:
        IntXX = 1.e-50

    if IntXX != IntXX: print(X1XX, 2.*MDM**2, 2.*Mmed**2)
        
    if(model == 1):
        return [IntSM, IntDM, IntXX]
    
#---------------------------------------------------------------------------------#
#                           <sigma.v>(T1,T2) computation                          #
#---------------------------------------------------------------------------------#

def sigma_v(T1, T2, paramsDM):

    #----------------------- Integrands -------------------------#
    def Int(t, pars):

        T1, T2, MDM, mf, Mmed, gV, gD, g_DM, model = pars
        
        paramsDM = [MDM, mf, Mmed, gV, gD, g_DM, model]
        
        z =  MDM*(T1+T2)/(T1*T2) + (1. - t)/t
        
        s = (T1**2*T2**2*z**2 - MDM**2 * (T1 - T2)**2)/(T1*T2)
        
        A = MDM**2 * (T2**2 - T1**2)/(T1*T2)
        C = s*sqrt(abs(1. - 4.*MDM**2/s))

        Fac = A*(1. + z)*exp(-z)/z + C*kn(1, z)
        
        return [i*Fac/t**2 for i in int_Omega_M2(s,paramsDM)]
    #-----------------------------------------------------------#

    MDM, mf, Mmed, gV, gD, g_DM, model = paramsDM

    pars = [T1, T2, MDM, mf, Mmed, gV, gD, g_DM, model]

    lims = [1.e-10, 1. - 1.e-10, 100]

    integ = Simp(Int, pars, lims).integral_1D()

    den = 8.*MDM**4*kn(2, MDM/T1)*kn(2, MDM/T2)

    return (1./(32.*pi**2))*abs(integ/den)

#---------------------------------------------------------------------------------#
#           <sigma.v>(T) using the Narrow-Width Approximation                     #
#---------------------------------------------------------------------------------#

def sigma_v_NWA(T, paramsDM):
    
    MDM, mf, Mmed, gV, gD, g_DM, model = paramsDM

    GamX = Gamma_med(paramsDM)

    ## DM + DM -> SM + SM
    neq=MDM**2*T/(np.pi**2)*kn(2, MDM/T)
    
    amp_factor = 3.*g_DM**2*gD**2*gV**2
    
    factor = T/(512.*np.pi**5*neq**2)

    IntSM = np.pi * np.sqrt(Mmed**2-4*MDM**2)*(Mmed**2+2*(MDM**2))*(Mmed**2+2*(mf**2))/(Mmed*GamX)*kn(1, Mmed/T)
    
    return amp_factor*factor*IntSM

#---------------------------------------------------------------------------------#
#           Interaction rate in the Narrow-Width Approximation                    #
#---------------------------------------------------------------------------------#


def Gamma_NWA(paramsDM):
    
    MDM, mf, Mmed, gV, gD, g_DM, model = paramsDM

    GamX = Gamma_med(paramsDM)
    
    T=0.785939*Mmed

    ## DM + DM -> SM + SM
    neq=MDM**2*T/(np.pi**2)*kn(2, MDM/T)
    
    amp_factor = 3.*g_DM**2*gD**2*gV**2
    
    factor = T/(512.*np.pi**5*neq**2)

    IntSM = np.pi * np.sqrt(Mmed**2-4*MDM**2)*(Mmed**2+2*(MDM**2))*(Mmed**2+2*(mf**2))/(Mmed*GamX)*kn(1, Mmed/T)
    
    return amp_factor*factor*IntSM*neq

#########################################################
#           Average DM and Mediator Momentum
#########################################################

def p_average_DM(Mi, asi, MDM, tau, Sol_t):

    def Integ_p(t, pars):

        MDM, sol = pars

        M = sol(t)[0]
        a = sol(t)[1]

        return 10.**t * log(10.) * bh.fDM(M, a, MDM)/M**2

    def Integ_n(t, pars):

        MDM, sol = pars

        M = sol(t)[0]
        a = sol(t)[1]

        return 10.**t * log(10.) * bh.Gamma_F(M, a, MDM)

    pars = [MDM, Sol_t]

    integ_p = integrate.quad(Integ_p, -10., tau, args=(pars))
    integ_n = integrate.quad(Integ_n, -10., tau, args=(pars))

    if(integ_n[0]!=0):

        return (bh.kappa * integ_p[0]/bh.GeV_in_g)/integ_n[0]
    else:
        return 0
    
def p_average_med(Mi, asi, MX, tau, Sol_t):

    def Integ_p(t, pars):

        MX, sol = pars

        M = sol(t)[0]
        a = sol(t)[1]

        return 10.**t * log(10.) * bh.fX(M, a, MX)/M**2

    def Integ_n(t, pars):

        MX, sol = pars

        M = sol(t)[0]
        a = sol(t)[1]

        return 10.**t * log(10.) * bh.Gamma_V(M, a, MX)

    pars = [MX, Sol_t]

    integ_p = integrate.quad(Integ_p, -10., tau, args=(pars))
    integ_n = integrate.quad(Integ_n, -10., tau, args=(pars))

    return (bh.kappa * integ_p[0]/bh.GeV_in_g)/integ_n[0]


#----------------------------------------------------------#
#    Warm Dark Matter limit -> True if DM is too warm      #
#----------------------------------------------------------#

def Is_DM_hot(Mi, asi, MDM, tau, Sol_t, Tev, NDMp, NDMbh): #NDMp and NDMbh are the values of the number densities at T=T0
    
    #extract the limits from data file (extracted from 1706.03118 Fig. 6)
    lim_f_WDM_Tab = pd.read_table("./Data/lim_f_WDM_new.dat",  names=['x','lim_f_WDM'])
    
    # x = ( keV / MDM )
    x_tab         = lim_f_WDM_Tab.iloc[:,0]
    lim_f_WDM_tab = lim_f_WDM_Tab.iloc[:,1]
    x_WDM_lim     = interp1d(lim_f_WDM_tab, x_tab)

    T0 = 2.34865e-13  # Temperature today in GeV

    # average momentum from evaporation of DM
    p_average = p_average_DM(Mi, asi, MDM, tau, Sol_t)
    
    # velocity of DM at evaporation
    vbh = p_average/MDM
    
    # velocity of DM today
    vbh0 = T0 / Tev * (bh.gstarS(T0)/bh.gstarS(Tev))**(1./3.) * vbh    
    
    # fraction of DM evaporated as compared to produced from FI/FO
    f_WDM=NDMbh/(NDMp+NDMbh)
    
    # If the fraction is small enough --> no constraint
    if(f_WDM < 0.024031133452245834):
        return False
    
    # If the fraction is large enough --> check constraints
    else:
        x_lim=x_WDM_lim(f_WDM)
        v_lim=3.9 * (1e-8) * (x_lim)**(4./3.)
        
        if(vbh0<=v_lim):
            return False
        else:
            return True
        
#------------------------------------------------------#
#              non-relativistic parameterization       #
#------------------------------------------------------#

def find_gV_gD(sv, Br_DM, mDM, mX, mf): # Br_DM is the branching fraction of X -> DM
    
    if(mDM>=mf):
        
        if(mX>2*mf):
            
            if(mX>2*mDM):
                
                factorD = 12.*Br_DM*mX*np.pi*(mX**2-4.*mDM**2)
                numD=sqrt(mDM*(2.*mf**2+mX**2)*sv/(2.*mDM**2+mX**2))
                A=72.*np.pi*Br_DM*(2.*mDM**2 + mf**2)*mX**2 * np.sqrt( (-mDM**2 + mf**2)*(-4. *mDM**2 + mX**2)/(4.*mf**2 - mX**2) )
                B=-72.*np.pi*Br_DM**2 * (2.*mDM**2 + mf**2)*mX**2 * np.sqrt( (-mDM**2 + mf**2)*(-4.* mDM**2 + mX**2)/(4.*mf**2 - mX**2) )
                C=mDM*(2.*mf**2 + mX**2)*(8.*mDM**4 + 2.*mDM**2*mX**2 - mX**4) * sv
                denD=np.sqrt( A+B+C )
                
                gD=sqrt(factorD*numD/denD)
                
                factorV=12.*(1-Br_DM)*mX*np.pi*(mX**2-4.*mDM**2)**(3./2.)
                numV=sqrt(mDM*(2.*mf**2+mX**2)*(2.*mDM**2+mX**2)*sv/(-4.*mf**2+mX**2)) / (2.*mf**2+mX**2)
                denV=denD
                
                gV=sqrt(factorV*numV/denV)
                
                return [gV,gD]

#########################################################

#----------------------------------------#
#   Diff. Equations before evaporation   #
#----------------------------------------#

def FBEqs(a, v, nphi, paramsDM, GammaX, p_DM, p_X, Br_DM, FO):

    M     = v[0] # PBH mass
    ast   = v[1] # PBH ang mom
    rRAD  = v[2] # Radiation energy density
    rPBH  = v[3] # PBH energy density
    Tp    = v[4] # Temperature
    NDMT  = v[5] # Thermal DM number density
    NDMB  = v[6] # PBH-induced DM number density
    NX    = v[7] # X number density

    NDMC  = v[8] # Thermal DM number density w/o PBH contribution
    NDMH  = v[9] # PBH-induced DM number density w/o thermal contact

    #----------------#
    #   Parameters   #
    #----------------#
    
    mDM, mf, mX, gV, gD, g_DM, model = paramsDM
    
    FSM = bh.fSM(M, ast)      # SM contribution
    FDM = bh.fDM(M, ast, mDM) # DM contribution
    FX  = bh.fX(M, ast, mX)   # Mediator contribution
    FT  = FSM + FDM + FX      # Total Energy contribution

    GSM = bh.gSM(M, ast)      # SM contribution
    GDM = bh.gDM(M, ast, mDM) # DM contribution
    GX  = bh.gX(M, ast, mDM)  # Mediator contribution
    GT  = GSM + GDM + GX      # Total Angular Momentum contribution
    
    H   = np.sqrt(8 * pi * bh.GCF * (rPBH * 10.**(-3*a) + rRAD * 10.**(-4*a))/3.) # Hubble parameter
    Del = 1. + Tp * bh.dgstarSdT(Tp)/(3. * bh.gstarS(Tp)) # Temperature parameter

    Br_SM = 1. - Br_DM # SM Branching ratio

    from ulysses.ulsbase import my_kn2, my_kn1
    
    TH  = bh.TBH(M, ast) # Hawking Temperature
    z   = mDM/Tp       
    zBH = mDM/TH
    
    # thermally averaged decay width of X
    GXt = min([GammaX * my_kn1(zBH)/my_kn2(zBH), 1e5*H])
    p_X_0  = p_X  * 10.**(-a)
    E_X  = sqrt(mX**2  + p_X_0**2)
    
    rX = E_X*NX*nphi # Mediator Energy density
    
    #----------------------------------------------#
    #    Radiation + PBH + Temperature equations   #
    #----------------------------------------------#

    dMda    = - bh.kappa * FT * M**-2/H     # (1 g/M)^2
    dastda  = - ast * bh.kappa * M**-3 * (GT - 2.*FT)/H
    drRADda = - (FSM/FT) * (dMda/M) * 10.**a * rPBH + 2.*(Br_SM*GXt/H) * 10.**a * rX
    drPBHda = + (dMda/M) * rPBH
    dTda    = - (Tp/Del) * (1.0 + (bh.gstarS(Tp)/bh.gstar(Tp))*((FSM/FT)*(dMda/M)* 10.**a * rPBH -
                                                                2.*(Br_SM*GXt/H) * 10.**a * rX)/(4.*rRAD))
    
    #---------------------------------------------#
    #           Temperature averaged <sv>         #
    #---------------------------------------------#

    NDMeq = (10.**(3*a) * mDM**2 * Tp * kn(2,z))/(pi**2)/nphi
    
    svTT = sigma_v_NWA(Tp, paramsDM) # <sigma v>(T_plasma, T_plasma)

    #-----------------------------------------#
    #           Dark Matter Equations         #
    #-----------------------------------------#
    
    dNDMTda = -(NDMT**2 - NDMeq**2)*svTT*nphi/(H*10.**(3.*a))
        
    dNDMBda = 2.*(Br_DM*GXt/H)*NX + (bh.Gamma_F(M, ast, mDM)/H)*(rPBH/(M/bh.GeV_in_g))/nphi

    dNXda  = -NX*GXt/H + (bh.Gamma_V(M, ast, mX)/H)*(rPBH/(M/bh.GeV_in_g))/nphi

    dNDMCda = -(NDMC**2 - NDMeq**2)*svTT*nphi/(H*10.**(3.*a))                               # Thermal Contribution w/o PBH evap

    dNDMHda = (bh.Gamma_F(M, ast, mDM)/H)*(rPBH/(M/bh.GeV_in_g))/nphi + 2.*(Br_DM*GXt/H)*NX # PBH-induced contribution w/o contact
    
    ##########################################################    
    
    dEqsda = [dMda, dastda, drRADda, drPBHda, dTda, dNDMTda, dNDMBda, dNXda, dNDMCda, dNDMHda]

    return [x * log(10.) for x in dEqsda]

#----------------------------------------#
#   Diff. Equations after evaporation    #
#----------------------------------------#

def FBEqs_aBE(a, v, nphi, paramsDM, GammaX, a_evap, T_bh_in, p_DM, p_X, Br_DM, FO):

    rRAD = v[0] # Radiation energy density
    Tp   = v[1] # Temperature
    NDMT = v[2] # Total DM number density
    NDMB = v[3] # PBH-induced DM number density
    NX   = v[4] # X number density

    NDMC = v[5] # Thermal DM number density w/o PBH contribution
    NDMH = v[6] # Thermal DM number density w/o PBH contribution
    
    #----------------#
    #   Parameters   #
    #----------------#
    
    mDM, mf, mX, gV, gD, g_DM, model = paramsDM
    
    H   = sqrt(8 * pi * bh.GCF * (rRAD * 10.**(-4*a))/3.)    # Hubble parameter
    Del = 1. + Tp * bh.dgstarSdT(Tp)/(3. * bh.gstarS(Tp))          # Temperature parameter
            
    z = mDM/Tp
    
    ####
    # Boost factor of the mediator and DM
    ####
    p_X_0  = p_X  * 10.**(a_evap - a)
    p_DM_0 = p_DM * 10.**(a_evap - a)
    
    E_DM = sqrt(mDM**2 + p_DM_0**2)
    E_X  = sqrt(mX**2  + p_X_0**2)
    
    if NX >= 0.: GXt = min([GammaX * mX/E_X, 1e3*H]) # boosted X width (saturated at 10^3 H in order to avoid stiffness)
    else: GXt = 0.
    
    rX = E_X*NX*nphi # Mediator energy density

    Br_SM = 1. - Br_DM # SM Branching ratio
    
    #----------------------------------------#
    #    Radiation + Temperature equations   #
    #----------------------------------------#

    drRADda = 2.*(Br_SM*GXt/H) * 10.**a * rX
    dTda    = - Tp/Del * (1.0 - (bh.gstarS(Tp)/bh.gstar(Tp))*(2.*(Br_SM*GXt/H) * 10.**a * rX/(4.*rRAD)))

    #---------------------------------------------#
    #           Temperature averaged <sv>         #
    #---------------------------------------------#

    #from ulysses.ulsbase import my_kn2, my_kn1

    NDMeq = (10.**(3*a) * mDM**2 * Tp * kn(2, z))/(pi**2)/nphi # Equilibrium DM number density

    TH = E_DM # We take the BH temp to be the energy average of the DM

    svTT = sigma_v_NWA(Tp, paramsDM) # <sigma v>(T_plasma, T_plasma)


    #-----------------------------------------#
    #           Dark Matter Equations         #
    #-----------------------------------------#

    dNDMTda = -(NDMT**2 -  NDMeq**2)*svTT*nphi/(H*10.**(3.*a)) 
        
    dNDMBda = 2*(Br_DM*GXt/H)*NX

    dNXda   = -NX*GXt/H
    
    dNDMCda = -(NDMC**2 -  NDMeq**2)*svTT*nphi/(H*10.**(3.*a))  # Thermal Contribution w/o PBH evap
    
    dNDMHda = 2*(Br_DM*GXt/H)*NX                              # PBH-induced contribution w/o contact
    
    dEqsda  = [drRADda, dTda, dNDMTda, dNDMBda, dNXda, dNDMCda, dNDMHda]

    return [x * log(10.) for x in dEqsda]

#------------------------------------------------------------------------------------------------------------------#
#                                            Input parameters                                                      #
#------------------------------------------------------------------------------------------------------------------#
class FrInPBH:

    def __init__(self, MPBHi, aPBHi, bPBHi, mDM, mX, mf, sv, BR, g_DM, model):

        self.MPBHi  = MPBHi # Log10[M/1g]
        self.aPBHi  = aPBHi # a_star
        self.bPBHi  = bPBHi # Log10[beta']
        self.mDM    = mDM
        self.mX     = mX
        self.mf     = mf
        self.sv     = sv
        self.BR     = BR
        self.g_DM   = g_DM
        self.model  = model
    
#------------------------------------------------------------------------------------------------------------------------------------#
#                                                       Input parameters                                                             #
#------------------------------------------------------------------------------------------------------------------------------------#
    
    def Omegah2(self):
        
        Mi     = 10**(self.MPBHi) # PBH initial Mass in grams
        asi    = self.aPBHi       # PBH initial rotation a_star factor
        bi     = 10**(self.bPBHi) # Initial PBH fraction
        Ti     = ((45./(16.*106.75*(pi*bh.GCF)**3.))**0.25) * sqrt(bh.gamma * bh.GeV_in_g/Mi) # Initial Universe temperature
        rRadi  = (pi**2./30.) * bh.gstar(Ti) * Ti**4  # Initial radiation energy density -- assuming a radiation dominated Universe
        rPBHi  = abs(bi/(sqrt(bh.gamma) -  bi))*rRadi # Initial PBH energy density
        nphi   = (2.*zeta(3)/pi**2)*Ti**3             # Initial photon number density
        
    
        TBHi   = bh.TBH(Mi, asi)  # Initial BH temperature
        
        mDM     = 10**self.mDM    # DM mass in GeV
        mX      = 10**self.mX     # Mediator Mass
        mf      = 10**self.mf     # SM fermion mass
        sv      = 10**self.sv     # <sigma*v>
        BR      = self.BR         # Branching ratio X -> DM
        g_DM    = self.g_DM       # DM d.o.f.
        model   = self.model      # DM Model
                
        # Derive fraction evaporated
        
        FSM_test = bh.fSM(Mi, asi)      # SM contribution
        FDM_test = bh.fDM(Mi, asi, mDM) # DM contribution
        FX_test  = bh.fX(Mi, asi, mX)   # Mediator contribution
        FT_test  = FSM_test + FDM_test + FX_test      # Total Energy contribution
        
        frac_SM=FSM_test/(FT_test)
        print('fraction evaporated into SM = ',frac_SM*100., ' %')
        
        
        if mX < 2.*mDM:
            print("Mediator mass should be bigger than 2*DM mass")
            exit()

        gV, gD  = find_gV_gD(sv, BR, mDM, mX, mf)  # gV, gD couplings obtained from sv and BR

        print('g_V = ', gV)
        print('g_D = ', gD)

        paramsDM=[mDM, mf, mX, gV, gD, g_DM, model]
        
        G_X = Gamma_med(paramsDM)

        #+++++++++++++++++++++++++++++++++++++++++++++++++++#
        #          Table for <sv>(T1,T2) interpolation      #
        #+++++++++++++++++++++++++++++++++++++++++++++++++++#

        Tmin = min(1.e-2*mDM, 0.01*TBHi)

        ns = 50
        nx = 15

        #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
        #                                           Solving the equations                                                   #
        #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

        #--------------------------------------------------------------------------------#
        #         Computing PBH lifetime and scale factor in which BHs evaporate         #
        #--------------------------------------------------------------------------------#
        
        MPL = lambda t, x:PlanckMass(t, x, Mi)
        MPL.terminal  = True
        MPL.direction = -1.

        tau_sol = solve_ivp(fun=lambda t, y: bh.ItauFI(t, y, mDM, mX), t_span = [-10., 40.], y0 = [Mi, asi], 
                            events=MPL, rtol=1.e-10, atol=1.e-20, dense_output=True)

        Sol_t = tau_sol.sol # Solutions for obtaining <p>
        tau = tau_sol.t[-1] # Log10@PBH lifetime in inverse GeV

        if bi > 1.e-19*(1.e9/Mi):
            af = root(bh.afin, [40.], args = (rPBHi, rRadi, 10.**tau, 0.), method='lm', tol=1.e-10) # Scale factor 
            aflog10 = af.x[0]            
        else:
            afw = np.sqrt(1. + 4.*10.**tau*np.sqrt(2.*np.pi*bh.GCF*rRadi/3.))
            aflog10 = np.log10(afw)

        #+++++++++++++++++++++++++++#
        #      Average momentum     #
        #+++++++++++++++++++++++++++#

        p_DM = p_average_DM(Mi, asi, mDM, tau, Sol_t)
        p_X  = p_average_med(Mi, asi, mDM, tau, Sol_t)
        
        #-----------------------------------------#
        #          Before BH evaporation          #
        #-----------------------------------------#
    
        v0 = [Mi, asi, rRadi, rPBHi, Ti, 0., 0., 0., 0., 0.]

        FO = [Ti, 0., 0., 0., 0., 0., 0., 0., 10.] # Temp, a,  <sv's> at DM decoupling, neq, nDM_BH, H

        # solve ODE
        solFBE = solve_ivp(lambda t, z: FBEqs(t, z, nphi, paramsDM,
                                              G_X, p_DM, p_X, BR, FO),
                           [0., 1.25*aflog10], v0, method='BDF', events=MPL, rtol=1.e-6, atol=1.e-10)

        aflog10 = solFBE.t[-1] # We update the value of log(a) at which PBHs evaporate

        
        #-----------------------------------------#
        #           After BH evaporation          #
        #-----------------------------------------#
        
        Tfin = 1.e-2*mDM # Final plasma temp in GeV
        
        azmax = aflog10 + np.log10(np.cbrt(bh.gstarS(solFBE.y[4,-1])/bh.gstarS(Tfin))*(solFBE.y[4,-1]/Tfin))
        afmax = max(aflog10, azmax)

        v0aBE = [solFBE.y[2,-1], solFBE.y[4,-1], solFBE.y[5,-1], solFBE.y[6,-1], solFBE.y[7,-1], solFBE.y[8,-1], solFBE.y[9,-1]]
        
        # solve ODE        
        solFBE_aBE = solve_ivp(lambda t, z: FBEqs_aBE(t, z, nphi, paramsDM,
                                                      G_X, aflog10, bh.TBH(solFBE.y[0,-1],solFBE.y[1,-1]), p_DM, p_X, BR, FO),
                               [aflog10, afmax], v0aBE, method='Radau', rtol=1.e-7, atol=1.e-10)

        npaf = solFBE_aBE.t.shape[0]


        #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
        #       Joining the solutions before and after evaporation       #
        #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

        t    = np.concatenate((solFBE.t[:], solFBE_aBE.t[:]), axis=None)
 
        MBH  = np.concatenate((solFBE.y[0,:], np.full(npaf, solFBE.y[0,0])), axis=None)
        ast  = np.concatenate((solFBE.y[1,:], np.zeros(npaf)), axis=None)
        Rad  = np.concatenate((solFBE.y[2,:], solFBE_aBE.y[0,:]), axis=None)    
        PBH  = np.concatenate((solFBE.y[3,:], np.zeros(npaf)),  axis=None)
        T    = np.concatenate((solFBE.y[4,:], solFBE_aBE.y[1,:]), axis=None)
        NDMT = np.concatenate((solFBE.y[5,:], solFBE_aBE.y[2,:]), axis=None)
        NDMB = np.concatenate((solFBE.y[6,:], solFBE_aBE.y[3,:]), axis=None)
        NX   = np.concatenate((solFBE.y[7,:], solFBE_aBE.y[4,:]), axis=None)
        NDMC = np.concatenate((solFBE.y[8,:], solFBE_aBE.y[5,:]), axis=None)
        NDMH = np.concatenate((solFBE.y[9,:], solFBE_aBE.y[6,:]), axis=None)

        from ulysses.ulsbase import my_kn1, my_kn2

        NDMeq = (mDM**2 * T * kn(2, mDM/T))/(pi**2)

        H   = np.sqrt(8 * pi * bh.GCF * (PBH * 10.**(-3*t) + Rad * 10.**(-4*t))/3.) # Hubble parameter
        npt = T.shape[0]
        TBH = np.concatenate((bh.TBH(solFBE.y[0,:],asi), np.sqrt(mDM**2 + p_DM**2*10.**(2.*(aflog10 - solFBE_aBE.t[:])))), axis=None)

        
        TDM = np.zeros((npt))
        GXt = np.zeros((npt))
                
        Tev=solFBE.y[4,-1]
                
        #------------------------------------------------------------#
        #                                                            #
        #                     Conversion to Oh^2                     #
        #                                                            #
        #------------------------------------------------------------#

        rc = 1.053672e-5*bh.cm_in_invkeV**-3*1.e-18   # Critical density in GeV^3

        T0 = 2.34865e-13  # Temperature today in GeV

        Conf = (bh.gstarS(T0)/bh.gstarS(T[-1]))*(T0/T[-1])**3*(1/rc)

        Oh2   =  (NDMT + NDMB) * nphi * 10.**(-3.*t) * mDM * Conf
        Oh2Th =  NDMT * nphi * 10.**(-3.*t) * mDM * Conf
        Oh2BH =  NDMB * nphi * 10.**(-3.*t) * mDM * Conf

        Oh2C =  NDMC * nphi * 10.**(-3.*t) * mDM * Conf # DM density w/o PBH contribution
        
        a0= 10**t[-1] * (T[-1]/T0)*(bh.gstarS(T[-1])/bh.gstarS(T0))**(1./3.)
        
        
        Tab_T=np.logspace(np.log10(T0),np.log10(Ti),num=1000)
        
        Tab_gT3=bh.gstarS(Tab_T)*Tab_T**3
        
        f_T=interp1d(Tab_gT3,Tab_T)
        
        def find_T(T_1,a_1,a_2):
            return f_T(bh.gstarS(T_1)*T_1**3*(a_1/a_2)**3)
        
        fig, ax = plt.subplots(1, 1, figsize=(8.,6))

        ax.plot(1/T, 10**-t*Rad, label='R')
        ax.plot(1/T, PBH, label='PBH')
        ax.plot(1/T, NDMT * nphi * (mDM), label='DM-TH')
        ax.plot(1/T, NDMB * nphi * np.sqrt(mDM**2), label='DM-BH')
        ax.plot(1/T, NDMT * nphi * (mDM ) + NDMB * nphi * np.sqrt(mDM**2  + p_DM**2 * 10.**(-2.*t)), label='DM-TOT')
        ax.plot(1/T, NX * nphi * np.sqrt(mX**2), label='X')
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.axvline(x=1/Tev, alpha=0.5, color = '#4E2A84', linestyle='--')
                
        ax.set_ylim(max(10**-t*Rad)*1.e-31, max(10**-t*Rad)*1.e-6)
        ax.set_xlim(1/max([100*Tev,100*mDM,100*mX]), 1/T[-1])
        ax.set_xlabel(r"$1/T$")
        ax.set_ylabel(r"$\rho_{i} a^3$")
        ax.legend(loc="lower left", fontsize = "small")
       
        print('Oh2 =', Oh2[-1])
        
        fig.tight_layout()

        plt.show()
               
        # ####### Is DM hot? 
        
        test_hot=Is_DM_hot(Mi, asi, mDM, tau, Sol_t, Tev, max([NDMT[-1],0]), max([NDMB[-1],0]))
        
        print('Is DM hot? -->', test_hot)
                        
        return Oh2[-1]

#-----------------------------------------------------------------------------------
#               Analytics
#-----------------------------------------------------------------------------------
    
    def Omegah2_analytics_FI(self):
        
        Mi     = 10**(self.MPBHi) # PBH initial Mass in grams
        asi    = self.aPBHi       # PBH initial rotation a_star factor
        bi     = 10**(self.bPBHi) # Initial PBH fraction
        Ti     = ((45./(16.*106.75*(pi*bh.GCF)**3.))**0.25) * sqrt(bh.gamma * bh.GeV_in_g/Mi) # Initial Universe temperature
        rRadi  = (pi**2./30.) * bh.gstar(Ti) * Ti**4  # Initial radiation energy density -- assuming a radiation dominated Universe
        rPBHi  = abs(bi/(sqrt(bh.gamma) -  bi))*rRadi # Initial PBH energy density
        

        TBHi   = bh.TBH(Mi, asi)  # Initial BH temperature
        
        mDM     = 10**self.mDM    # DM mass in GeV
        mX      = 10**self.mX     # Mediator Mass
        mf      = 10**self.mf     # SM fermion mass
        sv      = 10**self.sv     # <sigma*v>
        BR      = self.BR         # Branching ratio X -> DM
        g_DM    = self.g_DM       # DM d.o.f.
        model   = self.model      # DM Model
        
        
        # Derive fraction evaporated
        
        FSM_test = bh.fSM(Mi, asi)      # SM contribution
        FDM_test = bh.fDM(Mi, asi, mDM) # DM contribution
        FX_test  = bh.fX(Mi, asi, mX)   # Mediator contribution
        FT_test  = FSM_test + FDM_test + FX_test      # Total Energy contribution
        
        frac_SM=FSM_test/(FT_test)
        
        epsilon_SM=FSM_test
        epsilon=FT_test
        
        if mX < 2.*mDM:
            print("Mediator mass should be bigger than 2*DM mass")
            exit()
        
        gV, gD  = find_gV_gD(sv, BR, mDM, mX, mf)  # gV, gD couplings obtained from sv and BR
        
        #print([gV,gD])
        
        paramsDM=[mDM, mf, mX, gV, gD, g_DM, model]
        
        G_X = Gamma_med(paramsDM)
        
        G_X_SM = Gamma_med(paramsDM) * (1-BR)

                
        #------------------------------------------------------------#
        #                 useful parameters                          #
        #------------------------------------------------------------#

        rc = 1.053672e-5*bh.cm_in_invkeV**-3*1.e-18   # Critical density in GeV^3

        T0 = 2.34865e-13  # Temperature today in GeV
        
        mPL_red=bh.mPL/np.sqrt(8.*np.pi) # reduced Planck mass
        
        mBH_GeV = Mi/bh.GeV_in_g # PBH mass in GeV
        
        T_PBH=bh.TBH(Mi,0) # TBH_i
        
        
        # rename energy densities
        rho_rad_i=rRadi
        rho_PBH_i=rPBHi
        
        #energy fraction at t=ti
        beta=rho_PBH_i/rho_rad_i 
        
        #------------------------------------------------------------#
        #                 useful functions                           #
        #------------------------------------------------------------#

        Tab_T=np.logspace(np.log10(T0),np.log10(Ti),num=1000)
        
        Tab_gT3=bh.gstarS(Tab_T)*Tab_T**3
        
        Tab_gT4=bh.gstar(Tab_T)*Tab_T**4
        
        f_T=interp1d(Tab_gT3,Tab_T)
        find_gT4=interp1d(Tab_gT4,Tab_T)
        
        def find_T(T_1,a_1,a_2):
            return f_T(bh.gstarS(T_1)*T_1**3*(a_1/a_2)**3)
        
        
        # parameter nu^2 in order to get the evaporation temperature
        
        nu2=0.15**2
        
        nu2_bis=0.3**2
        
        #-----
        
        # PBH evaporation inverse lifetime
        Gamma_PBH= bh.kappa*epsilon / (3.*Mi**3)
        
        Tev_theo_2 = find_gT4(epsilon_SM/epsilon*90./(nu2*np.pi**2)*Gamma_PBH**2*mPL_red**2)
        
        Tev_theo_3 = find_gT4(epsilon_SM/epsilon*90./(nu2_bis*np.pi**2)*Gamma_PBH**2*mPL_red**2)
        
        
        # ratio (a_ev / a_in)^(1/3)
        aev_ai_3= nu2*rho_PBH_i/(3.*mPL_red**2*Gamma_PBH**2)#
        
        # a_ev
        aev=(aev_ai_3)**(1./3.)
        
        # a_c, T_c
        ac=aev/(beta*aev*frac_SM)**(2./5.)
        Tc=find_T(Ti,1,ac)
        
        # a_eq, T_eq
        aeq=1/beta
        Teq = find_T(Ti,1,aeq)
        
        # energy density at T_eq and T_c
        rho_PBH_ev = rho_PBH_i*aev_ai_3**(-1)
        rho_PBH_c  = rho_PBH_i*(1/ac)**3.

        # scale factor today
        a0_theo= aev * (Tev_theo_2/T0)*(bh.gstarS(Tev_theo_2)/bh.gstarS(T0))**(1./3.)
        
        
        #### relic ####
        
        # amplitude of the FI annihilation process in the NWA
        amp=2*2*(mX**2+2.*mDM**2)*mX*np.sqrt(mX**2-4.*mDM**2)/G_X
        
        # create variables for the relic density in the different regimes
        r_I = 0
        r_II = 0
        r_III = 0
        r_IV = 0
        
        # assume we start the calculation at a_in=1
        a_i=1
        
        #---------------------
        # REGIME I 
        #---------------------
        factor = (bh.gstarS(Ti)/bh.gstarS(0.786*mX))*a_i**3*(Ti**3 * mPL_red / mX**4)*27.*np.sqrt(10.) / np.sqrt(bh.gstar(0.786*mX)) * np.pi**2*gV**2*gD**2/(1024.*np.pi**6)

        rho_relic=mDM*factor*amp

        r_I = rho_relic / rc * a0_theo**(-3)

        #---------------------
        # REGIME II 
        #---------------------
        
        def Meij2(x):
            
            return 4. * 2.90467 * (mX/(Tc)*(bh.gstarS(Tc)/bh.gstarS(0.786*mX))**(-1./3.))**(-7./2.)

        factor=np.sqrt(3.*mPL_red**2/rho_PBH_c)*(3*gV**2*gD**2)/(2048.*np.pi**4)*Tc*ac**3 * (bh.gstarS(Tc)/bh.gstarS(0.786*mX))**(1./3.)

        rho_relic = mDM*factor*amp*Meij2(1.)

        r_II = rho_relic / rc * a0_theo**(-3)
        
        #---------------------
        # REGIME III 
        #---------------------
        
        def Meij3(x):

            return 4.*1.4746*1e6 * (mX/(Tev_theo_3))**(-11)


        factor=np.sqrt(3.*mPL_red**2/rho_PBH_ev)*(gV**2*gD**2)/(256.*np.pi**4)*Tev_theo_3*aev**3

        rho_relic=mDM*factor*amp*Meij3(1.)
        
        r_III = rho_relic / rc * a0_theo**(-3)
            
        #---------------------
        # REGIME IV 
        #---------------------
        
        factor = (bh.gstarS(T0)/bh.gstarS(0.786*mX))*a0_theo**3*(T0**3 * mPL_red / mX**4)*27.*np.sqrt(10.) / np.sqrt(bh.gstar(0.786*mX)) * np.pi**2*gV**2*gD**2/(1024.*np.pi**6)

        rho_relic = mDM*factor*amp

        r_IV = rho_relic / rc * a0_theo**(-3)
        
        # test thermalization
        thermalization_rate = Gamma_NWA(paramsDM)
        
        
        
        if(Tev_theo_2<=Tc and Tc<=Teq):
            
            # Regime I : radiation domination before BH domination
            if(r_III<= r_I and r_I<= r_II and r_II<=r_IV):
                
                H_max = np.sqrt( np.pi**2/30.*bh.gstar(0.786*mX)*(0.786*mX)**4 / (3.*mPL_red**2) )
                
                # thermal width of the mediator
                therm_X = G_X_SM/H_max
                
                return [r_I, max(thermalization_rate/H_max, therm_X) ]
            
            # Regime II : BH domination before significant injection
            elif(r_III<= r_II and r_II<= r_I and r_I<= r_IV):
                
                a_res = (bh.gstarS(Ti)/bh.gstarS(0.786*mX))**(1./3.)*(Ti/(0.786*mX))
                
                H_max = np.sqrt( rho_PBH_i * a_res **(-3.) / (3.*mPL_red**2) )
                
                # thermal width of the mediator
                therm_X = G_X_SM/H_max
                
                return [r_II, max(thermalization_rate/H_max, therm_X) ]
            
            # Regime III : BH domination and significant injection
            elif(r_II<= r_III and r_III<=r_IV):
                
                a_res = ac * ((0.786*mX)/Tc)**(-8./3.)
                
                H_max = np.sqrt( rho_PBH_i * a_res **(-3.) / (3.*mPL_red**2) )
                
                # thermal width of the mediator
                therm_X = G_X_SM/H_max
                
                return [r_III, max(thermalization_rate/H_max, therm_X)]
            
            # Regime IV : radiation domination after evaporation
            elif(r_IV <= r_III and r_II<=r_I):
                #print('Regime IV')
                H_max = np.sqrt( np.pi**2/30.*bh.gstar(0.786*mX)*(0.786*mX)**4 / (3.*mPL_red**2) )
                
                return [r_IV, max(thermalization_rate/H_max, G_X_SM/H_max) ]
            else:
                
                if(0.786*mX<Tev_theo_2):
                    
                    H_max = np.sqrt( np.pi**2/30.*bh.gstar(0.786*mX)*(0.786*mX)**4 / (3.*mPL_red**2) )
                    
                    # thermal width of the mediator
                    therm_X = G_X_SM/H_max
                
                    
                    return [r_IV, max(thermalization_rate/H_max, therm_X) ]
                
                else:
                    H_max = np.sqrt( np.pi**2/30.*bh.gstar(0.786*mX)*(0.786*mX)**4 / (3.*mPL_red**2) )
                    
                    # thermal width of the mediator
                    therm_X = G_X_SM/H_max
                
                    
                    return [r_I, max(thermalization_rate/H_max, therm_X) ]
                
        else:
            #there was no PBH domination, just calculate the relic density in the standard FI case
            a_i=1
            a0_no_injection= a_i * (Ti/T0)*(bh.gstarS(Ti)/bh.gstarS(T0))**(1./3.)
            
            factor_no_injection = (bh.gstarS(T0)/bh.gstarS(0.786*mX))*a0_no_injection**3*(T0**3 * mPL_red / mX**4)*27.*np.sqrt(10.) / np.sqrt(bh.gstar(0.786*mX)) * np.pi**2*gV**2*gD**2/(1024.*np.pi**6)

            rho_no_injection=mDM*factor_no_injection*amp
            
            Conf_no_injection = (1/rc)/a0_no_injection**3

            H_max = np.sqrt( np.pi**2/30.*bh.gstar(0.786*mX)*(0.786*mX)**4 / (3.*mPL_red**2) )
            
            # thermal width of the mediator
            therm_X = G_X_SM/H_max
                
                    
            return [ rho_no_injection * Conf_no_injection , max(thermalization_rate/H_max, therm_X) ]
   
        
