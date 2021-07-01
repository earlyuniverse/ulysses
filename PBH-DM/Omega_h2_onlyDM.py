###################################################################################################
#                                                                                                 #
#                         Primordial Black Hole + Dark Matter Generation.                         #
#                                    Only DM from evaporation                                     #
#                                                                                                 #
#         Authors: Andrew Cheek, Lucien Heurtier, Yuber F. Perez-Gonzalez, Jessica Turner         #
#                    Based on: arXiv:2107.xxxxx (P1) and  arXiv:2107.xxxxx (P2)                   #
#                                                                                                 #
###################################################################################################

import ulysses
import numpy as np
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
#          - 'Mi'   : Primordial BH initial Mass in grams
#
#          - 'ai'   : Primordial BH initial angular momentum a*
#
#          - 'bi'   : Primordial BH initial fraction beta^prim
#
#          - 'mDM'  : DM Mass in GeV
#
#          - 'g_DM' : DM degrees of freedom
#
#-----------------------------------------------------------------

#--------------------------   Credits  -----------------------------#
#
#      If using this code, please cite:
#
#      - Arxiv:2107.xxxxx    and    Arxiv:2107.xxxxx
#
#      - JCAP 12 (2017) 013 • e-Print: 1706.03118 (WDM constraints)
#
#-------------------------------------------------------------------#

def PlanckMass_A(t, v, Mi):

    eps = 1.e-2

    if (eps*Mi > bh.MPL): Mst = eps*Mi
    else: Mst = bh.MPL
    
    return v[0] - Mst # Function to stop the solver if the BH is equal or smaller than the Planck mass


def PlanckMass_B(t, v, Mi):

    eps = 1.e-2

    if (eps*Mi > bh.MPL): Mst = eps*Mi
    else: Mst = bh.MPL
    
    return Mi*10.**v[0] - Mst # Function to stop the solver if the BH is equal or smaller than the Planck mass


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

    integ_p = integrate.quad(Integ_p, -10., tau, args=(pars), epsabs=1.e-10, epsrel=1.e-5)
    integ_n = integrate.quad(Integ_n, -10., tau, args=(pars), epsabs=1.e-10, epsrel=1.e-5)

    return (bh.kappa * integ_p[0]/bh.GeV_in_g)/(integ_n[0] + 1.e-100)

    


#########################################################

#----------------------------------#
#   Equations before evaporation   #
#----------------------------------#

def FBEqs(a, v, nphi, mDM, Mi, ailog10):

    Mt    = v[0] # PBH mass
    ast   = v[1] # PBH ang mom
    rRAD  = v[2] # Radiation energy density
    rPBH  = v[3] # PBH energy density
    Tp    = v[4] # Temperature
    NDMH  = v[5] # PBH-induced DM number density

    M = Mi * 10.**Mt

    #print(M)

    #----------------#
    #   Parameters   #
    #----------------#
    
    FSM = bh.fSM(M, ast)      # SM contribution
    FDM = bh.fDM(M, ast, mDM) # DM contribution
    FT  = FSM + FDM           # Total Energy contribution

    GSM = bh.gSM(M, ast)      # SM contribution
    GDM = bh.gDM(M, ast, mDM) # DM contribution
    GT  = GSM + GDM           # Total Angular Momentum contribution
    
    H   = np.sqrt(8 * pi * bh.GCF * (rPBH * 10.**(-3*(a + ailog10)) + rRAD * 10.**(-4*(a + ailog10)))/3.) # Hubble parameter
    Del = 1. + Tp * bh.dgstarSdT(Tp)/(3. * bh.gstarS(Tp)) # Temperature parameter

    from ulysses.ulsbase import my_kn2, my_kn1
    
    TH  = bh.TBH(M, ast) # Hawking Temperature
    z   = mDM/Tp       
    zBH = mDM/TH   
    
    #----------------------------------------------#
    #    Radiation + PBH + Temperature equations   #
    #----------------------------------------------#

    dMtda   = - bh.kappa * FT * M**-3/H/log(10.)   
    dastda  = - ast * bh.kappa * M**-3 * (GT - 2.*FT)/H
    drRADda = - (FSM/FT) * (dMtda * log(10.)) * 10.**(a + ailog10) * rPBH
    drPBHda = + (dMtda * log(10.)) * rPBH
    dTda    = - (Tp/Del) * (1.0 + (bh.gstarS(Tp)/bh.gstar(Tp))*((FSM/FT)*(dMtda * log(10.))* 10.**(a + ailog10) * rPBH)/(4.*rRAD))

    #-----------------------------------------#
    #           Dark Matter Equations         #
    #-----------------------------------------#
    
    dNDMHda = (bh.Gamma_F(M, ast, mDM)/H)*(rPBH/(M/bh.GeV_in_g))/nphi # PBH-induced contribution w/o contact
    
    ##########################################################    
    
    dEqsda = [dMtda, dastda, drRADda, drPBHda, dTda, dNDMHda]

    return [x * log(10.) for x in dEqsda]

#----------------------------------#
#    Equations after evaporation   #
#----------------------------------#

def FBEqs_aBE(a, v, nphi, mDM, a_evap, T_bh_in, p_DM):

    rRAD = v[0] # Radiation energy density
    Tp   = v[1] # Temperature
    NDMH = v[2] # Thermal DM number density w/o PBH contribution
    
    #----------------#
    #   Parameters   #
    #----------------#

    H   = sqrt(8 * pi * bh.GCF * (rRAD * 10.**(-4*a))/3.)    # Hubble parameter
    Del = 1. + Tp * bh.dgstarSdT(Tp)/(3. * bh.gstarS(Tp))          # Temperature parameter
            
    z = mDM/Tp
    
    ####
    # Boost factor of DM
    ####
    p_DM_0 = p_DM * 10.**(a_evap - a)
    E_DM = sqrt(mDM**2 + p_DM_0**2)
    
    #----------------------------------------#
    #    Radiation + Temperature equations   #
    #----------------------------------------#

    drRADda = 0.
    dTda    = - Tp/Del
        
    #-----------------------------------------#
    #           Dark Matter Equations         #
    #-----------------------------------------#

    dNDMHda = 0.                              # PBH-induced contribution w/o contact
        
    dEqsda = [drRADda, dTda, dNDMHda]

    return [x * log(10.) for x in dEqsda]

#------------------------------------------------------------------------------------------------------------------#
#                                            Input parameters                                                      #
#------------------------------------------------------------------------------------------------------------------#
class FrInPBH:

    def __init__(self, MPBHi, aPBHi, bPBHi, mDM, g_DM):

        self.MPBHi  = MPBHi # Log10[M/1g]
        self.aPBHi  = aPBHi # a_star
        self.bPBHi  = bPBHi # Log10[beta']
        self.mDM    = mDM
        self.g_DM   = g_DM
    
#-------------------------------------------------------------------------------------------------------------------------------------#
#                                                       Input parameters                                                              #
#-------------------------------------------------------------------------------------------------------------------------------------#
    
    def Omegah2(self):
        
        Mi     = 10**(self.MPBHi) # PBH initial Mass in grams
        asi    = self.aPBHi       # PBH initial rotation a_star factor
        bi     = 10**(self.bPBHi) # Initial PBH fraction
        Ti     = ((45./(16.*106.75*(pi*bh.GCF)**3.))**0.25) * sqrt(bh.gamma * bh.GeV_in_g/Mi) # Initial Universe temperature
        rRadi  = (pi**2./30.) * bh.gstar(Ti) * Ti**4  # Initial radiation energy density -- assuming a radiation dominated Universe
        rPBHi  = abs(bi/(sqrt(bh.gamma) -  bi))*rRadi # Initial PBH energy density
        nphi   = (2.*zeta(3)/pi**2)*Ti**3             # Initial photon number density

        NDMHi  = 0.0

        TBHi   = bh.TBH(Mi, asi)  # Initial BH temperature
        
        mDM     = 10**self.mDM    # DM mass in GeV
        mX      = 3*10**self.mDM  # Mediator Mass
        g_DM    = self.g_DM       # DM d.o.f.

        paramsDM=[mDM, mX, g_DM] 
    
        #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
        #                                           Solving the equations                                                   #
        #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

        ailog10 = 0.

        Min  = Mi
        asin = asi

        tBE    = []
        MBHBE  = []
        astBE  = []
        RadBE  = []
        PBHBE  = []
        TBE    = []
        NDMHBE = []

        i = 0
        
        while Mi >= 10. * bh.MPL:# Loop on the solver such that BH mass reaches M_Planck

            #--------------------------------------------------------------------------------#
            #         Computing PBH lifetime and scale factor in which BHs evaporate         #
            #--------------------------------------------------------------------------------#
            
            MPL_A = lambda t, x:PlanckMass_A(t, x, Mi)
            MPL_A.terminal  = True
            MPL_A.direction = -1.
            
            tau_sol = solve_ivp(fun=lambda t, y: bh.ItauFO(t, y, mDM), t_span = [-80, 40.], y0 = [Mi, asi], 
                                 events=MPL_A, rtol=1.e-5, atol=1.e-20, dense_output=True)
            
            if i == 0:
                Sol_t = tau_sol.sol # Solutions for obtaining <p>
                tau = tau_sol.t[-1] # Log10@PBH lifetime in inverse GeV
            
            if bi > 1.e-19*(1.e9/Mi):
                af = root(bh.afin, [40.], args = (rPBHi, rRadi, 10.**tau, 0.), method='lm', tol=1.e-40) # Scale factor 
                aflog10 = af.x[0]            
            else:
                afw = np.sqrt(1. + 4.*10.**tau*np.sqrt(2.*np.pi*bh.GCF*rRadi/3.))
                aflog10 = np.log10(afw)

            #print(aflog10)
            
            #-----------------------------------------#
            #          Before BH evaporation          #
            #-----------------------------------------#

            MPL_B = lambda t, x:PlanckMass_B(t, x, Mi)
            MPL_B.terminal  = True
            MPL_B.direction = -1.
            
            v0 = [0., asi, rRadi, rPBHi, Ti, NDMHi]
            
            # solving ODE
            solFBE = solve_ivp(lambda t, z: FBEqs(t, z, nphi, mDM, Mi, ailog10),
                               [0., 1.05*abs(aflog10)], v0, method='BDF', events=MPL_B, rtol=1.e-5, atol=1.e-20)

            if solFBE.t[-1] < 0.:
                print(solFBE)
                print(afw, tau, 1.05*aflog10)
                break

            # Concatenating solutions
            
            tBE    = np.append(tBE,    np.log10(10.**solFBE.t[:] + 10.**ailog10))
            MBHBE  = np.append(MBHBE,  10.**(solFBE.y[0,:])*Mi)
            astBE  = np.append(astBE,  solFBE.y[1,:])
            RadBE  = np.append(RadBE,  solFBE.y[2,:])
            PBHBE  = np.append(PBHBE,  solFBE.y[3,:])
            TBE    = np.append(TBE,    solFBE.y[4,:])
            NDMHBE = np.append(NDMHBE, solFBE.y[5,:])

            # Updating values of initial parameters
            
            Mi    = 10.**(solFBE.y[0,-1])*Mi
            asi   = solFBE.y[1,-1]
            rRadi = solFBE.y[2,-1]
            rPBHi = solFBE.y[3,-1]
            Ti    = solFBE.y[4,-1]
            NDMHi = solFBE.y[5,-1]
            
            ailog10 += solFBE.t[-1]

            i += 1

            if i > 100:
                aflog10 = ailog10
                print("I'm stuck!", Mi, bi)
                print()
                break

        else:
            aflog10 = ailog10# We update the value of log(a) at which PBHs evaporate


        #+++++++++++++++++++++++++++#
        #      Average momentum     #
        #+++++++++++++++++++++++++++#

        p_DM = p_average_DM(Min, asin, mDM, tau, Sol_t)
        
        #-----------------------------------------#
        #           After BH evaporation          #
        #-----------------------------------------#
        
        Tfin = 1.e-3 # Final plasma temp in GeV
        
        azmax = aflog10 + np.log10(np.cbrt(bh.gstarS(TBE[-1])/bh.gstarS(Tfin))*(TBE[-1]/Tfin))
        afmax = max(aflog10, azmax)

        v0aBE = [RadBE[-1], TBE[-1], NDMHBE[-1]]
        
        # solve ODE        
        solFBE_aBE = solve_ivp(lambda t, z: FBEqs_aBE(t, z, nphi, mDM, aflog10, bh.TBH(MBHBE[-1], astBE[-1]), p_DM),
                               [aflog10, afmax], v0aBE, method='Radau', rtol=1.e-7, atol=1.e-20)

        npaf = solFBE_aBE.t.shape[0]

        #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
        #       Joining the solutions before and after evaporation       #
        #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

        t    = np.concatenate((tBE, solFBE_aBE.t[:]), axis=None)
 
        MBH  = np.concatenate((MBHBE,  np.full(npaf, solFBE.y[0,0])), axis=None)
        ast  = np.concatenate((astBE,  np.zeros(npaf)), axis=None)
        Rad  = np.concatenate((RadBE,  solFBE_aBE.y[0,:]), axis=None)    
        PBH  = np.concatenate((PBHBE,  np.zeros(npaf)),  axis=None)
        T    = np.concatenate((TBE,    solFBE_aBE.y[1,:]), axis=None)
        NDMH = np.concatenate((NDMHBE, solFBE_aBE.y[2,:]), axis=None)
        
        Tev  = TBE[-1]

        plt.plot(mDM/T, 10**-t*Rad, label='R', lw = 2)
        plt.plot(mDM/T, PBH, label='PBH', lw = 2, color = '#1e1f26')
        plt.plot(mDM/T, NDMH * nphi * np.sqrt(mDM**2 + p_DM**2 * 10.**(-2.*t)), label='DM-BH', lw = 2 , color = '#00b159')
        plt.xscale('log')
        plt.yscale('log')
        plt.title(r"$a_\star = 0.$")
        plt.axvline(x=mDM/Tev, alpha=0.5, color = '#4E2A84', linestyle='--')
        plt.ylim(max(10**-t*Rad)*1.e-30, max(10**-t*Rad)*1.e1)
        plt.xlabel(r"$m_{\rm DM}/T$")
        plt.ylabel(r"$\rho_{i} a^4~[GeV^4]$")
        plt.legend(loc='upper right', ncol=2, fontsize = "small")
        #plt.xlim(1e-13,1e2)

        plt.tight_layout()
        plt.show()
                
        #------------------------------------------------------------#
        #                                                            #
        #                     Conversion to Oh^2                     #
        #                                                            #
        #------------------------------------------------------------#

        rc = 1.053672e-5*bh.cm_in_invkeV**-3*1.e-18   # Critical density in GeV^3

        T0 = 2.34865e-13  # Temperature today in GeV

        Conf = (bh.gstarS(T0)/bh.gstarS(T[-1]))*(T0/T[-1])**3*(1/rc) 

        Oh2  =  NDMH * nphi * 10.**(-3.*t) * mDM * Conf

        return Oh2[-1] # Final relic abundance
