##################################################################################
#                                                                                #
#      Primordial Black Hole induced leptogenesis. Equations from 2010.XXXXX     #
#                                                                                #
##################################################################################

import ulysses
import numpy as np
from odeintw import odeintw
import pandas as pd
from scipy import interpolate
import matplotlib.pyplot as plt
import scipy.integrate as integrate
from scipy.integrate import quad, ode, solve_ivp, odeint
from scipy.optimize import root
from scipy.special import zeta
from numba import jit

import mpmath

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
mg   = 0.6  # Ficticious gluon mass ---> indicates the QCD phase transition

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

gf = 3.*gnu + 3.*gl + 6.*gq
gs = gH
gv = gW + gZ + gp + gg + ggl

# Constants

c     = 299792.458       # in km/s
gamma = np.sqrt(3.)**-3. # Collapse factor
GCF   = 6.70883e-39      # Gravitational constant in GeV^-2
mPL   = GCF**-0.5        # Planck mass in GeV
v     = 174              # Higgs vev
csp   = 0.35443          # sphaleron conversion factor

# Conversion factors

GeV_in_g     = 1.782661907e-24
Mpc_in_cm    = 3.085677581e24
cm_in_invkeV = 5.067730938543699e7
year_in_s    = 3.168808781402895e-8
GeV_in_invs  = cm_in_invkeV * c * 1.e11

TBH_in_g = (GeV_in_g/(8.*np.pi*GCF))*1.e-13  # Black hole temperature in GeV normalized to a mass of 1.e-13 g
kappa  = 5.34e25 * (1./GeV_in_invs)          # Evaporation constant in g * GeV -- from PRD41(1990)3052

# Contribution to evaporation function per degree of freedom and spin type -- from PRD44(1991)376

f0    = 0.267
f1    = 0.06
f12l  = 0.142
f12n  = 0.147
f2    = 0.007

#----------------------------------------------------#
#      Evaporation function for Schwarzschild BHs    #
#   Equations from PRD 44 (1991) 376 and 1910.07864  #
#----------------------------------------------------#

def fnuW(M):
    
    # Emitted power
    
    b0  = 2.66
    b1  = 6.04
    b12 = 4.53
    b2  = 9.56

    # Equivalent BH masses for each SM particle

    MW = TBH_in_g/mW
    MZ = TBH_in_g/mZ
    MH = TBH_in_g/mH

    Me   = TBH_in_g/me
    Mmu  = TBH_in_g/mmu
    Mtau = TBH_in_g/mtau

    Mu = TBH_in_g/mu
    Md = TBH_in_g/md
    Ms = TBH_in_g/ms
    Mc = TBH_in_g/mc
    Mb = TBH_in_g/mb
    Mt = TBH_in_g/mt

    Mg = TBH_in_g/0.6

    # Contribution from each particle

    fgr =  gg * f2                              # Graviton
    fp  =  gp * f1                              # Photon
    fgl = ggl * f1 * np.exp(-M*1.e-13/(b1*Mg))  # Gluon
    fW  =  gW * f1 * np.exp(-M*1.e-13/(b1*MW))  # W
    fZ  =  gZ * f1 * np.exp(-M*1.e-13/(b1*MZ))  # Z
    fH  =  gH * f0 * np.exp(-M*1.e-13/(b1*MH))  # Higgs

    fnu = 3. * gnu * f12n                       # Active neutrinos
    
    fl  = gl * f12l * (np.exp(-M*1.e-13/(b12*Me))  +
                       np.exp(-M*1.e-13/(b12*Mmu)) +
                       np.exp(-M*1.e-13/(b12*Mtau)))  # Charged leptons

    fq  = gq * f12l * (np.exp(-M*1.e-13/(b12*Mu)) +
                       np.exp(-M*1.e-13/(b12*Md)) +
                       np.exp(-M*1.e-13/(b12*Ms)) +
                       np.exp(-M*1.e-13/(b12*Mc)) +
                       np.exp(-M*1.e-13/(b12*Mb)) +
                       np.exp(-M*1.e-13/(b12*Mt)))    # Quarks
    
    return fgr + fp + fgl + fW + fZ + fH + fnu + fl + fq

# Right handed neutrino contribution

def fnuRH(M, mrh):

    b12 = 4.53

    Mrh = TBH_in_g/mrh

    return 2. * f12n * np.exp(-M*1.e-13/(b12*Mrh))

#----------------------------------------------------#
#        Momentum integrated Hawking rate            #
#----------------------------------------------------#

def GammaRH(M, MRH):

    GM = GCF*M/GeV_in_g # GM in GeV**-1

    al = GM * MRH # Dimensionless gravitational coupling

    zBH = 8. * np.pi * al

    In1 = float(zBH * mpmath.polylog(2, -np.exp(-zBH)))
    In2 = float(mpmath.polylog(3, -np.exp(-zBH)))
    
    In = - In1 - In2

    return  (gnu/(gf + gs + gv + gnu)) * (27/(256. * np.pi**3 * GM)) * In #

#-----------------------------#
#        PBHs lifetime        #
#-----------------------------#

def Itau(M, M1, M2, M3): return M**2./(kappa*(fnuW(M) + fnuRH(M,M1) + fnuRH(M,M2) + fnuRH(M,M3))) 

# Determining the scale fator where PBHs evaporate

def afin(aexp, rPBHi, rRadi, tau):

    a = [10**aexp[0]]
    
    A = -rPBHi * np.sqrt(GCF * (rPBHi + rRadi))
    B = a[0] * rPBHi * np.sqrt(GCF * (a[0]*rPBHi + rRadi))
    C = 2. * rRadi * (np.sqrt(GCF*(rPBHi + rRadi)) - np.sqrt(GCF*(a[0]*rPBHi + rRadi)))
    D = GCF * np.sqrt(6.*np.pi) * rPBHi**2
    
    return [A + B + C - D*tau]

#-------------------------------------#
#    g*(T) and g*S(T) interpolation   #
#-------------------------------------#

gTab = pd.read_table("./Data/gstar.dat",  names=['T','gstar'])

Ttab = gTab.iloc[:,0]
gtab = gTab.iloc[:,1]
tck  = interpolate.splrep(Ttab, gtab, s=0)

def gstar(T): return interpolate.splev(T, tck, der=0)

gSTab = pd.read_table("./Data/gstarS.dat",  names=['T','gstarS'])

TStab = gSTab.iloc[:,0]
gstab = gSTab.iloc[:,1]
tckS  = interpolate.splrep(TStab, gstab, s=0)

def gstarS(T): return interpolate.splev(T, tckS, der = 0)

def dgstarSdT(T): return interpolate.splev(T, tckS, der = 1)

#+++++++++++++++++++++++++++++++++++++++++++++++++#
#             FLRW-Boltzmann Equations            #
#+++++++++++++++++++++++++++++++++++++++++++++++++#

#----------------------------------#
#   Equations before evaporation   #
#----------------------------------#

#@jit
def FBEqs(v, a, rRADi, rPBHi, nphi, M1, M2, M3, epstt, epsmm, epsee, epstm, epste, epsme,c1t,c1m,c1e, widtht, widthm, d, w1, N1Req, dPBH):

    M     = v[0] # PBH mass
    rRAD  = v[1] # Radiation energy density
    rPBH  = v[2] # PBH       energy density
    N1R   = v[3] # nu_R number density
    Tp    = v[4] # Temperature
    Ntt   = v[5]
    Nmm   = v[6]
    Nee   = v[7]
    Ntm   = v[8]
    Nte   = v[9]
    Nme   = v[10]
    c1tc  = np.conjugate(c1t)
    c1mc  = np.conjugate(c1m)
    c1ec  = np.conjugate(c1e)

    #----------------#
    #   Parameters   #
    #----------------#
    
    eSM  = fnuW(M)                  #SM Evaporation contribution
    eRH1 = fnuRH(M, M1)             #1 RH neutrino  Evaporation contribution
    eRH2 = fnuRH(M, M2)             #2 RH neutrino  Evaporation contribution
    eRH3 = fnuRH(M, M3)             #3 RH neutrino  Evaporation contribution
    eD   = eSM + eRH1 + eRH2 + eRH3 #Total Evaporation contribution
    
    H = np.sqrt(8 * np.pi * GCF * (rPBH * 10.**(-3*a) + rRAD * 10.**(-4*a))/3.) # Hubble parameter

    Del = 1. + Tp * dgstarSdT(Tp)/(3. * gstarS(Tp)) # Temperature parameter

    #----------------------------------------------#
    #    Radiation + PBH + Temperature equations   #
    #----------------------------------------------#

    dMda    = - kappa * eD * M**-2/H
    drRADda = - (eSM/eD) * (dMda/M) * 10**a * rPBH
    drPBHda = + (dMda/M) * rPBH
    dTda    = - (Tp/Del) * (1.0 + (eSM/eD) * (dMda/M) * (gstarS(Tp)/gstar(Tp)) * (10**a * rPBH/(4.*rRAD)))

    #----------------------------------------#
    #    RH neutrinos + Lepton asymmetries   #
    #----------------------------------------#

    NTH = (N1R - N1Req) * d/H                              # Thermal contribution
    NBH = (GammaRH(M, M1)/H) * ((GeV_in_g * rPBH/M)/nphi)  # PBH non-thermal contribution      

    dN1Rda = - NTH + NBH
    
    dNttda = epstt * (NTH + NBH*dPBH) - 0.5*w1*(2*c1t*c1tc*Ntt + c1m*c1tc*Ntm + c1e*c1tc*Nte +
                                                np.conjugate(c1m*c1tc*Ntm+c1e*c1tc*Nte))/H
    dNmmda = epsmm * (NTH + NBH*dPBH) - 0.5*w1*(2*c1m*c1mc*Nmm + c1m*c1tc*Ntm + c1e*c1mc*Nme +
                                                np.conjugate(c1m*c1tc*Ntm+c1e*c1mc*Nme))/H
    dNeeda = epsee * (NTH + NBH*dPBH) - 0.5*w1*(2*c1e*c1ec*Nee + c1e*c1mc*Nme + c1e*c1tc*Nte +
                                                np.conjugate(c1e*c1mc*Nme+c1e*c1tc*Nte))/H
    
    dNtmda = epstm * (NTH + NBH*dPBH) - 0.5*w1*(c1t*c1mc*Nmm + c1e*c1mc*Nte + c1m*c1mc*Ntm +
                                             c1mc*c1t*Ntt + c1t*c1tc*Ntm + c1t*c1ec*np.conjugate(Nme))/H - widtht*Ntm/H - widthm*Ntm/H
    
    dNteda = epste * (NTH + NBH*dPBH) - 0.5*w1*(c1t*c1ec*Nee + c1e*c1ec*Nte + c1m*c1ec*Ntm +
                                             c1t*c1ec*Ntt + c1t*c1mc*Nme + c1t*c1tc*Nte )/H - widtht*Nte/H
    
    dNmeda = epsme * (NTH + NBH*dPBH) - 0.5*w1*(c1m*c1ec*Nee + c1e*c1ec*Nme + c1m*c1ec*Nmm +
                                             c1t*c1ec*np.conjugate(Ntm)  + c1m*c1mc*Nme + c1m*c1tc*Nte)/H - widthm*Nme/H

    
    dEqsda = [dMda, drRADda, drPBHda, dN1Rda, dTda, dNttda,  dNmmda, dNeeda, dNtmda, dNteda, dNmeda]

    return [x * np.log(10.) for x in dEqsda]

#----------------------------------#
#    Equations after evaporation   #
#----------------------------------#

#@jit
def FBEqs_aBE(v, a, rRADi, nphi, epstt, epsmm, epsee, epstm, epste, epsme,c1t,c1m,c1e, widtht, widthm, d, w1, N1Req):

    rRAD  = v[0] # Radiation energy density
    N1R   = v[1] # nu_R number density
    Tp    = v[2] # Temperature
    Ntt   = v[3]
    Nmm   = v[4]
    Nee   = v[5]
    Ntm   = v[6]
    Nte   = v[7]
    Nme   = v[8]
    c1tc  = np.conjugate(c1t)
    c1mc  = np.conjugate(c1m)
    c1ec  = np.conjugate(c1e)

    #----------------#
    #   Parameters   #
    #----------------#

    H   = np.sqrt(4. * np.pi**3 * GCF/45.) * np.sqrt(gstar(Tp)) * Tp**2  # Hubble parameter
    Del = 1. + Tp * dgstarSdT(Tp)/(3. * gstarS(Tp))                      # Temperature parameter

    #----------------------------------------#
    #    Radiation + Temperature equations   #
    #----------------------------------------#
    
    drRADda =  0.
    dTda    = - Tp/Del

    #----------------------------------------#
    #    RH neutrinos + Lepton asymmetries   #
    #----------------------------------------#
    
    NTH = (N1R - N1Req) * d/H

    dN1Rda  = - NTH
    
    dNttda   =  epstt * NTH - 0.5*w1*(2*c1t*c1tc*Ntt + c1m*c1tc*Ntm + c1e*c1tc*Nte + np.conjugate(c1m*c1tc*Ntm+c1e*c1tc*Nte))/H
    dNmmda   =  epsmm * NTH - 0.5*w1*(2*c1m*c1mc*Nmm + c1m*c1tc*Ntm + c1e*c1mc*Nme + np.conjugate(c1m*c1tc*Ntm+c1e*c1mc*Nme))/H
    dNeeda   =  epsee * NTH - 0.5*w1*(2*c1e*c1ec*Nee + c1e*c1mc*Nme + c1e*c1tc*Nte + np.conjugate(c1e*c1mc*Nme+c1e*c1tc*Nte))/H

    dNtmda   = epstm * NTH - 0.5*w1*(c1t*c1mc*Nmm + c1e*c1mc*Nte + c1m*c1mc*Ntm +
                                     c1mc*c1t*Ntt + c1t*c1tc*Ntm + c1t*c1ec*np.conjugate(Nme) )/H - widtht*Ntm/H - widthm*Ntm/H
    
    dNteda   = epste * NTH - 0.5*w1*(c1t*c1ec*Nee + c1e*c1ec*Nte + c1m*c1ec*Ntm +
                                      c1t*c1ec*Ntt + c1t*c1mc*Nme + c1t*c1tc*Nte)/H - widtht*Nte/H
    
    dNmeda   = epsme * NTH - 0.5*w1*(c1m*c1ec*Nee + c1e*c1ec*Nme + c1m*c1ec*Nmm +
                                      c1t*c1ec*np.conjugate(Ntm)  + c1m*c1mc*Nme + c1m*c1tc*Nte)/H - widthm*Nme/H


    dEqsda = [drRADda, dN1Rda, dTda, dNttda,  dNmmda, dNeeda, dNtmda, dNteda, dNmeda]

    dEqsda = [x * np.log(10.) for x in dEqsda]
    
    return dEqsda

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#                           ULYSSES class for PBH-Leptogenesis                            #     
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

class EtaB_PBH(ulysses.ULSBase):
    """
    Primordial black hole driven leptogenesis.
    density matrix equation (1DME) non-zero thermal width (flavours charged leptons accounted for) with one decaying sterile.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._meff1 = None
        self._k1 = None

        self.MPBHi = None # Log10[M/1g]
        self.bPBHi = None # Log10[beta']

        self.pnames = ['m',  'M1', 'M2', 'M3', 'delta', 'a21', 'a31', 'x1', 'x2', 'x3', 'y1', 'y2', 'y3',
                       't12', 't13', 't23', 'MPBHi', 'bPBHi']


    def setParams(self, pdict):
        """
        This set the model parameters. pdict is expected to be a dictionary
        """
        super().setParams(pdict)
        self.MPBHi = pdict["MPBHi"]
        self.bPBHi = pdict["bPBHi"]
        
    #*************************************#
    #               Before                #
    #*************************************#
    def RHS(self, y0, a, epstt, epsmm, epsee, epstm, epste, epsme, c1t, c1m, c1e, rRadi, rPBHi, nphi):
        Tp  = y0[4]                              # Ambient Temperature
        TBH = 1./(8.*np.pi*GCF*(y0[0]/GeV_in_g)) # BH temperature
        z   = self.M1/Tp
        zBH = self.M1/TBH

        from ulysses.ulsbase import my_kn2, my_kn1
    
        self._d    = np.real(self.Gamma1 * my_kn1(z)/my_kn2(z))
        self._dPBH = np.real(self.Gamma1 * my_kn1(zBH)/my_kn2(zBH))
        self._w1   = self._d * (0.25 * my_kn2(z) * z**2)
        self._n1eq = (10**(3*a) * self.M1**2 * Tp * my_kn2(z))/(np.pi**2)/nphi

        # thermal widths are non-zero needed for flavour effects
        widtht = 485e-10 * self.MP/self.M1
        widthm = 1.7e-10 * self.MP/self.M1
        
        #Yukawa entries
        c1tc    = np.conjugate(c1t)
        c1mc    = np.conjugate(c1m)
        c1ec    = np.conjugate(c1e)
        
        return FBEqs(y0, a, rRadi, rPBHi, nphi, self.M1, self.M2, self.M3, epstt, epsmm, epsee, epstm, epste, epsme,
                     c1t, c1m, c1e, widtht, widthm, self._d, self._w1, self._n1eq, self._dPBH)

    #*************************************#
    #                After                #
    #*************************************#
    def RHS_aBE(self, y0, a, epstt, epsmm, epsee, epstm, epste, epsme, c1t, c1m, c1e, rRadi, nphi):
        Tp   = y0[2] # Temperature
        k    = np.real(self.k1)
        z    = np.real(self.M1/Tp)

        from ulysses.ulsbase import my_kn2, my_kn1
   
        self._d    =  np.real(self.Gamma1 * my_kn1(z)/my_kn2(z))
        self._w1   =  self._d * (0.25 * my_kn2(z) * z**2)
        self._n1eq =  (10**(3*a) * self.M1**2 * Tp * my_kn2(z))/(np.pi**2)/nphi #self.N1Eq(z)

        # thermal widths are  non-zero needed for flavour effects
        widtht = 485e-10*self.MP/self.M1
        widthm = 1.7e-10*self.MP/self.M1
        
        return FBEqs_aBE(y0, a, rRadi, nphi, epstt, epsmm, epsee, epstm, epste, epsme, c1t, c1m, c1e, widtht, widthm,
                         self._d, self._w1, self._n1eq)

    @property
    def EtaB(self):
        #Define fixed quantities for BEs
        epstt = np.real(self.epsilon1ab(2,2))
        epsmm = np.real(self.epsilon1ab(1,1))
        epsee = np.real(self.epsilon1ab(0,0))
        epstm =         self.epsilon1ab(2,1)
        epste =         self.epsilon1ab(2,0)
        epsme =         self.epsilon1ab(1,0)

        # Yukawa couplings
        c1t = self.c1a(2)
        c1m = self.c1a(1)
        c1e = self.c1a(0)
        
        Mi     = 10**(self.MPBHi) # PBH initial Mass in grams
        bi     = 10**(self.bPBHi) # Initial PBH fraction
        Ti     = ((45./(16.*106.75*(np.pi*GCF)**3.))**0.25) * np.sqrt(gamma * GeV_in_g/Mi) # Initial Universe temperature
        rRadi  = (np.pi**2./30.) * gstar(Ti) * Ti**4  # Initial radiation energy density -- assuming a radiation dominated Universe
        rPBHi  = (bi/(np.sqrt(gamma) -  bi))*rRadi    # Initial PBH energy density
        nphi   = (2.*zeta(3)/np.pi**2)*Ti**3          # Initial photon number density
        N1Ri   = 0.                                   # Initial RH neutrino number density
        Ntti   = 0.                                   # Initial B-L number densities
        Nmmi   = 0.
        Neei   = 0.
        Ntmi   = 0.
        Ntei   = 0.
        Nmei   = 0.

        parms1 = np.array([epstt,epsmm,epsee,epstm,epste,epsme,c1t,c1m,c1e, np.real(rRadi),np.real(rPBHi),np.real(nphi)],
                          dtype=np.complex128)

        v0 = [Mi, rRadi, rPBHi, N1Ri, Ti, Ntti, Nmmi, Neei, Ntmi, Ntei, Nmei]
    
        tau = integrate.quad(Itau, GeV_in_g*mPL, Mi, args=(self.M1, self.M2, self.M3)) # PBH lifetime in s
        af = root(afin, [20.], args = (rPBHi, rRadi, tau[0]), method='hybr')           # Scale factor in which BHs evaporate
        aflog10 = 0.995*af.x[0]
        
        #------------------------------------------------------------#
        #                                                            #
        #                     Before BH evaporation                  #
        #                                                            #
        #------------------------------------------------------------#        

        # time points
        t1 = np.linspace(0., aflog10, num=1000, endpoint=True)

        # solve ODE
        solFBE = odeintw(self.RHS, v0, t1, args=tuple(parms1), rtol=1.e-13, atol=1.e-13)

        #------------------------------------------------------------#
        #                                                            #
        #                      After BH evaporation                  #
        #                                                            #
        #------------------------------------------------------------#

        # time points

        npaf = 500

        azmax = aflog10 + 2.*np.log10(25.*solFBE[-1,4]/self.M1) # Determining the z=M1/T at which PBH evaporate 
        afmax = max(aflog10, azmax)
       
        t2 = np.linspace(aflog10, afmax, num=npaf)

        v0aBE = [solFBE[-1,1], solFBE[-1,3], solFBE[-1,4], solFBE[-1,5], solFBE[-1,6],
                 solFBE[-1,7], solFBE[-1,8], solFBE[-1,9], solFBE[-1,10]]

        parms2 = np.array([epstt,epsmm,epsee,epstm,epste,epsme,c1t,c1m,c1e, np.real(solFBE[-1,1]),np.real(nphi)], dtype=np.complex128)

        # solve ODE
        solFBE_aBE = odeintw(self.RHS_aBE, v0aBE, t2, args=tuple(parms2), rtol=1.e-13, atol=1.e-13)

        #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
        #       Joining the solutions before and after evaporation       #
        #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
        
        t = np.concatenate((t1, t2), axis=None)

        MBH = np.concatenate((solFBE[:,0], np.zeros(npaf)), axis=None)
        Rad = np.concatenate((solFBE[:,1], solFBE_aBE[:,0]), axis=None)    
        PBH = np.concatenate((solFBE[:,2], np.zeros(npaf)),  axis=None)
        T   = np.concatenate((solFBE[:,4], solFBE_aBE[:,2]), axis=None)
       
        NRH   = np.concatenate((solFBE[:,3], solFBE_aBE[:,1]), axis=None)
        NBLtt = np.concatenate((solFBE[:,5], solFBE_aBE[:,3]), axis=None)
        NBLmm = np.concatenate((solFBE[:,6], solFBE_aBE[:,4]), axis=None)
        NBLee = np.concatenate((solFBE[:,7], solFBE_aBE[:,5]), axis=None)

        
        #------------------------------------------------------------#
        #                                                            #
        #                     Conversion to eta_B                    #
        #                                                            #
        #------------------------------------------------------------#

        gstarSrec = gstarS(0.3e-9) # d.o.f. at recombination
        gstarSoff = gstarS(T[-1])  # d.o.f. at the end of leptogenesis
     
        SMspl       = 28./79.
        zeta3       = zeta(3)
        ggamma      = 2.
        coeffNgamma = ggamma*zeta3/np.pi**2
        Ngamma      = coeffNgamma*(10**t*T)**3
        coeffsph    = (SMspl * gstarSrec)/(gstarSoff * Ngamma) 

        nb = coeffsph*(NBLtt +  NBLmm +  NBLee)*nphi

        return np.real(nb[-1])
