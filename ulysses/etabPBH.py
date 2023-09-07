##################################################################################
#                                                                                #
#                    Primordial Black Hole induced leptogenesis.                 #
#                  High Scale Scenario, including DL = 2 scattering              #
#                                                                                #
##################################################################################

import ulysses
import numpy as np
from odeintw import odeintw
from scipy import interpolate
import scipy.integrate as integrate
from scipy.integrate import quad, ode, solve_ivp, odeint
from scipy.optimize import root
from scipy.special import zeta
from numba import njit

import ulysses.BHProp as bh #Schwarzschild and Kerr BHs library

from termcolor import colored

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#                                                     FLRW-Boltzmann Equations                                                       #
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#----------------------------------#
#   Equations before evaporation   #
#----------------------------------#

def FBEqs(x, v, nphi, M1, M2, M3, eps, d1, w1, N1Req, dPBH1, WashDL2, xilog10):

    M    = v[0] # PBH mass in g
    ast  = v[1] # PBH angular momentum
    rRAD = v[2] # Radiation energy density
    rPBH = v[3] # PBH       energy density
    Tp   = v[4] # Plasma Temperature
    
    N1RT = v[5] # nu_R thermal number density
    N1RB = v[6] # nu_R number density from PBH evaporatin
    
    NBL  = v[7] # single Flavor B-L asymmetry

    #----------------#
    #   Parameters   #
    #----------------#
    
    M_GeV = M/bh.GeV_in_g     # PBH mass in GeV
    nPBH  = rPBH/M_GeV        # PBH number density in GeV^3
    
    xff = x + xilog10         # Logarithm_10 of scale factor, xilog is the log10 of the initial scale factor

    eps1tt, eps1mm, eps1ee, eps1tm, eps1te, eps1me = eps   # CP violation elements

    # Evaporation functions for Black Hole mass and spin
    
    FSM  = bh.fSM(M, ast)           # SM Evaporation contribution
    FRH1 = bh.fRH(M, ast, M1)       # 1 RH neutrino  Evaporation contribution
    FRH2 = bh.fRH(M, ast, M2)       # 2 RH neutrino  Evaporation contribution
    FRH3 = bh.fRH(M, ast, M3)       # 3 RH neutrino  Evaporation contribution
    FT   = FSM + FRH1 + FRH2 + FRH3 # Total Evaporation contribution

    GSM  = bh.gSM(M, ast)           # SM Evaporation contribution
    GRH1 = bh.gRH(M, ast, M1)       # 1 RH neutrino  Evaporation contribution
    GRH2 = bh.gRH(M, ast, M2)       # 2 RH neutrino  Evaporation contribution
    GRH3 = bh.gRH(M, ast, M3)       # 3 RH neutrino  Evaporation contribution
    GT   = GSM + GRH1 + GRH2 + GRH3 # Total Evaporation contribution
    
    H = np.sqrt(8 * np.pi * bh.GCF * (rPBH*10.**(-3*xff) + rRAD*10.**(-4*xff))/3.) # Hubble parameter

    Del = 1. + Tp * bh.dgstarSdT(Tp)/(3. * bh.gstarS(Tp)) # Temperature parameter

    #----------------------------------------------#
    #    Radiation + PBH + Temperature equations   #
    #----------------------------------------------#   

    dM_GeVdx = - FT/(bh.GCF**2 * M_GeV**2)/H   
    dastdx   = - ast * (GT - 2.*FT)/(bh.GCF**2 * M_GeV**3)/H
    
    drRADdx  = - (FSM/FT) * (dM_GeVdx/M_GeV) * 10**xff * rPBH
    drPBHdx  = + (dM_GeVdx/M_GeV) * rPBH
    
    dTdx     = - (Tp/Del) * (1.0 - (bh.gstar(Tp)/bh.gstarS(Tp))*(drRADdx/(4.*rRAD)))

    #----------------------------------------#
    #              RH neutrinos              #
    #----------------------------------------#

    NTH1 = (N1RT - N1Req) * d1/H
    NBH1 = N1RB * dPBH1/H

    dN1RTdx = -NTH1                                          # Thermal contribution
    dN1RBdx = -NBH1 + (bh.Gamma_F(M, ast, M1)/H)*(nPBH/nphi) # PBH-induced contribution, normalized wrt the initial photon density nphi

    #----------------------------------------#
    #            Lepton asymmetries          #
    #----------------------------------------#

    dNBLdx = (eps1tt+eps1mm+eps1ee)*(NTH1 + NBH1) - (w1 + WashDL2)*NBL/H
  
    #Equations

    kappa = bh.GeV_in_g # Conversion factor to have Mass equation rate for PBH mass in g
    
    dEqsdx = [kappa * dM_GeVdx, dastdx, drRADdx, drPBHdx, dTdx, dN1RTdx, dN1RBdx, dNBLdx]

    return [xeq * np.log(10.) for xeq in dEqsdx]

#----------------------------------#
#    Equations after evaporation   #
#----------------------------------#

def FBEqs_aBE(x, v, nphi, M1,M2,M3, eps, d1, w1, N1Req, dPBH1, WashDL2):

    rRAD  = v[0] # Radiation energy density
    Tp    = v[1] # Plasma Temperature
    
    N1RT  = v[2] # nu_R therma number density
    N1RB  = v[3] # nu_R number density from PBH evaporation
    
    NBL   = v[4] # single Flavor B-L asymmetry

    #----------------#
    #   Parameters   #
    #----------------#

    eps1tt, eps1mm, eps1ee, eps1tm, eps1te, eps1me = eps # CP violation elements

    H   = np.sqrt(8 * np.pi * bh.GCF * ((N1RT + N1RB)*M1 * 10.**(-3*x) + rRAD * 10.**(-4*x))/3.)  # Hubble parameter
    Del = 1. + Tp * bh.dgstarSdT(Tp)/(3. * bh.gstarS(Tp))                                         # Temperature parameter

    #----------------------------------------#
    #    Radiation + Temperature equations   #
    #----------------------------------------#
    
    drRADdx = 0.
    dTdx    = - Tp/Del

    #----------------------------------------#
    #              RH neutrinos              #
    #----------------------------------------#

    NTH1 = (N1RT - N1Req) * d1/H  # Thermal contribution
    NBH1 = N1RB * dPBH1/H         # PBH non-thermal contribution

    dN1RTdx = -NTH1
    dN1RBdx = -NBH1
    
    #----------------------------------------#
    #            Lepton asymmetries          #
    #----------------------------------------#

    dNBLdx = ((eps1tt+eps1mm+eps1ee)* (NTH1 + NBH1) - (w1 + WashDL2)*NBL/H)
    #Equations
    
    dEqsdx = [drRADdx, dTdx, dN1RTdx, dN1RBdx, dNBLdx]

    dEqsdx = [x * np.log(10.) for x in dEqsdx]
    
    return dEqsdx

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#                           ULYSSES class for PBH-Leptogenesis                            #     
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

class EtaB_PBH(ulysses.ULSBase):
    """
    Primordial black hole assisted Leptogenesis.
    One-flavoured BE with 1 Right-handed Neutrino. Including the DL=2 washout term.
    See arXiv:2010.03565 and arXiv:2203.08823
    """

    def shortname(self): return "1BE1F_PBH"

    def evolname(self): return "a"

    def flavourindices(self): return [1, 2]

    def flavourlabels(self): return ["$N^{\\rm B-L}_{\\rm TH}$", "$N^{\\rm B-L}_{\\rm PBH}$"]

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        """
        This model requires three additional parameters, PBH mass and spin and PBH initial abundance
        """
        self.MPBHi = None # Log10[M/1g]
        self.aPBHi = None # Initial spin factor a_*
        self.bPBHi = None # Log10[beta']

        self.pnames = ['m',  'M1', 'M2', 'M3', 'delta', 'a21', 'a31', 'x1', 'x2', 'x3', 'y1', 'y2', 'y3',
                       't12', 't13', 't23', 'MPBHi', 'aPBHi', 'bPBHi']

        #------------------------------------------------------------#
        #            Inverse Time dilatation factor <M/E>            #
        #        Averaged wrt the Hawking emission spectrum          #
        #------------------------------------------------------------#

        import os
        data_dir = os.path.dirname(ulysses.__file__)

        MEav_f   = os.path.join(data_dir, "timedil.txt")
        MEavTab  = np.loadtxt(MEav_f)

        self.MEav_ = interpolate.splrep(MEavTab[:,0],  MEavTab[:,1],  s=0)

    def ME(self,zBH):

        LzBH = np.log10(zBH)
        
        if LzBH < -4.:
            return 10.**(-0.6267457 + 0.9999617*LzBH)
        elif LzBH >= -4. and LzBH <= 2.85:
            return interpolate.splev(LzBH, self.MEav_, der=0)
        else:
            return 1.

    #-------------------------------------------------------------------#
    #        Function to stop the solver if the BH is equal than        #
    #               1% the initial mass  or the Planck mass             #
    #-------------------------------------------------------------------#
    
    def StopMass(self, t, v, Mi):
        
        eps = 0.01
        
        if (eps*Mi > bh.MPL): Mst = eps*Mi
        else: Mst = bh.MPL
        
        return v[0] - Mst
    
    def setParams(self, pdict):
        super().setParams(pdict)
        self.MPBHi = pdict["MPBHi"]
        self.aPBHi = pdict["aPBHi"]
        self.bPBHi = pdict["bPBHi"]
        
    #********************************************************#
    #        Equations Before  PBH evaporation               #
    #********************************************************#
    def RHS(self, x, y0, eps, nphi, xilog10): # x is the Log10 of the scale factor
        
        Tp   = y0[4]                 # Plasma Temperature
        TBH  = bh.TBH(y0[0], y0[1])  # BH temperature
        z    = self.M1/Tp
        zBH  = np.real(self.M1/TBH)

        from ulysses.ulsbase import my_kn2, my_kn1
    
        self._d1    = np.real(self.Gamma1 * my_kn1(z)/my_kn2(z)) # Therm-av RH decay width
        self._dPBH1 = np.real(self.Gamma1 * self.ME(zBH))        # Therm-av RH decay width wrt to TBH -> using full greybody factors
        
        self._w1    = self._d1 * (my_kn2(z) * z**2/(3. * zeta(3.)))
        
        # RH neutrino equilibrium number density, normalized to initial photon density
        self._n1eq  = (10**(3*(x + xilog10)) * self.M1**2 * Tp * my_kn2(z))/(np.pi**2)/nphi 
        
        # Neutrino masses squared
        m1sq = self.SqrtDm[0,0]**4
        m2sq = self.SqrtDm[1,1]**4
        m3sq = self.SqrtDm[2,2]**4

        # Lepton number density in equilibrium, 2 factor corresponds to the number of degrees of freedom
        nleq = (3./4.) * 2. * (zeta(3)/np.pi**2) * Tp**3

        # Thermally averaged scattering 
        gD2 = (3.*Tp**6/(4.*np.pi**5*self.v**4))*(m1sq + m2sq + m3sq)

        # Washout term for the scattering DL = 2 term
        WashDL2 = gD2/nleq
            
        return FBEqs(x, y0, nphi, self.M1,self.M2,self.M3, eps, self._d1, self._w1, self._n1eq, self._dPBH1, np.real(WashDL2), xilog10)

    #******************************************************#
    #        Equations After PBH evaporation               #
    #******************************************************#
    def RHS_aBE(self, x, y0, eps, nphi, MBHi, asi):
        
        Tp   = y0[1]                                   # Plasma Temperature
        TBH  = bh.TBH(MBHi, asi) # BH final temperature
        k    = np.real(self.k1)
        z    = np.real(self.M1/Tp)
        zBH  = np.real(self.M1/TBH)

        from ulysses.ulsbase import my_kn2, my_kn1
   
        self._d1    = np.real(self.Gamma1 * my_kn1(z)/my_kn2(z)) # Therm-av RH decay width
        self._dPBH1 = np.real(self.Gamma1 * self.ME(zBH))        # Therm-av RH decay width wrt to TBH -> using full greybody factors
        
        self._w1    = self._d1 * (my_kn2(z) * z**2/(3. * zeta(3.)))

        # RH neutrino equilibrium number density, normalized to initial photon density
        self._n1eq  = (10**(3*x) * self.M1**2 * Tp * my_kn2(z))/(np.pi**2)/nphi
        
        # Neutrino masses squared
        m1sq = self.SqrtDm[0,0]**4
        m2sq = self.SqrtDm[1,1]**4
        m3sq = self.SqrtDm[2,2]**4

        # Lepton number density in equilibrium, 2 factor corresponds to the number of degrees of freedom
        nleq = (3./4.) * 2. * (zeta(3)/np.pi**2) * Tp**3

        # Thermally averaged scattering 
        gD2 = (3.*Tp**6/(4.*np.pi**5*self.v**4))*(m1sq + m2sq + m3sq)

        # Washout term for the scattering DL = 2 term
        WashDL2 = gD2/nleq

        return FBEqs_aBE(x, y0, nphi, self.M1,self.M2,self.M3, eps, self._d1, self._w1, self._n1eq, self._dPBH1, np.real(WashDL2))

    #******************************************************#
    #                     Main Program                     #
    #******************************************************#

    @property
    def EtaB(self):
        
        #Define fixed quantities for BEs
        eps1tt = np.real(self.epsilon1ab(2,2))
        eps1mm = np.real(self.epsilon1ab(1,1))
        eps1ee = np.real(self.epsilon1ab(0,0))
        eps1tm =         self.epsilon1ab(2,1)
        eps1te =         self.epsilon1ab(2,0)
        eps1me =         self.epsilon1ab(1,0)

        eps = [eps1tt, eps1mm, eps1ee, eps1tm, eps1te, eps1me] # Array for CP violation elements
        
        Mi    = 10**(self.MPBHi) # PBH initial Mass in grams
        asi   = self.aPBHi       # PBH initial rotation a_star factor
        bi    = 10**(self.bPBHi) # Reduced Initial PBH fraction, beta^prime

        assert 0. <= asi and asi < 1., colored('initial spin factor a* is not in the range [0., 1.)', 'red')
        assert bi < np.sqrt(bh.gamma), colored('initial PBH density is larger than the total Universe\'s budget', 'red')

        Ti    = ((45./(16.*106.75*(np.pi*bh.GCF)**3.))**0.25) * np.sqrt(bh.gamma * bh.GeV_in_g/Mi) # Initial Universe temperature in GeV

        rRadi = (np.pi**2./30.) * bh.gstar(Ti) * Ti**4  # Initial radiation energy density -- assuming a radiation dominated Universe
        rPBHi = (bi/(np.sqrt(bh.gamma) -  bi))*rRadi    # Initial PBH energy density, assuming collapse of density fluctuations
        nphi  = (2.*zeta(3)/np.pi**2)*Ti**3             # Initial photon number density

        N1RTi = 0.  # Initial condition for Thermal RH neutrino number density
        N1RBi = 0.  # Initial condition for PBH-emitted RH neutrino number
        NBLi  = 0.  # Initial condition for B-L asymmetry
        
        #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
        #                                           Solving the equations                                                   #
        #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

        Min  = Mi   # We save the initial PBH mass.
        asin = asi  # We save the initial PBH spin
        
        xilog10 = 0. # Fixing the initial scale factor to be 1

        # Defining arrays to save the solution for the different components
        
        xBE     = [] # Log10 of the scale factor
        MBHBE   = [] # PBH mass
        astBE   = [] # PBH spin
        RadBE   = [] # Radiation energy density
        PBHBE   = [] # PBH energy density
        TBE     = [] # Plasma temperature
        N1RTBE  = [] # Thermal RH neutrino number density
        N1RBBE  = [] # PBH-emitted RH neutrino number density
        NBLBE   = [] # B-L number densities

        """
        The PBH evolution is done iteratively since the solver cannot reach the Planck Mass directly.
        We evolve the Friedmann & Boltzmann equations from the initial mass to 1% of such initial mass.
        Then, we use the solutions as initial conditions for a new iteraction.
        This is repeated until the PBH mass reaches the Planck mass.
        """

        i = 0

        while Mi > 1.125*bh.MPL: # We solve the equations until the PBH mass is larger than Planck Mass

            #----------------------------------------------------------------------#
            #     Compute BH lifetime and scale factor in which PBHs evaporate     #
            #----------------------------------------------------------------------#
    
            tau_sol = solve_ivp(fun=lambda t, y: bh.ItauRH(t, y, self.M1, self.M2, self.M3), t_span = [-80, 40.], y0 = [Mi, asi], 
                                  rtol=1.e-5, atol=1.e-20, dense_output=True)

            if i == 0: # For the first iteration, we determine the PBH lifetime.
                tau = tau_sol.t[-1] # Log10@PBH lifetime in inverse GeV

            # We compute the Log10@scale factor, xflog10, when the PBH evaporation happens
            if bi > 1.e-19*(1.e9/Mi): # If the initial PBH density is large enough, we include all components
                xf = root(bh.afin, [40.], args = (rPBHi, rRadi, 10.**tau, 0.), method='lm', tol=1.e-40) # Scale factor
                xflog10 = xf.x[0]
                
            else: # If the initial PBH density is negligible, we consider a radiation dominated Universe to obtain xflog10
                xfw = np.sqrt(1. + 4.*10.**tau*np.sqrt(2.*np.pi*bh.GCF*rRadi/3.))
                xflog10 = np.log10(xfw)
                
            #------------------------------------------------------------#
            #                                                            #
            #               Equations Before BH evaporation              #
            #                                                            #
            #------------------------------------------------------------#

            StopM = lambda t, x:self.StopMass(t, x, Mi) # Event to stop when the mass is 1% of the initial mass
            StopM.terminal  = True
            StopM.direction = -1.

            y0 = [Mi, asi, rRadi, rPBHi, Ti, N1RTi, N1RBi, NBLi] # Initial condition

            # Solving Equations
            solFBE = solve_ivp(lambda t, z: self.RHS(t, z, eps, np.real(nphi), xilog10),
                               [0., xflog10], y0, method='BDF', events=StopM, rtol=1.e-7, atol=1.e-10)

            assert solFBE.t[-1] > 0., colored('Solution going backwards...', 'red')

            # Appending solutions to predefined arrays
            
            xBE    = np.append(xBE,    solFBE.t[:] + xilog10)
            MBHBE  = np.append(MBHBE,  solFBE.y[0,:])
            astBE  = np.append(astBE,  solFBE.y[1,:])
            RadBE  = np.append(RadBE,  solFBE.y[2,:])
            PBHBE  = np.append(PBHBE,  solFBE.y[3,:])
            TBE    = np.append(TBE,    solFBE.y[4,:])
            
            N1RTBE = np.append(N1RTBE, solFBE.y[5,:])
            N1RBBE = np.append(N1RBBE, solFBE.y[6,:])            
            NBLBE  = np.append(NBLBE, solFBE.y[7,:])

            # Updating the initial conditions for next iteration

            Mi    = solFBE.y[0,-1]
            asi   = solFBE.y[1,-1]
            rRadi = solFBE.y[2,-1]
            rPBHi = solFBE.y[3,-1]
            Ti    = solFBE.y[4,-1]
            N1RTi = solFBE.y[5,-1]
            N1RBi = solFBE.y[6,-1]
            
            NBLi  = solFBE.y[7,-1]
           
            xilog10 += solFBE.t[-1]

            i += 1

            assert i < 100, print(colored('Loop is stuck!', 'red'), Mi, bi)

        else:
            xflog10 = xilog10  # We update the value of log10(scale factor) at which PBHs evaporate
            
        #------------------------------------------------------------#
        #                                                            #
        #                Solution after BH evaporation               #
        #                                                            #
        #------------------------------------------------------------#

        """
        If Thermal Leptogenesis occurs after PBH evaporation, we solve the second set of equations with
        initial conditions from the last evolution step, and then we join the solutions
        """

        # We determine the Log10@ scale factor for maximum value of z = M/T after PBH evaporation
        xzmax = xflog10 + np.log10(self._zmax*TBE[-1]/self.M1) 
        
        if xflog10 < xzmax:

            # Initial condition for second set of equations, taken from last evolution step
            y0_aBE = [RadBE[-1], TBE[-1], N1RTBE[-1], N1RBBE[-1], NBLBE[-1]] 
            
            # Solving Equations
            solFBE_aBE = solve_ivp(lambda t, z: self.RHS_aBE(t, z, eps, np.real(nphi), np.real(Min), asin),
                                   [xflog10, xzmax], y0_aBE, method='BDF', rtol=1.e-7, atol=1.e-15)

            npaf = solFBE_aBE.t.shape[0] # Dimension of solution array from Eqs solution

            #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
            #       Joining the solutions before and after evaporation       #
            #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
            
            xBE    = np.append(xBE, solFBE_aBE.t[:])
            
            MBHBE  = np.append(MBHBE,  np.full(npaf, solFBE.y[0,-1]))
            astBE  = np.append(astBE,  np.full(npaf, solFBE.y[1,-1]))
            RadBE  = np.append(RadBE,  solFBE_aBE.y[0,:])
            PBHBE  = np.append(PBHBE,  np.zeros(npaf))
            TBE    = np.append(TBE,    solFBE_aBE.y[1,:])
            
            N1RTBE = np.append(N1RTBE, solFBE_aBE.y[2,:])
            N1RBBE = np.append(N1RBBE, solFBE_aBE.y[3,:])            
            NBLBE  = np.append(NBLBE, solFBE_aBE.y[4,:])
       
        #------------------------------------------------------------#
        #                                                            #
        #                     Conversion to eta_B                    #
        #                                                            #
        #------------------------------------------------------------#

        gstarSrec = bh.gstarS(0.3e-9) # d.o.f. at recombination
        gstarSoff = bh.gstarS(TBE[-1])  # d.o.f. at the end of leptogenesis
     
        SMspl       = 28./79.
        zeta3       = zeta(3)
        ggamma      = 2.
        coeffNgamma = ggamma*zeta3/np.pi**2
        Ngamma      = coeffNgamma*(10**xBE*TBE)**3
        coeffsph    = (SMspl * gstarSrec)/(gstarSoff * Ngamma)

        nb = coeffsph * NBLBE * nphi

        etaB = nb[-1]

        dat = np.array([10.**xBE, np.real(N1RTBE), np.real(N1RBBE), np.real(nb)])

        dat = dat.T
                        
        self.setEvolDataPBH(dat)
        return etaB
