import ulysses
import numpy as np
from scipy import interpolate
from scipy.integrate import odeint
from scipy.special import zeta
from ulysses.numba import jit


#+++++++++++++++++++++++++++++++++++++++++++++++++#
#             FLRW-Boltzmann Equations            #
#+++++++++++++++++++++++++++++++++++++++++++++++++#

@jit
def fast_RHS(y0, a, rRADi, log10_ain, d, w1, epstt, epsmm, epsee, rnuRda_eq, Del, GCF):
    rnuR         =      y0[0] # RHN
    Tp           =      y0[1] # Temperature
    NBL          =      y0[2] # B-L asymmetry
    H            =      np.sqrt(8 * np.pi * GCF * rRADi * 10.**(4.*(log10_ain - a))/3.) #Hubble parameter
    expression1  =      np.log(10.) * (rnuR -  rnuRda_eq) * d/H
    dTda         =    - (Tp/Del) * np.log(10.)
    drnuRda      =    - expression1
    dNBLa        =      (epstt + epsmm + epsee) * expression1 - np.log(10.) * (w1/H) * NBL
    return [drnuRda, dTda, dNBLa]

class EtaB_1BE1Fsf(ulysses.ULSBase):
    """
    Boltzmann equations with one decaying sterile. For detailed discussions of
    equation derivation see arxiv:1104.2750.  Note these kinetic equations do
    not include off diagonal flavour oscillations.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.GCF   = 6.71862e-39      # Gravitational constant in GeV^-2

        #-------------------------------------#
        #    g*(T) and g*S(T) interpolation   #
        #-------------------------------------#
        import os
        data_dir = os.path.dirname(ulysses.__file__)

        fg  = os.path.join(data_dir, "etab1BE1Fscalefactor_gstar.txt")
        fgS = os.path.join(data_dir, "etab1BE1Fscalefactor_gstarS.txt")

        Dg  = np.loadtxt(fg)
        DgS = np.loadtxt(fgS)

        self.tck_   = interpolate.splrep(Dg[:,0],  Dg[:,1],  s=0)
        self.tckS_  = interpolate.splrep(DgS[:,0], DgS[:,1], s=0)
        self.evolname="$a$"

    def ipol_gstar(self, T):
        return interpolate.splev(T, self.tck_, der = 0)

    def ipol_gstarS(self, T):
        return interpolate.splev(T, self.tckS_, der = 0)

    def ipol_dgstarSdT(self,T):
        return interpolate.splev(T, self.tckS_, der = 1)

    def shortname(self): return "1BE1Fsf"

    def flavourindices(self): return [1, 2]

    def flavourlabels(self): return ["$T$", "$NBL$"]

    def RHS(self, y0, a, epstt, epsmm, epsee, rRADi, log10_ain):
        Tp            = y0[1]
        z             = self.M1/Tp
        from ulysses.ulsbase import my_kn2, my_kn1
        _d            = np.real(self.Gamma1* my_kn1(z) / my_kn2(z))
        _w1           = _d * 0.25 * kn2 * z**2
        rnuRda_eq     = self.N1Eq(z)
        Del           =  1. + Tp * self.ipol_dgstarSdT(Tp)/(3. * self.ipol_gstar(Tp)) # Temperature parameter
        return fast_RHS(y0, a, rRADi, log10_ain, _d, _w1,  epstt, epsmm, epsee, rnuRda_eq, Del, self.GCF)

    @property
    def EtaB(self):
        #Define fixed quantities for BEs
        epstt = np.real(self.epsilon1ab(2,2))
        epsmm = np.real(self.epsilon1ab(1,1))
        epsee = np.real(self.epsilon1ab(0,0))
        Ti      = 100 * self.M1 # initial temp 100 greater than mass N1
        rRadi   = np.pi**2 * self.ipol_gstar(Ti) / 30. * Ti**4 # initial radiation domination rho_RAD = pi^2* gstar(T[0])/30*T^4
        y0      = np.array([0., Ti, 0.])
        nphi    = (2.*zeta(3)/np.pi**2) * Ti**3
        params  = np.array([epstt, epsmm, epsee, np.real(rRadi), 0.])
        aflog10 = 5.0
        t1 = np.linspace(0., aflog10, num=1000, endpoint=True)

        # solve equation
        ys          = odeint(self.RHS, y0, t1, args = tuple(params))
        # functions for converting to etaB using the solution to find temp
        T           = ys[:, 1]
        gstarSrec   = self.ipol_gstarS(0.3e-9) # d.o.f. at recombination
        gstarSoff   = self.ipol_gstarS(T[-1])  # d.o.f. at the end of leptogenesis
        SMspl       = 28./79.
        zeta3       = zeta(3)
        ggamma      = 2.
        coeffNgamma = ggamma*zeta3/np.pi**2
        Ngamma      = coeffNgamma*(10**t1*T)**3
        coeffsph    = SMspl * gstarSrec/gstarSoff
        self.ys = np.empty((len(T), 4))
        self.ys[:,0]=t1
        self.ys[:,1]=T
        self.ys[:,2]=ys[:,2]
        self.ys[:,-1] = coeffsph*( ys[:,2])*nphi/Ngamma
        return self.ys[-1][-1]
