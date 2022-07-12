# non-resonant leptogenesis with one decaying right handed neutrino (RHN) and neglecting flavour effects. This is Case D1 of 0907.0205 i.e. the assumption of kinetic equilibrium of the RHN has been dropped
import ulysses
import numpy as np
import math
from odeintw import odeintw
import matplotlib.pyplot as plt
from numpy.lib.function_base import meshgrid
from scipy.integrate import quad, solve_ivp, simpson, trapezoid
from scipy.special import zeta, kn
from numba import njit


##################################################################################
#     functions to establish the RHS ODEs                                       #
##################################################################################


class Calculator(object):
    
    def __init__(self, K, yn, tMin=0.1, tMax=50):
        self.y0_ = np.zeros_like(yn)
        self.yn_ = yn
        self.dlogy  = np.log10(self.yn_[1]) - np.log10(self.yn_[0])
        self.dlogy0 = np.log10(self.yn_[0])
        self.tMin_ = tMin
        self.tMax_ = tMax
        self.K_ = K
        self.lowlim_ = self.yn_[0]
        self.highlim_ = self.yn_[-1]
        self.currz_ = None
        self.fN_ = None
        self.solve()


    def D1Nrhs(self, z, Nn):
        """Returns the RHS of right-handed neutrino (rhn) evolution for Case 1"""
        D = z * self.K_ * (kn(1, z) / kn(2, z))
        Neq = 3. * z * z * kn(2, z) / 8.
        
        return - D * (Nn - Neq)


        
    def solve(self, fN0 = [0], method="RK45", max_step=1/300.):
        """Solving the differential equation for each case to get fN(z,yN)"""
        self.solD1_ = solve_ivp(self.D1Nrhs, [self.tMin_, self.tMax_],
                                fN0, max_step=max_step, method=method,
                                dense_output=True)

            

    def logindex(self, y):
        """
        Find the index in the ynvals array that corresponds to point
        just below sought after y-value.
        """
        return math.floor((math.log10(y) - self.dlogy0) / self.dlogy)

    def lininterp(self, y, y1, y2, fn1, fn2):
        """
        Linear approximation in 1D, use two point P1 and P2 with
        y1 < y < y2 to compute fn according to line between P1 and P2.
        """
        return fn1 + (y - y1) * (fn2 - fn1) / (y2 - y1)

  


def NLrhs(z, Nl, K, eps, case, highlim = 300, epsrel = 1e-10, epsabs = 1e-10):
    """Retrieves solutions of yn integration and integrates over yl.
       Returns the full RHS of N_{l-l} equation for each case"""
    llowerlim = 1e-10
    if case == 1:
        """Case 1"""
        D = z * K * (kn(1, z) / kn(2, z))
        alpha = calc.solD1_.sol(z) - 0.375 * z * z * kn(2, z)
        W = 0.25 * K * z * z * z * kn(1, z)
        return (-W * Nl + eps * D * alpha )* 2

def Nneq(z_eval, z, y):
    """Returns N_N^{eq}"""
    en = np.sqrt(z * z + y * y)
    func = 1 / (np.exp(en) + 1)
    sol = Normalise(func, z_eval, y)
    return sol


def Normalise(array, y_eval, y):
    """Integrates inputted array over normalised yn phase space"""
    integrand = np.multiply(array, y * y * (3 / 8))
    result = simpson(integrand, x=y_eval, axis=0)
    return result.ravel()


def rhNsol(z_eval):
    """Retrives the solutions from calc, normalises f_N solutions (cases 2 and 4)
    and plots the solutions for N_N(z)"""
    z, y = np.meshgrid(z_eval, calc.yn_)
    Neq = Nneq(z_eval, z, y)
    solD1 = calc.solD1_.sol(z_eval).ravel()
    return solD1
      
      
def Lsol(z_span, z_eval, K, eps, method="RK45", atol=1e-10, rtol=1e-10):
    """Solves differential equations for each case to get N_{l-l}(z), plots the absolute value
    against z"""
    solLD1 = solve_ivp(NLrhs, z_span, [0], t_eval=z_eval,
                       args=(K, eps, 1,), method=method, atol=atol,
                       rtol=rtol, dense_output=True)
                           
    return   solLD1.y[-1][-1]


# number of points to evaluate yn (three-momentum of RHN normalised to T) and number of z points to be evaluated.
nevals     =  500
# yn integration range
yn_vals    = np.logspace(-3., np.log10(350.), nevals)
# span of ODE in z (mass of RHN normalised to T)
z_span     =  [1e-1,10.]
z_eval     = np.logspace(np.log10(z_span[0]), np.log10(z_span[1]), nevals)

# THIS NEEDS TO BE FIXED HERE 10 IS THE K PARAMETER BUT THIS SHOULD BE SET INTERNALLY
calc = Calculator(2.2778530535805257, yn_vals, tMin=z_span[0], tMax=z_span[1])

class EtaB_1BE1F_Case1(ulysses.ULSBase):
    """
    Boltzmann equation (BE) with one decaying sterile Case D1 i.e. dropping assumption of kinetic equilibrium. See arxiv:0907.0205
    Eqns. 3.25 and 3.27.  Note these kinetic equations do not include off diagonal
    flavour oscillations.
    """

    def shortname(self): return "1BE1FCase2"

    def flavourindices(self): return [1]

    def flavourlabels(self): return ["$NBL$"]

    def RHS(self, y0,z,epstt,epsmm,epsee,k):

        if z != self._currz or z == self.zmin:
            self._d       = np.real(self.D1(k,z))
            self._w1      = np.real(self.W1(k,z))
            self._n1eq    = self.N1Eq(z)
            self._currz=z

        return Lsol(z_span, z_eval, K, eps, method="RK45", atol=1e-10, rtol=1e-10)
        
    @property
    def EtaB(self):
        #Define fixed quantities for BEs
        epstt = np.real(self.epsilon1ab(2,2))
        epsmm = np.real(self.epsilon1ab(1,1))
        epsee = np.real(self.epsilon1ab(0,0))
        eps   = epsee + epsmm + epstt
        K       = np.real(self.k1)
        y0      = np.array([0+0j,0+0j], dtype=np.complex128)
       
        solLD1 = solve_ivp(NLrhs, z_span, [0], t_eval=z_eval,
                       args=(K, eps, 1,), method="RK45", atol=1e-10,
                       rtol=1e-10, dense_output=True)
#        due to difference return structure of return statement we must normalise internally in this code rather than call on line 89 ulsbase.py
        normfact = 0.013
        return   solLD1.y[-1][-1] * normfact



