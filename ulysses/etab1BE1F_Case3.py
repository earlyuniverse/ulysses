# non-resonant leptogenesis with one decaying right handed neutrino (RHN) and neglecting flavour effects. This is Case D3 of 0907.0205 i.e.  kinetic equilibrium of the RHN has been assumed but quantum statistics (rather than the usual assumption of Maxwell-Boltzmann statistics) has been aodopted.
import ulysses
import numpy as np
import math
from odeintw import odeintw
import matplotlib.pyplot as plt
from numpy.lib.function_base import meshgrid
from scipy.integrate import quad, solve_ivp, simpson
from scipy.special import zeta, kn
from numba import njit


##################################################################################
#     functions to establish the RHS ODEs                                       #
##################################################################################




def D3Lintegrand(yn, yl, z, Nl, eps, Nn):
        """Returns the integrand for the lepton asymmetry evolution in Case 3"""
        en = math.sqrt(z * z + yn * yn)
        fneq = math.exp(-en) / (1. + math.exp(-en))
        fleq = math.exp(-yl) / (1. + math.exp(-yl))
        fphi = math.exp(-(en-yl)) / (1. - math.exp(-(en-yl)))
        p1   = (yn / en) * (fphi + Nn * fneq) * ((4/3) * Nl + 2. * eps) * fleq
        p2   =  (yn / en) * (- 2 * eps * Nn * fneq * (1 + fphi))
        return p1 + p2


def ynintegral(yl, z, Nl, eps, Nn,  case, highlim = 300, epsrel = 1e-10, epsabs = 1e-10):
    """Performs the yn integral for a given z"""
    nlowerlim = np.abs((-z * z + 4. * yl * yl) / (4. * yl))
    integrand = D3Lintegrand
    integral = quad(integrand, nlowerlim, highlim, args = (yl, z, Nl, eps, Nn,), epsrel = epsrel, epsabs = epsabs)
    return integral[0]
    


def NLrhs(z, Nl, K, eps,  calc, highlim = 300, epsrel = 1e-10, epsabs = 1e-10):
    """Retrieves solutions of yn integration and integrates over yl.
       Returns the full RHS of N_{l-l} equation for each case"""
    llowerlim = 1e-10
    Neq = 0.375 * z * z * kn(2, z)
    Nn = calc.solD3_.sol(z) / Neq
    integral1 = quad(ynintegral, llowerlim, highlim, args=(z, Nl, eps, Nn,  calc), epsabs=epsabs, epsrel=epsrel)
    int = integral1[0]
    return -z * z * K * int * 1/(4*zeta(3))#(3. / 16.)


def Nneq(zs, z, y):
    """Returns N_N^{eq}"""
    en = np.sqrt(z * z + y * y)
    func = 1 / (np.exp(en) + 1)
    sol = Normalise(func, zs, y)
    return sol


def Normalise(array, y_eval, y):
    """Integrates inputted array over normalised yn phase space"""
    integrand = np.multiply(array, y * y / (2*zeta(3)))
    result = simpson(integrand, x=y_eval, axis=0)
    return result.ravel()


def rhNsol(zs):
    """Retrives the solutions from calc, normalises f_N solutions (cases 3)
    and plots the solutions for N_N(z)"""
    z, y = np.meshgrid(zs, calc.yn_)
    Neq = Nneq(zs, z, y)
    D3array = calc.solD3_.sol(zs)
    solD3 = Normalise(D3array, zs, y)
    return solD3
      
      
def Lsol(z_span, zs, K, eps, method="RK45", atol=1e-10, rtol=1e-10):
    """Solves differential equations for each case to get N_{l-l}(z), plots the absolute value
    against z"""
    solLD3 = solve_ivp(NLrhs, z_span, [0], t_eval=zs,
                           args=(K, eps, 3,), method=method, atol=1e-10,
                           rtol=1e-10, dense_output=True)
                           
    return   solLD3.y[-1][-1]
    
    
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


      
    def D3Nintegrand(self, yn, z):
        """Returns the integrand of rhn evolution for Case 3"""
        en = math.sqrt(z * z + yn * yn)
        alpha = np.sinh((en-yn)/2)#0.5*(math.exp((en-yn)/2)-math.exp(-(en-yn)/2))
        beta = np.sinh((en+yn)/2)#0.5*(math.exp((en+yn)/2)-math.exp(-(en+yn)/2))

        return (yn/en) * (1. / (math.exp(en) + 1.)) * math.log(alpha / beta)
        
    def D3Nrhs(self, z, Nn):
        """Returns the RHS of rhn evolution for Case 3"""

        integral = quad(self.D3Nintegrand, self.lowlim_, self.highlim_, args=(z,))
        Neq = (3. / 8.) * z * z * kn(2, z)
        output = (self.K_ / kn(2,z)) * (Nn - Neq) * integral[0]

        return output

        
    def solve(self, fN0 = [0], method="RK45", max_step=1/300.):
        """Solving the differential equation for each case to get fN(z,yN)"""
        self.solD3_ = solve_ivp(self.D3Nrhs, [self.tMin_, self.tMax_],
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



class EtaB_1BE1F_Case3(ulysses.ULSBase):
    """
    Boltzmann equation (BE) with one decaying sterile Case D2 i.e. dropping assumption of kinetic equilibrium. See arxiv:0907.0205
    Eqns. 3.29 and 3.31.  Note these kinetic equations do not include off diagonal
    flavour oscillations.
    """
    def __init__(self, *args, **kwargs):
        super().__init__(self, *args, **kwargs)
        self.calc = None

    def shortname(self): return "1BE1F_Case3"

    def flavourindices(self): return [1]

    def flavourlabels(self): return ["$NBL$"]

    def RHS(self, y0,z,epstt,epsmm,epsee,k):

        if z != self._currz or z == self.zmin:
            self._d       = np.real(self.D1(k,z))
            self._w1      = np.real(self.W1(k,z))
            self._n1eq    = self.N1Eq(z)
            self._currz=z

        return Lsol(self.z_span, self.zs,self._K, self.eps, method="RK45", atol=1e-10, rtol=1e-10)
        
    @property
    def EtaB(self):
        #Define fixed quantities for BEs
        #Define fixed quantities for BEs
        nevals               =  500
        yn_vals              =  np.logspace(-3., np.log10(350.), nevals)
        self.z_span          =  [1e-1, 100.]
        epstt                =  np.real(self.epsilon1ab(2,2))
        epsmm                =  np.real(self.epsilon1ab(1,1))
        epsee                =  np.real(self.epsilon1ab(0,0))
        self.eps             =  epsee + epsmm + epstt
        self._K              =  np.real(self.k1)
        self.zs          =  np.logspace(np.log10(self.z_span[0]), np.log10(self.z_span[1]), nevals)
        self.calc            =  Calculator(self._K, yn_vals, tMin=self.z_span[0], tMax=self.z_span[1])
        y0      = np.array([0+0j,0+0j], dtype=np.complex128)
       
        solLD3 = solve_ivp(NLrhs, self.z_span, [0], t_eval=self.zs,
                       args=(self._K, self.eps, self.calc), method="RK45", atol=1e-10,
                       rtol=1e-10, dense_output=True)
        ys    = np.transpose(( solLD3.sol(self.zs)[0],  solLD3.sol(self.zs)[0]))
        self.setEvolData(ys)
        
        #Gives the number of RHNs per comoving volume N_N(z)        
        N_N = self.zs, self.calc.solD3_.sol(self.zs).ravel()
        
        return self.ys[-1][-1]


