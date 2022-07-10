import numpy as np
import math
import matplotlib.pyplot as plt
from numpy.lib.function_base import meshgrid
from scipy.integrate import quad, solve_ivp, simpson, trapezoid
from scipy.special import zeta, kn
from numba import njit
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import rc
import matplotlib.patches as mpatches
import matplotlib.lines as mlines

rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
font = {'family' : 'normal',
    'weight' : 'bold',
        'size'   : 25}
rc('font', **font)

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

    def D2Nrhs(self, z, fN):
        """Returns the RHS of rhn evolution for Case 2"""
        en = np.sqrt(z * z + self.yn_ * self.yn_)
        return ((z * z * self.K_)/en) * (np.exp(-en) - fN)


    def D3Nintegrand(self, yn, z):
        """Returns the integrand of rhn evolution for Case 3"""
        
        en = math.sqrt(z * z + yn * yn)
        alpha = 0.5*(math.exp((en-yn)/2)-math.exp(-(en-yn)/2))
        beta = 0.5*(math.exp((en+yn)/2)-math.exp(-(en+yn)/2))

        return (yn/en) * (1. / (math.exp(en) + 1.)) * math.log(alpha / beta)

    def D3Nrhs(self, z, Nn):
        """Returns the RHS of rhn evolution for Case 3"""

        integral = quad(self.D3Nintegrand, self.lowlim_, self.highlim_, args=(z,))
        Neq = (3. / 8.) * z * z * kn(2, z)
        output = (self.K_ / kn(2,z)) * (Nn - Neq) * integral[0] 

        return output

    def D4Nrhs(self, z, fN):
        """Returns the RHS of rhn evolution for Case 4"""
        
        en = np.sqrt(z * z + self.yn_ * self.yn_)
        yn = self.yn_
        gamma = (z * z * self.K_)/(en * yn)
        delta = (fN - 1. + np.exp(en)*fN)/(np.exp(en) + 1.)
        
        A = (en-yn)/2
        B = (en+yn)/2
        alpha = 0.5*(np.exp(A)-np.exp(-A))
        beta = 0.5*(np.exp(B)-np.exp(-B))

        return gamma * delta * np.log(alpha / beta)

    def solve(self, fN0 = [0], method="RK45", max_step=1/300.):
        """Solving the differential equation for each case to get fN(z,yN)"""
        print('D1 - N')
        self.solD1_ = solve_ivp(self.D1Nrhs, [self.tMin_, self.tMax_],
                                fN0, max_step=max_step, method=method,
                                dense_output=True)
        print('D2 - N')
        self.solD2_ = solve_ivp(self.D2Nrhs, [self.tMin_, self.tMax_],
                                self.y0_, max_step=max_step,
                                method=method, dense_output=True)
        print('D3 - N')
        self.solD3_ = solve_ivp(self.D3Nrhs, [self.tMin_, self.tMax_],
                                fN0, max_step=max_step, method=method,
                                dense_output=True)
        print('D4 - N')
        self.solD4_ = solve_ivp(self.D4Nrhs, [self.tMin_, self.tMax_],
                                self.y0_, max_step=max_step,
                                method=method, dense_output=True)
                     
            


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


    def eval(self, z, y, case):
        """Evaluates f_N for a given z and y for cases 2 and 4"""
        if z != self.currz_:
            if case == 2:
                self.fN_ = self.solD2_.sol(z) #(500,)
            else:
                self.fN_ = self.solD4_.sol(z)
            self.currz_ = z
        
        try:
            global yindex
            yindex = self.logindex(y)  
        except:
            pass

        if yindex < np.size(self.yn_) - 1:
            return self.lininterp(y, self.yn_[yindex], self.yn_[yindex + 1], self.fN_[yindex], self.fN_[yindex + 1])
        else:
            return 0.0



def D2Lintegrand(yn, yl, z, Nl, eps, Nn, part):
    """Returns the integrand for the lepton asymmetry evolution in Case 2"""
    en   = math.sqrt(z * z + yn * yn)
    fun  = calc.eval(z, np.abs(yn), 2)
    fNeq = math.exp(-en)

    if part == 1: 
        return (yn/en) * (4/3) * Nl * fNeq
    else: 
        return (yn/en) * (- 2. * eps * (fun - fNeq))


@njit
def D3Lintegrand(yn, yl, z, Nl, eps, Nn, part):
    """Returns the integrand for the lepton asymmetry evolution in Case 3"""

    en = math.sqrt(z * z + yn * yn)
    fneq = math.exp(-en) / (1. + math.exp(-en))
    fleq = math.exp(-yl) / (1. + math.exp(-yl))
    fphi = math.exp(-(en-yl)) / (1. - math.exp(-(en-yl)))

    if part == 1:
        return (yn / en) * (fphi + Nn * fneq) * ((4/3) * Nl + 2. * eps) * fleq
    else:
        return (yn / en) * (- 2 * eps * Nn * fneq * (1 + fphi))

def D4Lintegrand(yn, yl, z, Nl, eps, Nn, part):
    """Returns the integrand for the lepton asymmetry evolution in Case 4"""

    en = math.sqrt(z * z + yn * yn)
    funN = calc.eval(z, abs(yn), 4)
    fleq = math.exp(-yl)/(1. + math.exp(-yl))
    fphi = math.exp(-(en-yl)) / (1. - math.exp(-(en-yl)))

    if part == 1:
        return (yn / en) * ((fphi + funN) * ((4/3) * Nl + 2 * eps) * fleq)
    else:
        return (yn/en) * (-2 * eps * funN * (1 + fphi))


def ynintegral(yl, z, Nl, eps, Nn, part, case, highlim = 300, epsrel = 1e-10, epsabs = 1e-10):
    """Performs the yn integral for a given z"""
    nlowerlim = np.abs((-z * z + 4. * yl * yl) / (4. * yl))
    if case == 2:
        integrand = D2Lintegrand
    elif case == 3:
        integrand = D3Lintegrand
    else:
        integrand = D4Lintegrand
    integral = quad(integrand, nlowerlim, highlim, args = (yl, z, Nl, eps, Nn, part), epsrel = epsrel, epsabs = epsabs)
   
    return integral[0]
    

def NLrhs(z, Nl, K, eps, case, highlim = 300, epsrel = 1e-10, epsabs = 1e-10):
    """Retrieves solutions of yn integration and integrates over yl. 
       Returns the full RHS of N_{l-l} equation for each case"""
    llowerlim = 1e-10    
    if case == 1:
        """Case 1"""
        D = z * K * (kn(1, z) / kn(2, z))
        alpha = calc.solD1_.sol(z) - 0.375 * z * z * kn(2, z)
        W = 0.25 * K * z * z * z * kn(1, z)

        print("Nl - Case 1", "z:", z, "Nl:", Nl)
        return (-W * Nl + eps * D * alpha )* 2
 
    elif case == 2:
        """Case 2"""
        Nn = None

        integral1 = quad(ynintegral, llowerlim, highlim, args=(z, Nl, eps, Nn, 1, case,), epsabs=epsabs, epsrel=epsrel)
        integral2 = quad(ynintegral, llowerlim, highlim, args=(z, Nl, eps, Nn, 2, case,), epsabs=epsabs, epsrel=epsrel)
        int = integral1[0] + integral2[0]

        print("Nl - Case 2", "z:", z, "Nl:", Nl)
        return -z * z * K * int * (3/16)  * 2

    elif case == 3:
        """Case 3 (includes factor of 1/Neq missing from integrand)"""
        Neq = 0.375 * z * z * kn(2, z)
        Nn = calc.solD3_.sol(z) / Neq

        integral1 = quad(ynintegral, llowerlim, highlim, args=(z, Nl, eps, Nn, 1, 3,), epsabs=epsabs, epsrel=epsrel)
        integral2 = quad(ynintegral, llowerlim, highlim, args=(z, Nl, eps, Nn, 2, 3,), epsabs=epsabs, epsrel=epsrel)
        int = integral1[0] + integral2[0]

        print("Nl - Case 3","z:", z, "Nl:", Nl)
        return -z * z * K * int * (3. / 16.)  * 2

    elif case == 4:
        """Case 4"""
        Nn = None

        integral1 = quad(ynintegral, llowerlim, highlim, args=(z, Nl, eps, Nn, 1, 4,), epsabs=epsabs, epsrel=epsrel)
        integral2 = quad(ynintegral, llowerlim, highlim, args=(z, Nl, eps, Nn, 2, 4,), epsabs=epsabs, epsrel=epsrel)
        int = integral1[0] + integral2[0]

        print("Nl - Case 4","z:", z, "Nl:", Nl)
        return -z * z * K * int * (3 / 16)  * 2


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


def rhNsol(z_eval, ax):
    """Retrives the solutions from calc, normalises f_N solutions (cases 2 and 4)
    and plots the solutions for N_N(z)"""
    z, y = np.meshgrid(z_eval, calc.yn_)

    Neq = Nneq(z_eval, z, y)
    solD1 = calc.solD1_.sol(z_eval).ravel()
    D2array = calc.solD2_.sol(z_eval)
    solD2 = Normalise(D2array, z_eval, y)
    solD3 = calc.solD3_.sol(z_eval).ravel()
    D4array = calc.solD4_.sol(z_eval)
    solD4 = Normalise(D4array, z_eval, y)

    ax[0].loglog(z_eval, Neq, color='black', dashes=[2, 2]) 
    ax[0].loglog(z_eval, solD1, color=colors[0]) 
    ax[0].loglog(z_eval, solD2, color=colors[1])
    ax[0].loglog(z_eval, solD3, color=colors[2])
    ax[0].loglog(z_eval, solD4, color=colors[3]) 
    ax[0].loglog(z_eval, Neq, color='black', dashes=[2, 2]) 
    
    ax[0].set(box_aspect=1, xlim=(0.1, 50), ylim=(1e-11, 10))
    ax[0].set_xlabel("z", fontsize=20)
    ax[0].set_ylabel('$N_\mathrm{N}$', fontsize=20)
    ax[0].tick_params(axis='both', labelsize=18)
    ax[0].yaxis.set_ticks_position('both')
    ax[0].legend(['$N_\mathrm{N}^\mathrm{eq}$','Case 1', 'Case 2', 'Case 3', 'Case 4'], fontsize=18,loc='lower left',frameon=False)

    if K == 0.1:
        ax[0].text(0.11,1.5,"(a)",fontdict=dict(fontsize = 24))
        ax[0].text(11,1.5,"K = "+str(K),fontdict=dict(fontsize = 24))
    else:
        ax[0].text(0.11,1.5,"(c)",fontdict=dict(fontsize = 24))
        ax[0].text(12,1.5,"K = "+str(K),fontdict=dict(fontsize = 24))

def Lsol(z_span, z_eval, K, eps, ax, method="RK45", atol=1e-10,
         rtol=1e-10):
    """Solves differential equations for each case to get N_{l-l}(z), plots the absolute value
    against z"""
    solLD1 = solve_ivp(NLrhs, z_span, [0], t_eval=z_eval,
                       args=(K, eps, 1,), method=method, atol=atol,
                       rtol=rtol, dense_output=True)
    solLD2 = solve_ivp(NLrhs, z_span, [0], t_eval=z_eval,
                       args=(K, eps, 2,), method=method, atol=atol,
                       rtol=rtol, dense_output=True)
    solLD3 = solve_ivp(NLrhs, z_span, [0], t_eval=z_eval,
                       args=(K, eps, 3,), method=method, atol=atol,
                       rtol=rtol, dense_output=True)
    solLD4 = solve_ivp(NLrhs, z_span, [0], t_eval=z_eval,
                       args=(K, eps, 4,), method=method, atol=atol,
                       rtol=rtol, dense_output=True)


    ax[1].loglog(z_eval, np.abs(solLD1.sol(z_eval)[0]), color=colors[0]) 
    ax[1].loglog(z_eval, np.abs(solLD2.sol(z_eval)[0]), color=colors[1])
    ax[1].loglog(z_eval, np.abs(solLD3.sol(z_eval)[0]), color=colors[2])
    ax[1].loglog(z_eval, np.abs(solLD4.sol(z_eval)[0]), color=colors[3])
    ax[1].set(box_aspect=1, xlim=(0.1, 50), ylim=(1e-11, 1e-6))
    ax[1].set_xlabel("z", fontsize=20)
    #ax[1].set_ylabel('$\mathrm{N}_\mathrm{N}$', fontsize=20)
    ax[1].tick_params(axis='both', labelsize=18)
    ax[1].yaxis.set_ticks_position('both')

    
    if K == 0.1:
        ax[1].text(0.11,4.5e-7,"(b)",fontdict=dict(fontsize = 24))
        ax[1].text(11,4.5e-7,"K = "+str(K),fontdict=dict(fontsize = 24))
        ax[1].legend(['Case 1', 'Case 2', 'Case 3', 'Case 4'], fontsize=18,loc='lower center',frameon=False)
    else:
        ax[1].text(0.11,4.5e-7,"(d)",fontdict=dict(fontsize = 24))
        ax[1].text(12,4.5e-7,"K = "+str(K),fontdict=dict(fontsize = 24))
        ax[1].legend(['Case 1', 'Case 2', 'Case 3', 'Case 4'], fontsize=18,loc='lower left',frameon=False)
        

# Epsilon (amount of CP violation)
eps = 1e-6  

# Decay parameter
K = 10

# Setting different colours for weak and strong washout regime
colors1 = ['#1876cf','#4ac7e4','#93ca7b','#468b2f'] #blue to green
colors2 = ['#de294a','#f5b6b5','#cfb5e0','#8341bb'] #red to purple
if K < 1:
    colors = colors1
else:
    colors = colors2

# Upper and lower bounds of z and points to evaluate at
z_span =  [0.1, 50.]   
z_eval = np.logspace(np.log10(z_span[0]), np.log10(z_span[1]), 500) 

# Values of yn
yn_vals = np.logspace(-3., np.log10(350.), 500)


fig, ax = plt.subplots(1, 2)

# Calculate N_N equations
calc = Calculator(K, yn_vals, tMin=z_span[0], tMax=z_span[1])
# Plot N_N
rhNsol(z_eval, ax)
# Solve and Plot N_{l-l}
Lsol(z_span, z_eval, K, eps, ax)

plt.savefig('all.pdf')
plt.show()
