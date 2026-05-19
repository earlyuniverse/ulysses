# Non-resonant leptogenesis with one decaying RHN, neglecting flavour effects.
# Case D3 of arXiv:0907.0205 — kinetic equilibrium assumed, quantum statistics adopted.
import ulysses
import numpy as np
import math
from scipy.integrate import solve_ivp
from scipy.special import kn
from numba import njit
from NumbaQuadpack import quadpack_sig, dqags
import numba as nb


##################################################################################
#     Functions to establish the RHS ODEs                                        #
##################################################################################

@nb.cfunc(quadpack_sig)
def D3Nintegrand(yn, data_):
    """Integrand of the RHN evolution equation for Case 3."""
    data  = nb.carray(data_, (1,))
    z     = data[0]
    en    = math.sqrt(z * z + yn * yn)
    alpha = math.sinh((en - yn) / 2)
    beta  = math.sinh((en + yn) / 2)
    return (yn / en) * (1. / (math.exp(en) + 1.)) * math.log(alpha / beta)
D3Nintegrand_ptr = D3Nintegrand.address


def D3Nrhs(z, Nn, K, lowlim, highlim):
    """RHS of the RHN abundance ODE for Case 3."""
    data     = np.array([z], np.float64)
    integral = dqags(D3Nintegrand_ptr, lowlim, highlim, data,
                     epsrel=1e-10, epsabs=1e-10)[0]
    Neq = (3. / 8.) * z * z * kn(2, z)
    return (K / kn(2, z)) * (Nn - Neq) * integral


def solve(yn, K, zs, N0=0., method="RK45", max_step=1/300.):
    """Solve the RHN abundance ODE.

    N0: initial integrated number density at zs[0].
        0 → vanishing IC. For thermal IC pass N1Eq(zs[0]) = 3/8*zs[0]**2*kn(2,zs[0]).
    """
    sol = solve_ivp(D3Nrhs, [zs[0], zs[-1]], [N0],
                    args=(K, yn[0], yn[-1]),
                    max_step=max_step, method=method, dense_output=True)
    return sol.sol(zs)[0]


def Nl_lstore(fn_store, zs):
    """Build and return the RHS of the lepton asymmetry ODE for Case 3."""
    dlogz  = np.log10(zs[1]) - np.log10(zs[0])
    dlogz0 = np.log10(zs[0])

    @njit
    def logindexz(z):
        """Index in zs just below z."""
        return math.floor((math.log10(z) - dlogz0) / dlogz)

    @njit
    def lininterp(y, y1, y2, fn1, fn2):
        """Linear interpolation between (y1, fn1) and (y2, fn2)."""
        return fn1 + (y - y1) * (fn2 - fn1) / (y2 - y1)

    @njit
    def eval(z):
        try:
            zindex = logindexz(z)
        except:
            pass
        if zindex < len(zs) - 1 and zindex >= 0:
            return lininterp(z, zs[zindex], zs[zindex+1],
                             fn_store[zindex], fn_store[zindex+1])
        return 0.0

    @nb.cfunc(quadpack_sig)
    def D3Lintegrand(yn, data_):
        """Integrand for the lepton asymmetry evolution in Case 3."""
        data = nb.carray(data_, (5,))
        yl, z, Nl, eps, Nn = data
        en   = math.sqrt(z * z + yn * yn)
        fneq = math.exp(-en) / (1. + math.exp(-en))
        fleq = math.exp(-yl) / (1. + math.exp(-yl))
        fphi = math.exp(-(en - yl)) / (1. - math.exp(-(en - yl)))
        p1   = (yn / en) * (fphi + Nn * fneq) * ((4/3) * Nl + 2. * eps) * fleq
        p2   = (yn / en) * (-2 * eps * Nn * fneq * (1 + fphi))
        return p1 + p2
    D3Lintegrand_ptr = D3Lintegrand.address

    @nb.cfunc(quadpack_sig)
    def ynintegral(yl, data_):
        """yn integral for a given z."""
        data  = nb.carray(data_, (4,))
        z, Nl, eps, Nn = data
        data1 = np.array([yl, z, Nl, eps, Nn], np.float64)
        nlowerlim = np.abs((-z * z + 4. * yl * yl) / (4. * yl))
        return dqags(D3Lintegrand_ptr, nlowerlim, 300., data1,
                     epsrel=1e-10, epsabs=1e-10)[0]
    ynintegral_ptr = ynintegral.address

    def Nlrhs(z, Nl, eps, K):
        """Full RHS of the N_{l-l} ODE for Case 3."""
        Neq      = 0.375 * z * z * kn(2, z)
        Nn       = eval(z) / Neq
        data     = np.array([z, Nl[0], eps, Nn], np.float64)
        integral = dqags(ynintegral_ptr, 1e-10, 300., data,
                         epsrel=1e-10, epsabs=1e-10)[0]
        # 4.80822761 = 4 * zeta(3)
        return -z * z * K / 4.80822761 * integral

    return Nlrhs


class EtaB_1BE1F_Case3(ulysses.ULSBase):
    """
    Boltzmann equations with one decaying sterile, Case D3 of arXiv:0907.0205.
    Kinetic equilibrium of the RHN is assumed; quantum statistics adopted.
    See eqs. 3.29 and 3.31. Flavour oscillations not included.
    """

    def shortname(self): return "1BE1F_Case3"

    def flavourindices(self): return [1]

    def flavourlabels(self): return [r"$N_{B-L}$"]

    @property
    def EtaB(self):
        nevals      = 500
        yn_vals     = np.logspace(-3., np.log10(350.), nevals)
        self.z_span = [1e-2, 50.]
        epstt       = np.real(self.epsilon1ab(2, 2))
        epsmm       = np.real(self.epsilon1ab(1, 1))
        epsee       = np.real(self.epsilon1ab(0, 0))
        eps         = epsee + epsmm + epstt
        K           = np.real(self.k1)
        zs          = np.logspace(np.log10(self.z_span[0]), np.log10(self.z_span[1]), nevals)
        self.zs     = zs

        N0       = self.initial_abundance * self.N1Eq(zs[0])
        fn_store = solve(yn_vals, K, zs, N0=N0)
        Nlrhs    = Nl_lstore(fn_store, zs)

        sol  = solve_ivp(Nlrhs, [zs[0], zs[-1]], y0=[0], t_eval=zs,
                         args=(eps, K), method="RK45", atol=1e-10, rtol=1e-10,
                         dense_output=True)
        ys   = np.transpose((sol.sol(zs)[0], sol.sol(zs)[0]))
        self.setEvolData(ys)

        return self.ys[-1][-1]
