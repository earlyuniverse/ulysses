# Non-resonant leptogenesis with one decaying RHN, neglecting flavour effects.
# Case D4 of arXiv:0907.0205 — kinetic equilibrium dropped, quantum statistics adopted.
# This is the "full" kinetic BE.
import ulysses
import numpy as np
import math
from scipy.integrate import solve_ivp
from numba import njit
from NumbaQuadpack import quadpack_sig, dqags
import numba as nb


##################################################################################
#     Functions to establish the RHS ODEs                                        #
##################################################################################

def D4Nrhs(z, fN, yn, K):
    """RHS of the RHN distribution function ODE for Case 4."""
    en    = np.sqrt(z * z + yn * yn)
    gamma = (z * z * K) / (en * yn)
    delta = (fN - 1. + np.exp(en) * fN) / (np.exp(en) + 1.)
    alpha = 0.5 * (np.exp((en - yn) / 2) - np.exp(-(en - yn) / 2))
    beta  = 0.5 * (np.exp((en + yn) / 2) - np.exp(-(en + yn) / 2))
    return gamma * delta * np.log(alpha / beta)


def solve(yn, K, zs, fN0=None, method="RK45", max_step=1/300.):
    """Solve the RHN distribution function ODE for each yn.

    fN0: initial Fermi-Dirac distribution values at zs[0] (one per yn entry).
         None → vanishing IC (zeros). For thermal IC pass
         1.0 / (np.exp(np.sqrt(zs[0]**2 + yn**2)) + 1).
    """
    y0 = fN0 if fN0 is not None else np.zeros_like(yn)
    sol = solve_ivp(D4Nrhs, [zs[0], zs[-1]],
                    y0, args=(yn, K),
                    max_step=max_step, method=method, dense_output=True)
    return sol.sol(zs)


def Nl_lstore(fn_store, yn_vals, zs):
    """Build and return the RHS of the lepton asymmetry ODE for Case 4."""
    dlogy  = np.log10(yn_vals[1]) - np.log10(yn_vals[0])
    dlogy0 = np.log10(yn_vals[0])
    dlogz  = np.log10(zs[1]) - np.log10(zs[0])
    dlogz0 = np.log10(zs[0])

    @njit
    def logindexy(y):
        """Index in yn_vals just below y."""
        return math.floor((math.log10(y) - dlogy0) / dlogy)

    @njit
    def logindexz(z):
        """Index in zs just below z."""
        return math.floor((math.log10(z) - dlogz0) / dlogz)

    @njit
    def lininterp(y, y1, y2, fn1, fn2):
        """Linear interpolation between (y1, fn1) and (y2, fn2)."""
        return fn1 + (y - y1) * (fn2 - fn1) / (y2 - y1)

    @njit
    def eval(z, y):
        try:
            yindex = logindexy(y)
            zindex = logindexz(z)
        except:
            pass
        if zindex < len(zs) - 1 and zindex >= 0:
            fn_z  = lininterp(z, zs[zindex], zs[zindex+1],
                               fn_store[yindex][zindex], fn_store[yindex][zindex+1])
            fn_z1 = lininterp(z, zs[zindex], zs[zindex+1],
                               fn_store[yindex+1][zindex], fn_store[yindex+1][zindex+1])
        else:
            fn_z = fn_z1 = 0.0
        if yindex < len(yn_vals) - 1 and yindex >= 0:
            return lininterp(y, yn_vals[yindex], yn_vals[yindex + 1], fn_z, fn_z1)
        return 0.0

    @nb.cfunc(quadpack_sig)
    def D4Lintegrand(yn, data_):
        """Integrand for the lepton asymmetry evolution in Case 4."""
        data = nb.carray(data_, (4,))
        yl, z, Nl, eps = data
        en   = math.sqrt(z * z + yn * yn)
        funN = eval(z, abs(yn))
        fleq = math.exp(-yl) / (1. + math.exp(-yl))
        fphi = math.exp(-(en - yl)) / (1. - math.exp(-(en - yl)))
        p1   = (yn / en) * ((fphi + funN) * ((4/3) * Nl + 2 * eps) * fleq)
        p2   = (yn / en) * (-2 * eps * funN * (1 + fphi))
        return p1 + p2
    D4Lintegrand_ptr = D4Lintegrand.address

    @nb.cfunc(quadpack_sig)
    def ynintegral(yl, data_):
        data  = nb.carray(data_, (3,))
        z, Nl, eps = data
        data1 = np.array([yl, z, Nl, eps], np.float64)
        nlowerlim = np.abs((-z * z + 4. * yl * yl) / (4. * yl))
        return dqags(D4Lintegrand_ptr, nlowerlim, 300., data1,
                     epsrel=1e-10, epsabs=1e-10)[0]
    ynintegral_ptr = ynintegral.address

    @njit
    def Nlrhs(z, Nl, eps, K):
        data = np.array([z, Nl[0], eps], np.float64)
        integral = dqags(ynintegral_ptr, 1e-9, 300., data,
                         epsrel=1e-10, epsabs=1e-10)[0]
        # 4.80822761 = 4 * zeta(3)
        return -z * z * K / 4.80822761 * integral

    return Nlrhs


class EtaB_1BE1F_Case4(ulysses.ULSBase):
    """
    Boltzmann equations with one decaying sterile, Case D4 of arXiv:0907.0205.
    Kinetic equilibrium of the RHN is dropped; quantum statistics adopted.
    This is the full kinetic BE. See eqs. 3.32 and 3.34.
    Flavour oscillations not included.
    """

    def shortname(self): return "1BE1F_Case4"

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

        fN0      = self.initial_abundance / (np.exp(np.sqrt(zs[0]**2 + yn_vals**2)) + 1.)
        fn_store = solve(yn_vals, K, zs, fN0=fN0)
        Nlrhs    = Nl_lstore(fn_store, yn_vals, zs)

        sol  = solve_ivp(Nlrhs, [zs[0], zs[-1]], y0=[0], t_eval=zs,
                         args=(eps, K), method="RK45", atol=1e-10, rtol=1e-10,
                         dense_output=True)
        ys   = np.transpose((sol.sol(zs)[0], sol.sol(zs)[0]))
        self.setEvolData(ys)

        return self.ys[-1][-1]
