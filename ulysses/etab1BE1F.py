# non-resonant leptogenesis with one decaying sterile neutrino using the Boltzmann equations and neglecting flavour effects.
import ulysses
import numpy as np
from odeintw import odeintw


from ulysses.numba import jit
@jit
def fast_RHS(y0, d, w1, n1eq, epstt,epsmm,epsee):
    N1      = y0[0]
    NBL     = y0[1]

    rhs1 =                       -d*(N1-n1eq)

    rhs2 =    (epstt+epsmm+epsee)*d*(N1-n1eq)-w1*NBL

    return [rhs1, rhs2]

class EtaB_1BE1F(ulysses.ULSBase):
    """
    Boltzmann equation (BE) with one decaying sterile. See arxiv:1112.4528
    Eqns. 4 and 5.  Note these kinetic equations do not include off diagonal
    flavour oscillations.
    """

    def shortname(self): return "1BE1F"

    def flavourindices(self): return [1]

    def flavourlabels(self): return ["$NBL$"]

    def RHS(self, y0,z,epstt,epsmm,epsee,k):

        if z != self._currz or z == self.zmin:
            self._d       = np.real(self.D1(k,z))
            self._w1      = np.real(self.W1(k,z))
            self._n1eq    = self.N1Eq(z)
            self._currz=z

        return fast_RHS(y0, self._d, self._w1, self._n1eq, epstt,epsmm,epsee)

    @property
    def EtaB(self):
        #Define fixed quantities for BEs
        epstt = np.real(self.epsilon1ab(2,2))
        epsmm = np.real(self.epsilon1ab(1,1))
        epsee = np.real(self.epsilon1ab(0,0))

        k       = np.real(self.k1)
        y0      = np.array([0+0j,0+0j], dtype=np.complex128)

        params  = np.array([epstt,epsmm,epsee,k], dtype=np.complex128)

        ys      = odeintw(self.RHS, y0, self.zs, args = tuple(params))
        self.setEvolData(ys)

        return self.ys[-1][-1]
