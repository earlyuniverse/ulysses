# non-resonant leptogenesis with one decaying sterile neutrino using the Boltzmann equations and neglecting flavour effects.
import ulysses
import numpy as np
from odeintw import odeintw

from ulysses.numba import jit

@jit
def fast_RHS(y0,eps1tt,eps1mm,eps1ee,eps1tm,eps1te,eps1me,eps2tt,eps2mm,eps2ee,eps2tm,eps2te,eps2me,d1,d2,w1,w2,n1eq,n2eq):
    N1      = y0[0]
    N2      = y0[1]
    NBL     = y0[2]

    rhs1 =                       -d1*(N1-n1eq)
    rhs2 =                       -d2*(N2-n2eq)
    rhs3 = (eps1tt+eps1mm+eps1ee)*d1*(N1-n1eq)+(eps2tt+eps2mm+eps2ee)*d2*(N2-n2eq)-w1*NBL-w2*NBL

    return [rhs1, rhs2, rhs3]

class EtaB_2BE1F(ulysses.ULSBase):
    """
    Boltzmann equation (BE) with two decaying steriles. See arxiv:1112.4528 Eqns. 4 and 5.
    Note these kinetic equations do not include off diagonal flavour
    oscillations.
    """

    def shortname(self): return "2BE1F"

    def flavourindices(self): return [2]

    def flavourlabels(self): return ["$NBL$"]

    def RHS(self, y0, z, ETA, K):
        eps1tt,eps1mm,eps1ee,eps1tm,eps1te,eps1me,eps2tt,eps2mm,eps2ee,eps2tm,eps2te,eps2me = ETA
        k1term,k2term = K

        if z != self._currz or z == self.zmin:
            self._d1      = np.real(self.D1(k1term, z))
            self._w1      = np.real(self.W1(k1term, z))
            self._d2      = np.real(self.D2(k2term, z))
            self._w2      = np.real(self.W2(k2term, z))
            self._n1eq    = self.N1Eq(z)
            self._n2eq    = self.N2Eq(z)
            self._currz=z


        return fast_RHS(y0,eps1tt,eps1mm,eps1ee,eps1tm,eps1te,eps1me,eps2tt,eps2mm,eps2ee,eps2tm,eps2te,eps2me,self._d1,self._d2,self._w1,self._w2,self._n1eq,self._n2eq)

    @property
    def EtaB(self):
        #Define fixed quantities for BEs

        _HT = [
            np.real(self.hterm(2,0)),
            np.real(self.hterm(1,0)),
            np.real(self.hterm(0,0)),
            np.real(self.hterm(2,1)),
            np.real(self.hterm(1,1)),
            np.real(self.hterm(0,1))
            ]

        _K      = [np.real(self.k1), np.real(self.k2)]
        y0      = np.array([0+0j,0+0j,0+0j], dtype=np.complex128)

        _ETA = [
            np.real(self.epsilon1ab(2,2)),
            np.real(self.epsilon1ab(1,1)),
            np.real(self.epsilon1ab(0,0)),
                    self.epsilon1ab(2,1) ,
                    self.epsilon1ab(2,0) ,
                    self.epsilon1ab(1,0),
            np.real(self.epsilon2ab(2,2)),
            np.real(self.epsilon2ab(1,1)),
            np.real(self.epsilon2ab(0,0)),
                    self.epsilon2ab(2,1) ,
                    self.epsilon2ab(2,0) ,
                    self.epsilon2ab(1,0),
            ]
        _K = [np.real(self.k1), np.real(self.k2)]

        ys      = odeintw(self.RHS, y0, self.zs, args = tuple([_ETA, _K]))
        self.setEvolData(ys)
        return self.ys[-1][-1]
