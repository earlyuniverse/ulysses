# resonant leptogenesis with two  sterile neutrinos. Equations from 0705.2183
import ulysses
import numpy as np
from odeintw import odeintw
from numba import jit

@jit
def fast_RHS(y0, epstt,epsmm,epsee,C,d1,w1,n1eq):
    N1, Ntt, Nmm, Nee = y0
    c1t,c1m,c1e = C
    c1tc          = np.conjugate(c1t)
    c1mc          = np.conjugate(c1m)
    c1ec          = np.conjugate(c1e)

    #define the different RHSs for each equation
    rhs1           =      -d1*(N1-n1eq)
    rhs2           = (  2 * epstt * d1 * (N1-n1eq)
                                 -  2 * w1 * (  c1t * c1tc * Ntt))
    rhs3           = (  2 * epsmm * d1 * (N1-n1eq)
                                 -  2 * w1 * (  c1m * c1mc * Nmm))
    rhs4           = (  2 * epsee * d1 * (N1-n1eq)
                                 -  2 * w1 * (  c1e * c1ec * Nee))

    RHStemp = [rhs1, rhs2, rhs3, rhs4]
    return RHStemp

class EtaB_2Resonant(ulysses.ULSBase):
    """
    Resonant equations with two steriles and three lepton flavours. See arxiv:0705.2183.
    """

    def RHS(self, y0, zzz, ETA, C, K):
        k1term,k2term = K
        epstt,epsmm,epsee = ETA

        if zzz != self._currz or zzz == self.zmin:
            self._d1            = np.real(self.D1(k1term, zzz))
            self._w1            = np.real(self.W1(k1term, zzz))
            self._n1eq          = self.N1Eq(zzz)
            self._currz=zzz

        return fast_RHS(y0,epstt,epsmm,epsee, C,self._d1,self._w1,self._n1eq)

    @property
    def EtaB(self):

        _HT = [
            np.real(self.hterm(2,0)),
            np.real(self.hterm(1,0)),
            np.real(self.hterm(0,0)),
            np.real(self.hterm(2,1)),
            np.real(self.hterm(1,1)),
            np.real(self.hterm(0,1))
            ]

        _K      = [np.real(self.k1), np.real(self.k2)]
        y0      = np.array([0+0j,0+0j,0+0j,0+0j], dtype=np.complex128)

        _ETA = [
                self.epsilonaaRES(2),
                self.epsilonaaRES(1),
                self.epsilonaaRES(0)
            ]
        _C = [  self.c1a(2), self.c1a(1), self.c1a(0)]
        _K = [np.real(self.k1), np.real(self.k2)]

        ys      = odeintw(self.RHS, y0, self.zs, args = tuple([_ETA, _C, _K]))
        nb      = self.sphalfact*(ys[-1,1]+ys[-1,2]+ys[-1,3])

        self.ys = np.real(ys[:, [1,2,3]])

        return np.real(nb)
