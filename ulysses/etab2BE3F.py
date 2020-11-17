# non-resonant leptogenesis with two decaying sterile neutrino using the Boltzmann equations. Note these kinetic equations do not include off diagonal flavour oscillations. Equations from 1112.4528
import ulysses
import numpy as np
from odeintw import odeintw

from ulysses.numba import jit

@jit
def fast_RHS(y0,eps1tt,eps1mm,eps1ee,eps1tm,eps1te,eps1me,eps2tt,eps2mm,eps2ee,eps2tm,eps2te,eps2me,d1,d2,w1,w2,n1eq,n2eq,C):
    N1, N2, Ntt, Nmm, Nee = y0
    c1t,c1m,c1e,c2t,c2m,c2e = C
    c1tc    = np.conjugate(c1t)
    c1mc    = np.conjugate(c1m)
    c1ec    = np.conjugate(c1e)

    c2tc    = np.conjugate(c2t)
    c2mc    = np.conjugate(c2m)
    c2ec    = np.conjugate(c2e)

    #define the different RHSs for each equation
    rhs1    =      -d1*(N1-n1eq)
    rhs2    =      -d2*(N2-n2eq)
    rhs3    = eps1tt*d1*(N1-n1eq)+eps2tt*d2*(N2-n2eq)-0.5*w1*(2*c1t*c1tc*Ntt) -0.5*w2*(2*c2t*c2tc*Ntt)
    rhs4    = eps1mm*d1*(N1-n1eq)+eps2mm*d2*(N2-n2eq)-0.5*w1*(2*c1m*c1mc*Nmm) -0.5*w2*(2*c2m*c2mc*Nmm)
    rhs5    = eps1ee*d1*(N1-n1eq)+eps2ee*d2*(N2-n2eq)-0.5*w1*(2*c1e*c1ec*Nee) -0.5*w2*(2*c2e*c2ec*Nee)

    RHStemp = [rhs1, rhs2, rhs3, rhs4, rhs5]
    return RHStemp

class EtaB_2BE3F(ulysses.ULSBase):
    """
    Three-flavoured Boltzmann equation (BE) with two decaying steriles. See
    arxiv:1112.4528.  Note these kinetic equations do not include off diagonal
    flavour oscillations.
    """

    def shortname(self): return "2BE3F"
    def flavourindices(self): return [2, 3, 4]
    def flavourlabels(self): return ["$N_{\\tau\\tau}$", "$N_{\mu\mu}$", "$N_{ee}$"]

    def RHS(self, y0, z, ETA, _C, K):
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

        from ulysses.numba import List
        C=List()
        [C.append(c) for c in _C]

        return fast_RHS(y0,eps1tt,eps1mm,eps1ee,eps1tm,eps1te,eps1me,eps2tt,eps2mm,eps2ee,eps2tm,eps2te,eps2me,self._d1,self._d2,self._w1,self._w2,self._n1eq,self._n2eq, C)



    @property
    def EtaB(self):
        #Define fixed quantities for BEs
        _ETA = [
            np.real(self.epsilon(0,1,2,2)),
            np.real(self.epsilon(0,1,2,1)),
            np.real(self.epsilon(0,1,2,0)),
            np.real(self.epsilon(1,0,2,2)),
            np.real(self.epsilon(1,0,2,1)),
            np.real(self.epsilon(1,0,2,0))
            ]

        _HT = [
            np.real(self.hterm(2,0)),
            np.real(self.hterm(1,0)),
            np.real(self.hterm(0,0)),
            np.real(self.hterm(2,1)),
            np.real(self.hterm(1,1)),
            np.real(self.hterm(0,1))
            ]

        _K      = [np.real(self.k1), np.real(self.k2)]
        y0      = np.array([0+0j,0+0j,0+0j,0+0j,0+0j], dtype=np.complex128)

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
        _C = [  self.c1a(2), self.c1a(1), self.c1a(0),
                self.c2a(2), self.c2a(1), self.c2a(0)]
        _K = [np.real(self.k1), np.real(self.k2)]

        ys      = odeintw(self.RHS, y0, self.zs, args = tuple([_ETA, _C, _K]))
        self.setEvolData(ys)
        return self.ys[-1][-1]
