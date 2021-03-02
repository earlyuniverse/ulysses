# non-resonant leptogenesis with two decaying sterile neutrino using the density matrix equations. Equations from 1112.4528
import ulysses
import numpy as np
from odeintw import odeintw

from ulysses.numba import jit
@jit
def fast_RHS(y0, eps1tt,eps1mm,eps1ee,eps1tm,eps1te,eps1me,eps2tt,eps2mm,eps2ee,eps2tm,eps2te,eps2me,C,d1,d2,w1,w2,n1eq,n2eq,W):
    N1, N2, Ntt, Nmm, Nee, Ntm, Nte, Nme = y0
    widtht,widthm = W
    c1t,c1m,c1e,c2t,c2m,c2e = C

    c1tc          = np.conjugate(c1t)
    c1mc          = np.conjugate(c1m)
    c1ec          = np.conjugate(c1e)

    c2tc          = np.conjugate(c2t)
    c2mc          = np.conjugate(c2m)
    c2ec          = np.conjugate(c2e)

    #define the different RHSs for each equation
    rhs1           =      -d1*(N1-n1eq)
    rhs2           =      -d2*(N2-n2eq)
    rhs3           = (eps1tt * d1 * (N1-n1eq) + eps2tt * d2 * (N2-n2eq)
             - 0.5 * w1 * (2 * c1t * c1tc * Ntt + c1m * c1tc * Ntm + c1e * c1tc * Nte
             + np.conjugate(c1m * c1tc * Ntm + c1e * c1tc * Nte))
             - 0.5 * w2 * (2 * c2t * c2tc * Ntt + c2m * c2tc * Ntm + c2e * c2tc * Nte
             + np.conjugate(c2m * c2tc * Ntm + c2e * c2tc * Nte)))

    rhs4           = (eps1mm * d1 * (N1-n1eq) + eps2mm * d2 * (N2-n2eq)
         - 0.5 * w1 * (2 * c1m * c1mc * Nmm + c1m * c1tc * Ntm + c1e * c1mc * Nme
         + np.conjugate(c1m * c1tc * Ntm + c1e * c1mc * Nme))
         - 0.5 * w2 * (2 * c2m * c2mc * Nmm + c2m * c2tc * Ntm + c2e * c2mc * Nme
         + np.conjugate(c2m * c2tc * Ntm + c2e * c2mc * Nme)))

    rhs5           = (eps1ee * d1 * (N1-n1eq) + eps2ee * d2 * (N2-n2eq)
         - 0.5 * w1 * (2 * c1e * c1ec * Nee + c1e * c1mc * Nme + c1e * c1tc * Nte
         + np.conjugate(c1e * c1mc * Nme + c1e * c1tc * Nte))
         - 0.5 * w2 * (2 * c2e * c2ec * Nee + c2e * c2mc * Nme + c2e * c2tc * Nte
         + np.conjugate(c2e * c2mc * Nme + c2e * c2tc * Nte)))

    rhs6           = (eps1tm * d1 * (N1-n1eq) + eps2tm * d2 * (N2-n2eq)
         - 0.5 * w1 * (c1t * c1mc * Nmm + c1e * c1mc * Nte + c1m * c1mc * Ntm + c1mc * c1t * Ntt
         + c1t * c1tc * Ntm + c1t * c1ec * np.conjugate(Nme))
         - 0.5 * w2 * (c2t * c2mc * Nmm + c2e * c2mc * Nte + c2m * c2mc * Ntm + c2mc * c2t * Ntt
         + c2t * c2tc * Ntm + c2t * c2ec * np.conjugate(Nme))
         - widtht * Ntm - widthm * Ntm)

    rhs7           = (eps1te * d1 * (N1-n1eq) + eps2te * d2 * (N2-n2eq)
         - 0.5 * w1 * (c1t * c1ec * Nee + c1e * c1ec * Nte + c1m * c1ec * Ntm
         + c1t * c1ec * Ntt + c1t * c1mc * Nme + c1t * c1tc * Nte)
         - 0.5 * w2 * (c2t * c2ec * Nee + c2e * c2ec * Nte + c2m * c2ec * Ntm
         + c2t * c2ec * Ntt + c2t * c2mc * Nme + c2t * c2tc * Nte)
         - widtht * Nte)

    rhs8           = (eps1me * d1 * (N1-n1eq) + eps2me * d2 * (N2-n2eq)
         - 0.5 * w1 * (c1m * c1ec * Nee + c1e * c1ec * Nme + c1m * c1ec * Nmm
         + c1t * c1ec * np.conjugate(Ntm) + c1m * c1mc * Nme + c1m * c1tc * Nte)
         - 0.5 * w2 * (c2m * c2ec * Nee + c2e * c2ec * Nme + c2m * c2ec * Nmm
         + c2t * c2ec * np.conjugate(Ntm) + c2m * c2mc * Nme + c2m * c2tc * Nte)
         - widthm * Nme)

    RHStemp = [rhs1, rhs2, rhs3, rhs4, rhs5, rhs6, rhs7, rhs8]
    return RHStemp

class EtaB_2DME(ulysses.ULSBase):
    """
    Density matrix equation (DME) with two decaying steriles. See arxiv:1112.4528.
    """

    def shortname(self): return "2DME"
    def flavourindices(self): return [2, 3, 4]
    def flavourlabels(self): return ["$N_{\\tau\\tau}$", "$N_{\mu\mu}$", "$N_{ee}$"]

    def RHS(self, y0, zzz, ETA, _C, K, _W):
        eps1tt,eps1mm,eps1ee,eps1tm,eps1te,eps1me,eps2tt,eps2mm,eps2ee,eps2tm,eps2te,eps2me = ETA
        k1term,k2term = K

        if zzz != self._currz or zzz == self.zmin:
            self._d1      = np.real(self.D1(k1term, zzz))
            self._w1      = np.real(self.W1(k1term, zzz))
            self._d2      = np.real(self.D2(k2term, zzz))
            self._w2      = np.real(self.W2(k2term, zzz))
            self._n1eq    = self.N1Eq(zzz)
            self._n2eq    = self.N2Eq(zzz)
            self._currz=zzz
        from ulysses.numba import List
        C=List()
        W=List()
        [C.append(c) for c in _C]
        [W.append(w) for w in _W]

        return fast_RHS(y0, eps1tt,eps1mm,eps1ee,eps1tm,eps1te,eps1me,eps2tt,eps2mm,eps2ee,eps2tm,eps2te,eps2me, C,self._d1,self._d2,self._w1,self._w2,self._n1eq,self._n2eq,W)

    @property
    def EtaB(self):

        #Define fixed quantities for BEs
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
        _W = [ 485e-10*self.MP/self.M1, 1.7e-10*self.MP/self.M1]

        y0      = np.array([0+0j,0+0j,0+0j,0+0j,0+0j,0+0j,0+0j,0+0j], dtype=np.complex128)

        ys, _      = odeintw(self.RHS, y0, self.zs, args = tuple([_ETA, _C , _K, _W]), full_output=1)

        self.setEvolData(ys)
        return self.ys[-1][-1]
