# resonant leptogenesis with two sterile neutrinos. Equations from 0705.2183
import ulysses
import numpy as np
from odeintw import odeintw
from ulysses.numba import jit

@jit
def fast_RHS(y0, eps2tt,eps2mm,eps2ee,eps1tt,eps1mm,eps1ee,C,d1,d2,w1,w2,n1eq,n2eq):
    N1, N2, Ntt, Nmm, Nee = y0
    c1t,c1m,c1e,c2t,c2m,c2e = C
    c1tc          = np.conjugate(c1t)
    c1mc          = np.conjugate(c1m)
    c1ec          = np.conjugate(c1e)
    c2tc          = np.conjugate(c2t)
    c2mc          = np.conjugate(c2m)
    c2ec          = np.conjugate(c2e)


    #Very low temperature spectator process factors from https://arxiv.org/pdf/hep-ph/0601084.pdf
    CPhie, CPhimu, CPhitau    = 16./79., 16./79., 16./79.

    Clee, Clemu, Cletau       = 442./711., -32./711., -32./711.

    Clmue, Clmumu, Clmutau    = -32./711., 442./711., -32./711.

    Cltaue, Cltaumu, Cltautau = -32./711., -32./711., 442./711.

    #define the different RHSs for each equation
    rhs1           =      -d1*(N1-n1eq)

    rhs2           =      -d2*(N2-n2eq)

    rhs3           = (  eps1tt * d1 * (N1-n1eq) + eps2tt * d2 * (N2-n2eq)
                                 - w1 * (  c1t * c1tc * ((Cltautau + CPhitau)  * Ntt + (Cltaumu + CPhimu)  * Nmm + (Cltaue + CPhie)  * Nee))
                                 - w2 * (  c2t * c2tc * ((Cltautau + CPhitau)  * Ntt + (Cltaumu + CPhimu)  * Nmm + (Cltaue + CPhie)  * Nee)))

    rhs4           = (  eps1mm * d1 * (N1-n1eq) + eps2mm * d2 * (N2-n2eq)
                                 - w1 * (  c1m * c1mc * ((Clmutau + CPhitau)  * Ntt + (Clmumu + CPhimu)  * Nmm + (Clmue + CPhie)  * Nee))
                                 - w2 * (  c2m * c2mc * ((Clmutau + CPhitau)  * Ntt + (Clmumu + CPhimu)  * Nmm + (Clmue + CPhie)  * Nee)))

    rhs5           = (  eps1ee * d1 * (N1-n1eq) + eps2ee * d2 * (N2-n2eq)
                                 - w1 * (  c1e * c1ec * ((Cletau + CPhitau)  * Ntt + (Clemu + CPhimu)  * Nmm + (Clee + CPhie)  * Nee))
                                 - w2 * (  c2e * c2ec * ((Cletau + CPhitau)  * Ntt + (Clemu + CPhimu)  * Nmm + (Clee + CPhie)  * Nee)))

    RHStemp = [rhs1, rhs2, rhs3, rhs4, rhs5]
    return RHStemp

class EtaB_2RESsp(ulysses.ULSBase):
    """
    Resonant equations with two steriles and three lepton flavours including spectator processes. See arXiv:hep-ph/0601084.
    """

    def shortname(self): return "2RESsp"
    def flavourindices(self): return [2, 3, 4]
    def flavourlabels(self): return ["$N_{\\tau\\tau}$", "$N_{\mu\mu}$", "$N_{ee}$"]

    def RHS(self, y0, zzz, ETA, _C, K):
        k1term,k2term = K
        eps2tt,eps2mm,eps2ee,eps1tt,eps1mm,eps1ee = ETA

        if zzz != self._currz or zzz == self.zmin:
            self._d1            = np.real(self.D1(k1term, zzz))
            self._d2            = np.real(self.D2(k2term, zzz))
            self._w1            = np.real(self.W1(k1term, zzz))
            self._w2            = np.real(self.W2(k2term, zzz))
            self._n1eq          = self.N1Eq(zzz)
            self._n2eq          = self.N2Eq(zzz)
            self._currz=zzz

        from ulysses.numba import List
        C=List()
        [C.append(c) for c in _C]
        return fast_RHS(y0, eps2tt, eps2mm, eps2ee, eps1tt, eps1mm, eps1ee, C,self._d1,self._d2,self._w1,self._w2,self._n1eq,self._n2eq)

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
        y0      = np.array([0+0j,0+0j,0+0j,0+0j,0+0j], dtype=np.complex128)

        _ETA = [
                self.epsiloniaaRES(2,1,0),
                self.epsiloniaaRES(1,1,0),
                self.epsiloniaaRES(0,1,0),
                self.epsiloniaaRES(2,0,1),
                self.epsiloniaaRES(1,0,1),
                self.epsiloniaaRES(0,0,1)
            ]
        _C = [  self.c1a(2), self.c1a(1), self.c1a(0), self.c2a(2), self.c2a(1), self.c2a(0)]
        _K = [np.real(self.k1), np.real(self.k2)]

        ys      = odeintw(self.RHS, y0, self.zs, args = tuple([_ETA, _C, _K]))
        self.setEvolData(ys.real)
        return self.ys[-1][-1]
