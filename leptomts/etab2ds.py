import leptocalc
import numpy as np
from odeintw import odeintw

class EtaB_2DS(leptocalc.LeptoCalc):
    def RHS(self, y0, zzz, ETA, C, K, W):
        N1, N2, Ntt, Nmm, Nee, Ntm, Nte, Nme = y0
        eps1tt,eps1mm,eps1ee,eps1tm,eps1te,eps1me,eps2tt,eps2mm,eps2ee,eps2tm,eps2te,eps2me = ETA
        c1t,c1m,c1e,c2t,c2m,c2e = C
        k1term,k2term = K
        widtht,widthm = W
        d1            = np.real(self.D1(k1term, zzz))
        w1            = np.real(self.W1(k1term, zzz))
        d2            = np.real(self.D2(k2term, zzz))
        w2            = np.real(self.W2(k2term, zzz))
        n2eq          = self.N2Eq(zzz)
        n1eq          = self.N1Eq(zzz)


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

        zcrit   = 1e100
        ys, _      = odeintw(self.RHS, y0, self.xs, args = tuple([_ETA, _C , _K, _W]), full_output=1)
        nb      = np.real(self.sphalfact*(ys[-1,2]+ys[-1,3]+ys[-1,4]))

        return nb


if __name__ == "__main__":
    pars = {
            'delta'  :270,
            'a'      :0,
            'b'      :0,
            'theta23':48.7,
            'theta12':33.63,
            'theta13': 8.52,
            'x1'    :45,
            'y1'    :45,
            'x2'    :45,
            'y2'    :45,
            'x3'    :45,
            'y3'    :45,
            'ordering':0,
            'm1'     :-0.60206,
            'M1'     :8,
            'M2'     :9,
            'M3'     :11
            }
    ETA = EtaB_2DS()
    print(ETA(pars))

    import leptomts
    L=leptomts.LeptoCalc(nds=2)
    L.setParams(pars)
    print("Previous code gives etab = ",np.real(L.EtaB))
