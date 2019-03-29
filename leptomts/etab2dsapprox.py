import leptocalc
import numpy as np
from odeintw import odeintw

class EtaB_2DS_Approx(leptocalc.LeptoCalc):
    def RHS(self, y0, z, ETA, C, K):
        N1, N2, Ntt, Nmm, Nee = y0
        eps1tt,eps1mm,eps1ee,eps1tm,eps1te,eps1me,eps2tt,eps2mm,eps2ee,eps2tm,eps2te,eps2me = ETA
        c1t,c1m,c1e,c2t,c2m,c2e = C
        k1term,k2term = K

        d1      = np.real(self.D1(k1term, z))
        w1      = np.real(self.W1(k1term, z))
        d2      = np.real(self.D2(k2term, z))
        w2      = np.real(self.W2(k2term, z))
        n1eq    = self.N1Eq(z)
        n2eq    = self.N2Eq(z)

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

        # xs      = np.linspace(self.xmin, self.xmax, self.xsteps)
        _K      = [np.real(self.k1), np.real(self.k2)]
        y0      = np.array([0+0j,0+0j,0+0j,0+0j,0+0j], dtype=np.complex128)

        # KKK
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

        ys      = odeintw(self.RHS, y0, self.xs, args = tuple([_ETA, _C, _K]))
        nb      = self.sphalfact*(ys[-1,2]+ys[-1,3]+ys[-1,4])

        return np.real(nb)

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
    ETA = EtaB_2DS_Approx()
    print(ETA(pars))

    import leptomts
    L=leptomts.LeptoCalc(nds=2,approx=True)
    L.setParams(pars)
    print("Previous code gives etab = ",np.real(L.EtaB))
