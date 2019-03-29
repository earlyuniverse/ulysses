import leptocalc
import numpy as np
from odeintw import odeintw

class EtaB_1DS_Approx(leptocalc.LeptoCalc):

    def RHS(self, y0,z,epstt,epsmm,epsee,c1t,c1m,c1e,k):
        N1      = y0[0]
        Ntt     = y0[1]
        Nmm     = y0[2]
        Nee     = y0[3]

        d    = np.real(self.D1(k,z))
        w1   = np.real(self.W1(k,z))
        n1eq = self.N1Eq(z)

        c1tc    = np.conjugate(c1t)
        c1mc    = np.conjugate(c1m)
        c1ec    = np.conjugate(c1e)


        #define the different RHSs for each equation
        rhs1 =      -d*(N1-n1eq)

        rhs2 = epstt*d*(N1-n1eq)-0.5*w1*(2*c1t*c1tc*Ntt)
        rhs3 = epsmm*d*(N1-n1eq)-0.5*w1*(2*c1m*c1mc*Nmm)
        rhs4 = epsee*d*(N1-n1eq)-0.5*w1*(2*c1e*c1ec*Nee)

        return [rhs1, rhs2, rhs3, rhs4]

    @property
    def EtaB(self):
        #Define fixed quantities for BEs
        epstt = np.real(self.epsilonab(2,2))
        epsmm = np.real(self.epsilonab(1,1))
        epsee = np.real(self.epsilonab(0,0))

        c1t   =                 self.c1a(2)
        c1m   =                 self.c1a(1)
        c1e   =                 self.c1a(0)

        xs      = np.linspace(self.xmin, self.xmax, self.xsteps)
        k       = np.real(self.k1)
        y0      = np.array([0+0j,0+0j,0+0j,0+0j], dtype=np.complex128)

        params  = np.array([epstt,epsmm,epsee,c1t,c1m,c1e,k], dtype=np.complex128)

        ys      = odeintw(self.RHS, y0, self.xs, args = tuple(params))
        nb      = self.sphalfact*(ys[-1,1]+ys[-1,2]+ys[-1,3])

        return np.real(nb)

    def __call__(self, x):
        """
        Operator that returns etaB.

        x --- parameter dictionary
        """
        self.setParams(x)
        return self.EtaB


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
    ETA = EtaB_1DS_Approx()
    print(ETA(pars))

    import leptomts
    L=leptomts.LeptoCalc(approx=True)
    L.setParams(pars)
    print("Previous code gives etab = ",np.real(L.EtaB))
