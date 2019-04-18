import leptomts
import numpy as np
from odeintw import odeintw

from numba import jit
@jit
def fast_RHS(y0, d, w1, n1eq, epstt,epsmm,epsee,epstm,epste,epsme,c1t,c1m,c1e):

    N1      = y0[0]
    Ntt     = y0[1]
    Nmm     = y0[2]
    Nee     = y0[3]
    Ntm     = y0[4]
    Nte     = y0[5]
    Nme     = y0[6]
    c1tc    = np.conjugate(c1t)
    c1mc    = np.conjugate(c1m)
    c1ec    = np.conjugate(c1e)

# thermal widths are set to zero such that we are in the "one-flavoured regime"
    widtht = 0.0
    widthm = 0.0

    #define the different RHSs for each equation
    rhs1 =      -d*(N1-n1eq)

    rhs2 = epstt*d*(N1-n1eq)-0.5*w1*(2*c1t*c1tc*Ntt + c1m*c1tc*Ntm + c1e*c1tc*Nte + np.conjugate(c1m*c1tc*Ntm+c1e*c1tc*Nte)                  )
    rhs3 = epsmm*d*(N1-n1eq)-0.5*w1*(2*c1m*c1mc*Nmm + c1m*c1tc*Ntm + c1e*c1mc*Nme + np.conjugate(c1m*c1tc*Ntm+c1e*c1mc*Nme)                  )
    rhs4 = epsee*d*(N1-n1eq)-0.5*w1*(2*c1e*c1ec*Nee + c1e*c1mc*Nme + c1e*c1tc*Nte + np.conjugate(c1e*c1mc*Nme+c1e*c1tc*Nte)                  )
    rhs5 = epstm*d*(N1-n1eq)-0.5*w1*(  c1t*c1mc*Nmm + c1e*c1mc*Nte + c1m*c1mc*Ntm + c1mc*c1t*Ntt + c1t*c1tc*Ntm + c1t*c1ec*np.conjugate(Nme) ) - widtht*Ntm - widthm*Ntm
    rhs6 = epste*d*(N1-n1eq)-0.5*w1*(  c1t*c1ec*Nee + c1e*c1ec*Nte + c1m*c1ec*Ntm + c1t*c1ec*Ntt + c1t*c1mc*Nme + c1t*c1tc*Nte               ) - widtht*Nte
    rhs7 = epsme*d*(N1-n1eq)-0.5*w1*(  c1m*c1ec*Nee + c1e*c1ec*Nme + c1m*c1ec*Nmm + c1t*c1ec*np.conjugate(Ntm)  + c1m*c1mc*Nme + c1m*c1tc*Nte) - widthm*Nme

    return [rhs1, rhs2, rhs3, rhs4, rhs5, rhs6, rhs7]

class EtaB_1DS(leptomts.LeptoCalc):
    """
    density matrix equation (DME) zero thermal width (implies one-flavoured regime (1F)) with one decaying sterile.
    """

    def RHS(self, y0,z,epstt,epsmm,epsee,epstm,epste,epsme,c1t,c1m,c1e,k):

        if z != self._currx or z == self.xmin:
            self._d       = np.real(self.D1(k,z))
            self._w1      = np.real(self.W1(k,z))
            self._n1eq    = self.N1Eq(z)
            self._currx=z

        return fast_RHS(y0,self._d, self._w1, self._n1eq,epstt,epsmm,epsee,epstm,epste,epsme,c1t,c1m,c1e)


    @property
    def EtaB(self):
        #Define fixed quantities for BEs
        epstt = np.real(self.epsilonab(2,2))
        epsmm = np.real(self.epsilonab(1,1))
        epsee = np.real(self.epsilonab(0,0))
        epstm =         self.epsilonab(2,1)
        epste =         self.epsilonab(2,0)
        epsme =         self.epsilonab(1,0)

        c1t   =                 self.c1a(2)
        c1m   =                 self.c1a(1)
        c1e   =                 self.c1a(0)

        xs      = np.linspace(self.xmin, self.xmax, self.xsteps)
        k       = np.real(self.k1)
        y0      = np.array([0+0j,0+0j,0+0j,0+0j,0+0j,0+0j,0+0j], dtype=np.complex128)

        params  = np.array([epstt,epsmm,epsee,epstm,epste,epsme,c1t,c1m,c1e,k], dtype=np.complex128)

        ys, _      = odeintw(self.RHS, y0, self.xs, args = tuple(params), full_output=True)
        nb      = self.sphalfact*(ys[-1,1]+ys[-1,2]+ys[-1,3])

        self.ys  = np.real(ys[:, [1,2,3]])

        return np.real(nb)
