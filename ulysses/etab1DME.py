# non-resonant leptogenesis with one decaying sterile neutrino using the density matrix equations. Equations from 1112.4528
import ulysses
import numpy as np
from odeintw import odeintw

from ulysses.numba import jit
@jit
def fast_RHS(y0, d, w1, n1eq, epstt,epsmm,epsee,epstm,epste,epsme,c1t,c1m,c1e, widtht, widthm):

    N1      = y0[0]
    Ntt     = y0[1]
    Nmm     = y0[2]
    Nee     = y0[3]
    Ntm     = y0[4]
    Nte     = y0[5]
    Nme     = y0[6]
    c1tc    = c1t.conjugate()
    c1mc    = c1m.conjugate()
    c1ec    = c1e.conjugate()


    #define the different RHSs for each equation
    rhs1 =      -d*(N1-n1eq)

    rhs2 = epstt*d*(N1-n1eq)-0.5*w1*(2*c1t*c1tc*Ntt + c1m*c1tc*Ntm + c1e*c1tc*Nte + (c1m*c1tc*Ntm+c1e*c1tc*Nte).conjugate()                  )
    rhs3 = epsmm*d*(N1-n1eq)-0.5*w1*(2*c1m*c1mc*Nmm + c1m*c1tc*Ntm + c1e*c1mc*Nme + (c1m*c1tc*Ntm+c1e*c1mc*Nme).conjugate()                  )
    rhs4 = epsee*d*(N1-n1eq)-0.5*w1*(2*c1e*c1ec*Nee + c1e*c1mc*Nme + c1e*c1tc*Nte + (c1e*c1mc*Nme+c1e*c1tc*Nte).conjugate()                  )
    rhs5 = epstm*d*(N1-n1eq)-0.5*w1*(  c1t*c1mc*Nmm + c1e*c1mc*Nte + c1m*c1mc*Ntm + c1mc*c1t*Ntt + c1t*c1tc*Ntm + c1t*c1ec*(Nme.conjugate()) ) - widtht*Ntm - widthm*Ntm
    rhs6 = epste*d*(N1-n1eq)-0.5*w1*(  c1t*c1ec*Nee + c1e*c1ec*Nte + c1m*c1ec*Ntm + c1t*c1ec*Ntt + c1t*c1mc*Nme + c1t*c1tc*Nte               ) - widtht*Nte
    rhs7 = epsme*d*(N1-n1eq)-0.5*w1*(  c1m*c1ec*Nee + c1e*c1ec*Nme + c1m*c1ec*Nmm + c1t*c1ec*(Ntm.conjugate())  + c1m*c1mc*Nme + c1m*c1tc*Nte) - widthm*Nme

    return [rhs1, rhs2, rhs3, rhs4, rhs5, rhs6, rhs7]

class EtaB_1DME(ulysses.ULSBase):
    """
    Density matrix equation (DME) with one decaying sterile. See arxiv:1112.4528.
    """

    def shortname(self): return "1DME"

    def flavourindices(self): return [1, 2, 3]

    def flavourlabels(self): return ["$N_{\\tau\\tau}$", "$N_{\mu\mu}$", "$N_{ee}$"]

    def RHS(self, y0,z,epstt,epsmm,epsee,epstm,epste,epsme,c1t,c1m,c1e,k):

        if z != self._currz or z == self.zmin:
            n1,d1,w1 = self.NDW1(k,z)
            self._d       = d1.real
            self._w1      = w1.real
            self._n1eq    = n1

            self._currz=z

        # thermal widths are set to zero such that we are in the "one-flavoured regime"
        widtht = 485e-10*self.MP/self.M1
        widthm = 1.7e-10*self.MP/self.M1
        return fast_RHS(y0,self._d, self._w1, self._n1eq,epstt,epsmm,epsee,epstm,epste,epsme,c1t,c1m,c1e, widtht, widthm)


    @property
    def EtaB(self):
        #Define fixed quantities for BEs
        epstt = np.real(self.epsilon1ab(2,2))
        epsmm = np.real(self.epsilon1ab(1,1))
        epsee = np.real(self.epsilon1ab(0,0))
        epstm =         self.epsilon1ab(2,1)
        epste =         self.epsilon1ab(2,0)
        epsme =         self.epsilon1ab(1,0)

        c1t   =                 self.c1a(2)
        c1m   =                 self.c1a(1)
        c1e   =                 self.c1a(0)

        k       = np.real(self.k1)
        y0      = np.array([0+0j,0+0j,0+0j,0+0j,0+0j,0+0j,0+0j], dtype=np.complex128)

        params  = np.array([epstt,epsmm,epsee,epstm,epste,epsme,c1t,c1m,c1e,k], dtype=np.complex128)

        ys, _      = odeintw(self.RHS, y0, self.zs, args = tuple(params), full_output=True)

        self.setEvolData(ys)

        return self.ys[-1][-1]
