from scipy.special import kn
import numpy as np
from odeintw import odeintw
import timeit
from numpy import linalg

from numba import jit

@jit
def fast_isPerturbative(y):
    """
    Check perturbativity of Yukawas
    """
    #limit of perturbativity for y
    limit = np.power(4*np.pi,0.5)
    #check if any element of column 1, 2 or 3 is larger than limit
    col1               = (y[0,0] < limit)*(y[1,0] < limit)*(y[2,0] < limit)
    col2               = (y[0,1] < limit)*(y[1,1] < limit)*(y[2,1] < limit)
    col3               = (y[0,2] < limit)*(y[1,2] < limit)*(y[2,2] < limit)
    return col1*col2*col3

def my_kn1(x):
    """
    Convenience wrapper for kn(1, x)
    """
    return kn(1, x) if x<=600 else 1e-100#3. + 8.*x

def my_kn2(x):
    """
    Convenience wrapper for kn(2, x)
    """
    return kn(2, x) if x<=600 else 1e-100#15. + 8.*x



# This is the base class
class LeptoCalc(object):
    def __init__(self, *args, **kwargs):
        """
        Set global constants here.
        """
        #Higgs vev, mass and Z-mass in GeV
        self.v = 174.
        self.MH = 125.
        self.MZ = 91.1876
        #relativistic degrees of freedom at high temperature
        self.gstar = 106.75
        #Planck mass in GeV
        self.MP = 1.22e+19
        #neutrino cosmological mass in GeV
        self.mstar = 1.0e-12
        # Flags
        self.debug   = kwargs["debug"]  if kwargs.get("debug")  is not None else False
        # Parameters of the solver
        self._zmin   = kwargs["zmin"]   if kwargs.get("zmin")   is not None else 0.1
        self._zmax   = kwargs["zmax"]   if kwargs.get("zmax")   is not None else 1000
        self._zsteps = kwargs["zsteps"] if kwargs.get("zsteps") is not None else 1000
        self._currz = self.zmin

        # Model switches
        self.ordering = kwargs["ordering"] if kwargs.get("ordering") is not None else 0
        self.loop = kwargs["loop"] if kwargs.get("loop") is not None else False



        self.zs=None
        self.ys=None
        self.setZS()
        self.sphalfact = 0.01

    def __call__(self, x):
        """
        Operator that returns etaB.

        NOTE --- this will work only with derived classes where EtaB is actually implemented

        x --- parameter dictionary
        """
        self.setParams(x)
        return self.EtaB

    def __str__(self):
        s="Model:\n{}".format(self.__doc__)
        s+= "\nNormal ordering\n" if self.ordering==0 else "\nInverted ordering\n"
        s+= "Loop-corrected Yukawa\n" if self.loop else "Tree-level Yukawa\n"
        s+="Integration in [{}, {}] in {} steps\n".format(self._zmin, self._zmax, self._zsteps)
        return s

    def setZMin(self, x):
        self._zmin=x
        self.setZS()

    def setZMax(self, x):
        self._zmax=x
        self.setZS()

    def setZSteps(self, x):
        self._zsteps=x
        self.setZS()

    def setZS(self):
        self.zs = np.geomspace(self.zmin, self.zmax, self.zsteps)
        if self.debug:
            print("Integration range:",self.zs.min(),self.zs.max())

    @property
    def zmin(self): return self._zmin

    @property
    def zmax(self): return self._zmax

    @property
    def zsteps(self): return self._zsteps

    @property
    def evolData(self):
        """
        Return a 4-D array of the evolution data.
        The first column is the evolution variable
        The second column corresponds to Ntautau, the
        third to Nmumu and the last columnd to Nee
        """
        pd = np.empty((self.zsteps, 4))
        pd[:,      0] = self.zs
        pd[:,[1,2,3]] = self.ys
        return pd

    def setParams(self, pdict):
        """
        This set the model parameters. pdict is expected to be a dictionary
        """
        self.delta    = pdict['delta']/180*np.pi
        self.a        = pdict['a21']/180*np.pi
        self.b        = pdict['a31']/180*np.pi
        self.theta12  = pdict['t12']/180*np.pi
        self.theta23  = pdict['t23']/180*np.pi
        self.theta13  = pdict['t13']/180*np.pi
        self.x1       = pdict['x1']/180*np.pi
        self.y1       = pdict['y1']/180*np.pi
        self.x2       = pdict['x2']/180*np.pi
        self.y2       = pdict['y2']/180*np.pi
        self.x3       = pdict['x3']/180*np.pi
        self.y3       = pdict['y3']/180*np.pi
        self.m1       = 10**pdict['m'] * 1e-9 # NOTE input is in log10(m1) in eV --- we convert here to the real value in GeV
        self.M1       = 10**pdict['M1']  #
        self.M2       = 10**pdict['M2']  #
        self.M3       = 10**pdict['M3']  #

    # Some general calculators purely based on input parameters
    @property
    def R(self):
        """
        Orthogonal matrix R = R1.R2.R3
        """

        R1 = np.array([[1., 0., 0.],
                       [0.,  np.cos(self.x1+self.y1*1j), np.sin(self.x1+self.y1*1j)],
                       [0., -np.sin(self.x1+self.y1*1j), np.cos(self.x1+self.y1*1j)]], dtype=np.complex128)

        R2 = np.array([[ np.cos(self.x2+self.y2*1j), 0., np.sin(self.x2+self.y2*1j)],
                       [0., 1. , 0.],
                       [-np.sin(self.x2+self.y2*1j), 0., np.cos(self.x2+self.y2*1j)]], dtype=np.complex128)

        R3 = np.array([[ np.cos(self.x3+self.y3*1j), np.sin(self.x3+self.y3*1j), 0.],
                       [-np.sin(self.x3+self.y3*1j), np.cos(self.x3+self.y3*1j), 0.],
                       [-2., 0., 1.]], dtype=np.complex128)

        return R1 @ R2 @ R3

    @property
    def SqrtDM(self):
        """
        Matrix square root of heavy masses
        """
        return np.array([[np.sqrt(self.M1), 0., 0.],
                         [0., np.sqrt(self.M2), 0.],
                         [0., 0., np.sqrt(self.M3)]], dtype=np.complex128)

    @property
    def DM(self):
        """
        Heavy mass matrix
        """
        return np.array([[self.M1, 0., 0.],
                         [0., self.M2, 0.],
                         [0., 0., self.M3]], dtype=np.complex128)

    @property
    def SqrtDm(self):
        """
        Matrix square root of light masses.
        Everything is in GeV.
        """
        msplit2_solar       =  7.40e-5*1e-18 # 2017
        msplit2_athm_normal = 2.515e-3*1e-18 # Values
        msplit2_athm_invert = 2.483e-3*1e-18 # from nu-fit 3.1

        if self.ordering==0:
            m11 = np.sqrt(self.m1)
            m22 = np.sqrt(np.sqrt(msplit2_solar       + self.m1*self.m1))
            m33 = np.sqrt(np.sqrt(msplit2_athm_normal + self.m1*self.m1))
        elif self.ordering==1:
            m11 = np.sqrt(np.sqrt(msplit2_athm_invert + self.m1*self.m1 - msplit2_solar))
            m22 = np.sqrt(np.sqrt(msplit2_athm_invert + self.m1*self.m1))
            m33 = np.sqrt(self.m1)
        else:
            raise Exception("ordering %i not implemented"%self.ordering)

        return np.array([ [m11,  0.,  0.],
                          [ 0., m22,  0.],
                          [ 0.,  0., m33] ], dtype=np.complex128)

    @property
    def U(self):
        """
        PMNS
        """
        s12     = np.sin(self.theta12)
        s23     = np.sin(self.theta23)
        s13     = np.sin(self.theta13)
        c12     = np.power(1-s12*s12,0.5)
        c23     = np.power(1-s23*s23,0.5)
        c13     = np.power(1-s13*s13,0.5)
        return np.array([ [c12*c13,c13*s12*np.exp(self.a*1j/2.), s13*np.exp(self.b*1j/2-self.delta*1j)],
                           [-c23*s12 - c12*np.exp(self.delta*1j)*s13*s23,np.exp((self.a*1j)/2.)*(c12*c23 - np.exp(self.delta*1j)*s12*s13*s23) , c13*np.exp((self.b*1j)/2.)*s23],
                           [-c12*c23*np.exp(self.delta*1j)*s13 + s12*s23,np.exp((self.a*1j)/2.)*(-c23*np.exp(self.delta*1j)*s12*s13 - c12*s23) ,c13*c23*np.exp((self.b*1j)/2.)]], dtype=np.complex128)


    @property
    def fMR(self):
        """
        This returns the sqrt of the inverse of what is called f(mR) or so
        Equation (3) --- the inverse
        Note: the prefactor differs from the paper due to the fact that in
        this code we use v=174 (i.e.sqrt(2)*vev(higgs))

        """
        return np.sqrt(linalg.inv(self.fMRcore))

    @property
    def fMRcore(self):
        a = 1./self.M1
        b = 1./self.M2
        c = 1./self.M3

        prefactor = -1./(32*np.pi**2*self.v**2)
        d = self.fMLoop(self.M1)
        e = self.fMLoop(self.M2)
        f = self.fMLoop(self.M3)
        A = np.diag([a,b,c])
        B = prefactor*np.diag([d,e,f])
        return A+B

    @property
    def helper(self):
        prefactor = 1./(32*np.pi**2*self.v**2)
        d = self.fMLoop(self.M1)
        e = self.fMLoop(self.M2)
        f = self.fMLoop(self.M3)
        B = prefactor*np.diag([d,e,f])
        return B

    def fMLoop(self, x):
        """
        The loop correction bit in Eq. 3 --- x is M1, M2 or M3
        """
        rH2 = (x/self.MH)**2
        rZ2 = (x/self.MZ)**2
        return x*(np.log(rH2)/(rH2-1) + 3*np.log(rZ2)/(rZ2-1) )

    # TODO check logic here after change of h_loop
    @property
    def isTreeDominant(self):
        min_tree = abs(self.h_tree).min()
        max_loop = abs(self.h_loop).max()
        print("Want (Loop, tree)", max_loop, "to be <", min_tree)
        return max_loop < .1 * min_tree

    # TODO check logic here after change of h_loop
    @property
    def isLoopDominant(self):
        max_tree = abs(self.h_tree).max()
        min_loop = abs(self.h_loop).min()
        print("Want (Loop, tree)", min_loop, "to be >", max_tree)
        return min_loop > 10 * max_tree

    @property
    def m_tree(self):
        """
        Tree-level mass matrix
        """
        return self.v**2 * self.h @ np.linalg.inv(self.DM) @ np.transpose(self.h)

    @property
    def m_tot(self):
        """
        (Tree + loop)  mass matrix
        """
        return self.v**2 * self.h @ self.fMRcore @ np.transpose(self.h)

    @property
    def m_loop(self):
        """
        loop  mass matrix
        """
        return -1 * self.v**2 * self.h @ self.helper @ np.transpose(self.h)

    @property
    def m_2loop(self):
        """
        loop  mass matrix
        """
        return 1./(16*np.pi**2)*(np.max(np.abs(self.h))**2)*self.m_loop

    # @property
    # def h_loop(self):
        # return self.h +  self.h_tree

    @property
    def h(self):
        return self.h_loop if self.loop else self.h_tree

    @property
    def h_loop(self):
        """
        Yukawa matrix (LOOP + Tree)
        """
        return (1./self.v)*(self.U @ self.SqrtDm @ np.transpose(self.R) @ self.fMR)

    @property
    def h_tree(self):
        """
        Yukawa matrix, tree-level
        """
        return (1./self.v)*(self.U @ self.SqrtDm @ np.transpose(self.R) @ self.SqrtDM)

    def FTmeasure(self, debug=False):
        U_tree, S_tree, V_tree = linalg.svd(self.m_tree)
        U_loop, S_loop, V_loop = linalg.svd(self.m_loop)
        U_tot,  S_tot,  V_tot  = linalg.svd(self.m_tot)
        U_2loop,S_2loop,V_2loop= linalg.svd(self.m_2loop)

        meas = sum([abs(x) for x in S_loop])/sum([abs(x) for x in S_tot])

        if debug:
            print("Total of singular values tree:", np.sum(S_tree))
            print("Total of singular values loop:", np.sum(S_loop))
            print("Total of singular values tot:" , np.sum(S_tot))
            print("Total of singular values 2loop:" , np.sum(S_2loop))
            print("measure:", meas)

        return meas

    @property
    def meff1(self):
        """
        Effective mass 1
        """
        return np.dot(np.conjugate(np.transpose(self.h)),self.h)[0,0]*(self.v**2)/self.M1

    @property
    def meff2(self,):
        """
        Effective mass 2
        """
        return np.dot(np.conjugate(np.transpose(self.h)),self.h)[1,1]*(self.v**2)/self.M2

    @property
    def meff3(self,):
        """
        Effective mass 3
        """
        return np.dot(np.conjugate(np.transpose(self.h)),self.h)[2,2]*(self.v**2)/self.M3

    @property
    def k1(self):
        """
        Decay parameter 1
        """
        return self.meff1/self.mstar

    @property
    def k2(self):
        """
        Decay parameter 2
        """
        return self.meff2/self.mstar

    @property
    def k3(self):
        """
        Decay parameter 3
        """
        return self.meff3/self.mstar

    def D1(self, k1, z):
        """
        Decay term for Boltzmann equation
        """
        return  k1*z*my_kn1( z)/my_kn2( z)

    def j(self, z):
        """
        Function that multiplies washouts to incorporate scatterings
        """
        M         = self.DM
        prefactor = 0.1*(1+15/(8*z))
        t         = np.log(M[0,0]/self.MH)
        bracket   = z*t*np.log(1+(1/z)*8.77298/t)+1/z
        return prefactor*bracket

    def DS(self, k1, z):
        """
        Function that replaces decay term when scattering is also present
        """
        M         = self.DM
        t         = np.log(M[0,0]/self.MH)
        DStemp    = 0.1*k1*(1+t*(z**2)*np.log(1+8.77298/(t*z)))
        return DStemp

    def KNR(self, x):
        """
        Test function to see the numerical stability behaviour of the ratio
        of two modified Bessel functions.
        """
        return my_kn1(x)/my_kn2(x)

    def D2(self, k,z):
        """
        Decay term for Boltzmann equation
        """
        r =self.M2/self.M1
        a = r*r
        x=np.real(r*z)
        b = my_kn1(x)
        c = my_kn2(x)
        return k*z*a*b/c

    def D3(self, k,z):
        """
        Decay term for Boltzmann equation
        """
        r =self.M3/self.M1
        a = r*r
        x=np.real(r*z)
        b = my_kn1(x)
        c = my_kn2(x)
        return k*z*a*b/c

    def N1Eq(self, z):
        """
        Equilibrium N1 number density
        """
        n1 = 3./8.*(z**2)*my_kn2(z)
        return n1

    def N2Eq(self, z):
        """
        Equilibrium N2 number density
        For numerical reasons, cut off the return value if there are more than 5 orders of
        magnitude between N1Eq and N2Eq.
        """
        r = self.M2/self.M1
        # n1 = self.N1Eq(z)
        n2 = 3./8.*np.power(r*z,2)*my_kn2(r*z)

        # RATIO = n1/n2
        return n2
        # return n2 if RATIO < 8e40 and RATIO >1e-5  else 0 # NOTE these are ad-hoc magic values

    def N3Eq(self, z):
        """
        Equilibrium N2 number density
        For numerical reasons, cut off the return value if there are more than 5 orders of
        magnitude between N1Eq and N3Eq.
        """
        r = self.M3/self.M1
        # n1 = self.N1Eq(z)
        n3 = 3./8.*np.power(r*z,2)*my_kn2(r*z)

        # RATIO = n1/n3
        return n3
        # return n3 if RATIO < 8e40 and RATIO >1e-5  else 0 # NOTE these are ad-hoc magic values

    def W1(self, k1, z):
        """
        Washout parameter 1
        """
        w1 = 1./4*(z**3)*k1*my_kn1(z)
        return w1

    def W2(self, k, z):
        """
        Washout parameter 2
        """
        w1=self.W1(k,z)
        r = self.M2/self.M1
        w2 = k*r/4*np.power(r*z,3) * my_kn1(r*z)
        return w2
        # RATIO = w1/w2
        # return w2 if RATIO < 8e4 and RATIO >1e-5  else 0 # NOTE these are ad-hoc magic values

    def W2p(self, k, z):
        """
        Washout parameter 2 --- derivative w.r.t. r*z
        """
        # w1=self.W1(k,z)
        r = self.M2/self.M1
        w2p = k*r/4 *(3* np.power(r*z,2)*my_kn1(r*z) + np.power(r*z,3)*kvp(1,r*z))
        return w2p
        # RATIO = w1/w2
        # return w2 if RATIO < 8e4 and RATIO >1e-5  else 0 # NOTE these are ad-hoc magic values

    def W3(self, k, z):
        """
        Washout parameter 3
        """
        r = self.M2/self.M1
        w3 = k*r/4*np.power(r*z,3) * my_kn1(r*z)
        return w3

    def hterm(self, a, b):
        """
        Probability coefficient

        a ... [0,1,2]
              0 = e
              1 = mu
              2 = tau
        b ... [0,1]
              0 = term1
              1 = term2

        """
        norm          = 1./((np.dot(np.conjugate(np.transpose(self.h)), self.h))[0,0])
        return norm*(np.abs(self.h[a,b])**2)

    def c1a(self, a):
        """
        Probability coefficient for 1 a
        """
        norm          = np.sqrt(1./((np.dot(np.conjugate(np.transpose(self.h)), self.h))[0,0]))
        return norm*(self.h[a,0])

    def c2a(self, a):
        """
        Probability coefficient for 1 a
        """
        norm          = np.sqrt(1./((np.dot(np.conjugate(np.transpose(self.h)), self.h))[1,1]))
        return norm*(self.h[a,1])

    def c3a(self, a):
        """
        Probability coefficient for 1 a
        """
        norm          = np.sqrt(1./((np.dot(np.conjugate(np.transpose(self.h)), self.h))[2,2]))
        return norm*(self.h[a,2])

    def f1(self, x):
        """
        f1(x) appears in the expression for epsilon but is numerically unstable so approx. with pw definition
        """
        r2=np.power(x,2)

        f1temp = 2./3.*r2*( (1.+r2) * np.log( (1.+r2) / r2 ) - (2.-r2)/(1.-r2) )

        return f1temp if x<10000 else 1

    def f2(self, x):
        """
        f2(x) appears in the expression for epsilon
        """
        return (2./3.)*(1/(np.power(x,2)-1.))

    def epsilon(self, i, j, k, m):
        """
        CP asymmetry parameter
        """
        l         = self.h
        ldag      = np.conjugate(np.transpose(l))
        lcon      = np.conjugate(l)
        M         = self.DM
        lsquare   = np.dot(ldag,l)

        #define terms of epsilon: prefactor and first term (first), second term (second) etc.
        prefactor   = (3/(16*np.pi))*(1/(lsquare[i,i]))
        first       = np.imag(lsquare[i,j]*l[m,j]*lcon[m,i])*(M[i,i]/M[j,j])*self.f1(M[j,j]/M[i,i])

        second      = np.imag(lsquare[j,i]*l[m,j]*lcon[m,i])*(2./3.)*(1/(np.power(M[j,j]/M[i,i],2)-1.))
        third       = np.imag(lsquare[i,k]*l[m,k]*lcon[m,i])*(M[i,i]/M[k,k])*self.f1(M[k,k]/M[i,i])
        fourth      = np.imag(lsquare[k,i]*l[m,k]*lcon[m,i])*(2./3.)*(1/(np.power(M[k,k]/M[i,i],2)-1.))
        return prefactor*(first+second+third+fourth)

    def epsilonab(self, a, b):
        """
        CP asymmetry parameter. a and b are NOT the model parameters of the same name
        """
        l         = self.h
        ldag      = np.conjugate(np.transpose(l))
        lcon      = np.conjugate(l)
        M         = self.DM
        lsquare   = np.dot(ldag,l)

        #define terms of epsilon: prefactor and first term (first), second term (second) etc.
        prefactor   = (3/(32*np.pi))*(1/(lsquare[0,0]))
        first       = 1j*(lsquare[1,0]*l[a,0]*lcon[b,1]-lsquare[0,1]*l[a,1]*lcon[b,0]) * (M[0,0]/M[1,1])*self.f1(M[1,1]/M[0,0])
        third       = 1j*(lsquare[2,0]*l[a,0]*lcon[b,2]-lsquare[0,2]*l[a,2]*lcon[b,0])*(M[0,0]/M[2,2])*self.f1(M[2,2]/M[0,0])
        second      = 1j*(2./3.)*(1/(np.power(M[1,1]/M[0,0],2)-1.))*(l[a,0]*lcon[b,1]*lsquare[0,1]-lcon[b,0]*l[a,1]*lsquare[1,0])
        fourth      = 1j*(2./3.)*(1/(np.power(M[2,2]/M[0,0],2)-1.))*(l[a,0]*lcon[b,2]*lsquare[0,2]-lcon[b,0]*l[a,2]*lsquare[2,0])
        return prefactor*(first+second+third+fourth)

    def epsilon1ab(self,a,b):
        l         = self.h
        ldag      = np.conjugate(np.transpose(l))
        lcon      = np.conjugate(l)
        M         = self.DM
        lsquare   = np.dot(ldag,l)

        #define terms of epsilon: prefactor and first term (first), second term (second) etc.
        prefactor   = (3/(32*np.pi))*(1/(lsquare[0,0]))
        first       = 1j*(lsquare[1,0]*l[a,0]*lcon[b,1]-lsquare[0,1]*l[a,1]*lcon[b,0])*(M[0,0]/M[1,1])*self.f1(M[1,1]/M[0,0])
        third       = 1j*(lsquare[2,0]*l[a,0]*lcon[b,2]-lsquare[0,2]*l[a,2]*lcon[b,0])*(M[0,0]/M[2,2])*self.f1(M[2,2]/M[0,0])
        second      = 1j*(2./3.)*(1/(np.power(M[1,1]/M[0,0],2)-1.))*(l[a,0]*lcon[b,1]*lsquare[0,1]-lcon[b,0]*l[a,1]*lsquare[1,0])
        fourth      = 1j*(2./3.)*(1/(np.power(M[2,2]/M[0,0],2)-1.))*(l[a,0]*lcon[b,2]*lsquare[0,2]-lcon[b,0]*l[a,2]*lsquare[2,0])
        epsilon1abtemp = prefactor*(first+second+third+fourth)
        return epsilon1abtemp

    #CP asymmetry parameter
    def epsilon2ab(self,a,b):
        l         = self.h
        ldag      = np.conjugate(np.transpose(l))
        lcon      = np.conjugate(l)
        M         = self.DM
        lsquare   = np.dot(ldag,l)

        #define terms of epsilon: prefactor and first term (first), second term (second) etc.
        prefactor   = (3/(32*np.pi))*(1/(lsquare[1,1]))
        first       = 1j*(lsquare[0,1]*l[a,1]*lcon[b,0]-lsquare[1,0]*l[a,0]*lcon[b,1])*(M[1,1]/M[0,0])*self.f1(M[0,0]/M[1,1])
        third       = 1j*(lsquare[2,1]*l[a,1]*lcon[b,2]-lsquare[1,2]*l[a,2]*lcon[b,1])*(M[1,1]/M[2,2])*self.f1(M[2,2]/M[1,1])
        second      = 1j*(2./3.)*(1/(np.power(M[0,0]/M[1,1],2)-1.))*(l[a,1]*lcon[b,0]*lsquare[1,0]-lcon[b,1]*l[a,0]*lsquare[0,1])
        fourth      = 1j*(2./3.)*(1/(np.power(M[2,2]/M[1,1],2)-1.))*(l[a,1]*lcon[b,2]*lsquare[1,2]-lcon[b,1]*l[a,2]*lsquare[2,1])
        epsilon2abtemp = prefactor*(first+second+third+fourth)
        return epsilon2abtemp
        #################################

    def epsilon3ab(self,a,b):
        l         = self.h
        ldag      = np.conjugate(np.transpose(l))
        lcon      = np.conjugate(l)
        M         = self.DM
        lsquare   = np.dot(ldag,l)

        #define terms of epsilon: prefactor and first term (first), second term (second) etc.
        prefactor   = (3/(32*np.pi))*(1/(lsquare[2,2]))
        first       = 1j*(lsquare[0,2]*l[a,2]*lcon[b,0]-lsquare[2,0]*l[a,0]*lcon[b,2])*(M[2,2]/M[0,0])*self.f1(M[0,0]/M[2,2])
        third       = 1j*(lsquare[1,2]*l[a,2]*lcon[b,1]-lsquare[2,1]*l[a,1]*lcon[b,2])*(M[2,2]/M[1,1])*self.f1(M[1,1]/M[2,2])
        second      = 1j*(2./3.)*(1/(np.power(M[0,0]/M[2,2],2)-1.))*(l[a,2]*lcon[b,0]*lsquare[2,0]-lcon[b,2]*l[a,0]*lsquare[0,2])
        fourth      = 1j*(2./3.)*(1/(np.power(M[1,1]/M[2,2],2)-1.))*(l[a,2]*lcon[b,1]*lsquare[2,1]-lcon[b,2]*l[a,1]*lsquare[1,2])
        epsilon3abtemp = prefactor*(first+second+third+fourth)
        return epsilon3abtemp

    @property
    def isPerturbative(self):
        """
        Check perturbativity of Yukawas
        """
        return fast_isPerturbative(self.h)
        y = self.h
        #limit of perturbativity for y
        limit = np.power(4*np.pi,0.5)
        #check if any element of column 1, 2 or 3 is larger than limit
        col1               = (y[0,0] < limit)*(y[1,0] < limit)*(y[2,0] < limit)
        col2               = (y[0,1] < limit)*(y[1,1] < limit)*(y[2,1] < limit)
        col3               = (y[0,2] < limit)*(y[1,2] < limit)*(y[2,2] < limit)
        return col1*col2*col3

    def loopContTest(self):
        """
        Something tree bla loop contribution compare with lightest neutrino mass
        """
        lhs=0
        return lhs < self.m1

    #################################
    #Check we are not in resonance  #
    #################################

    @property
    def eps1(self):
        epsilon100 = np.real(self.epsilon1ab(0,0))
        epsilon111 = np.real(self.epsilon1ab(1,1))
        epsilon122 = np.real(self.epsilon1ab(2,2))
        epstot  = epsilon100+epsilon111+epsilon122
        return epstot

    @property
    def eps1(self):
        epsilon100 = np.real(self.epsilon1ab(0,0))
        epsilon111 = np.real(self.epsilon1ab(1,1))
        epsilon122 = np.real(self.epsilon1ab(2,2))
        epstot  = epsilon100+epsilon111+epsilon122
        return epstot


    @property
    def eps2(self):
        epsilon200 = np.real(self.epsilon2ab(0,0))
        epsilon211 = np.real(self.epsilon2ab(1,1))
        epsilon222 = np.real(self.epsilon2ab(2,2))
        epstot  = epsilon200+epsilon211+epsilon222
        return epstot

    @property
    def eps3(self):
        epsilon300 = np.real(self.epsilon3ab(0,0))
        epsilon311 = np.real(self.epsilon3ab(1,1))
        epsilon322 = np.real(self.epsilon3ab(2,2))
        epstot  = epsilon300+epsilon311+epsilon322
        return epstot

    @property
    def eps100(self):
        return self.epsilon1ab(0,0)

    @property
    def eps111(self):
        return self.epsilon1ab(1,1)

    @property
    def eps122(self):
        return self.epsilon1ab(2,2)

    @property
    def eps200(self):
        return self.epsilon2ab(0,0)

    @property
    def eps211(self):
        return self.epsilon2ab(1,1)

    @property
    def eps222(self):
        return self.epsilon2ab(2,2)

    @property
    def eps300(self):
        return self.epsilon3ab(0,0)

    @property
    def eps311(self):
        return self.epsilon3ab(1,1)

    @property
    def eps322(self):
        return self.epsilon3ab(2,2)

    def resonance(self, z):
        """
        calculate decay rate Gamma in terms of the total epsilon (epstot)
        """
        eps1tau = np.real(self.epsilon(0,1,2,2))
        eps1mu  = np.real(self.epsilon(0,1,2,1))
        eps1e   = np.real(self.epsilon(0,1,2,0))
        epstot  = eps1tau+eps1mu+eps1e
        k       = self.k1
        d       = self.D1(k,z)
        Gamma   = (self.M1**2/self.MP)*np.sqrt(2*np.pi/3)*np.sqrt((np.pi**2)*self.gstar/30)*(1+epstot)*d/z

        #calculate (decay rate)/(mass splitting)
        return Gamma/(M2-M1)

    @property
    def deltaMu2(self):
        M         = self.DM
        r         = np.transpose(self.R)
        rstar     = np.transpose(np.conjugate(r))

        coeff     = (1/self.v**2)*(1/(4*(np.pi**2)))

        deltaMu2  = coeff*(self.SqrtDm @ self.SqrtDm @ r @ M @ M @ M @ rstar)

        return np.trace(deltaMu2)

    @property
    def HiggsNaturalness(self):
        l         = self.h
        ldag      = np.conjugate(np.transpose(l))
        lcon      = np.conjugate(l)
        M         = self.DM
        lsquare   = np.dot(ldag,l)

        coeff     = 2.9*1e7

        paran     = [(M[i,i]*0.05*1e-9)/((self.v**2)*lsquare[i,i]) for i in range(3)]
        bounds    = [coeff*np.power(paran[i],1./3.) for i in range(3)]
        return bounds

    @property
    def Gamma1(self):
        l         = self.h
        ldag      = np.conjugate(np.transpose(l))
        lcon      = np.conjugate(l)
        M         = self.DM
        lsquare   = np.dot(ldag,l)

        return (M[0,0]/(8*np.pi))*lsquare[0,0]

    @property
    def Gamma2(self):
        l         = self.h
        ldag      = np.conjugate(np.transpose(l))
        lcon      = np.conjugate(l)
        M         = self.DM
        lsquare   = np.dot(ldag,l)

        return (M[1,1]/(8*np.pi))*lsquare[1,1]

    @property
    def Gamma3(self):
        l         = self.h
        ldag      = np.conjugate(np.transpose(l))
        lcon      = np.conjugate(l)
        M         = self.DM
        lsquare   = np.dot(ldag,l)

        return (M[2,2]/(8*np.pi))*lsquare[2,2]

    ################################################
    # Define functions for CPV plots               #
    ################################################
    @property
    def effectiveMass(self):
        u = self.U
        if self.ordering == 0:
            m1   = self.SqrtDm[0,0]**2
            m2   = self.SqrtDm[1,1]**2
            m3   = self.SqrtDm[2,2]**2
            meff = np.abs(m1*u[0,0]**2+m2*u[0,1]**2+m3*u[0,2]**2)
            return meff

    @property
    def S1(self):
        u      = self.U
        s1temp = np.imag(np.conj(u[2,0])*u[2,1])
        return s1temp

    @property
    def S2(self):
        u      = self.U
        s2temp = np.imag(np.conj(u[2,1])*u[2,2])
        return s2temp

    @property
    def JCP(self):
        u       = self.U
        JCPtemp = np.imag(u[0,0]*u[1,1]*np.conj(u[0,1])*np.conj(u[1,0]))
        return JCPtemp

   #####   HERE we insert the code for resonant leptogenesis #######
      #CP asymmetry parameter for flavoured resonant leptogenesis from 0705.218averaged for large  washout
    def epsilonaaRES(self,a):
        l         = self.h
        ldag      = np.conjugate(np.transpose(l))
        lcon      = np.conjugate(l)
        M         = self.DM
        lsquare   = np.dot(ldag,l)
        gamma1    = (lsquare[0,0]/(8*np.pi))*M[0,0]
        gamma2    = (lsquare[1,1]/(8*np.pi))*M[1,1]
        DeltaM    = M[1,1]-M[0,0]
        sum1      = np.imag(ldag[0,a]*l[a,1]*lsquare[0,1])+np.imag(ldag[0,a]*l[a,1]*lsquare[1,0])
        epsbar1   = sum1/(lsquare[0,0]*lsquare[1,1])
        fmix      =  -2*(DeltaM/gamma2)/(1+(2*DeltaM/gamma2)**2)
        fosc      = -0.5*(DeltaM/gamma2)/(1+(DeltaM/gamma2)**2)
        epsbar    = -0.5*epsbar1 * (fosc + fmix)
        return epsbar


    def findZmax(self, zmin=0.1, zmax=100, steps=1000):
        k1, k2 = [np.real(self.k1), np.real(self.k2)]
        ZTEST=np.linspace(zmin, zmax,steps)
        W2      = [np.real(self.W2( k2, x)) for x in ZTEST]
        W2p     = [np.real(self.W2p(k2, x)) for x in ZTEST]

        for num, i in enumerate(W2p):
            if i>0:
                break
        return ZTEST[num-10]

    def Control(self):
        k1, k2 = [np.real(self.k1), np.real(self.k2)]

        XXX=self.xs

        D1      = [np.real(self.D1(k1, x)) for x in XXX]
        W1      = [np.real(self.W1(k1, x)) for x in XXX]
        D2      = [np.real(self.D2(k2, x)) for x in XXX]
        W2      = [np.real(self.W2(k2, x)) for x in XXX]
        N1eq    = [np.real(self.N1Eq(x))   for x in XXX]
        N2eq    = [np.real(self.N2Eq(x))   for x in XXX]

        import pylab
        pylab.clf()
        pylab.plot(XXX, D1, label="D1")
        pylab.plot(XXX, D2, label="D2")
        pylab.plot(XXX, W1, label="W1")
        pylab.plot(XXX, W2, label="W2")
        pylab.plot(XXX, N1eq, label="N1eq")
        pylab.plot(XXX, N2eq, label="N2eq")
        pylab.legend()
        pylab.xscale("log")
        pylab.yscale("log")
        pylab.xlim((self.xmin,self.xmax))
        pylab.savefig(self.plotprefix+"control.pdf")



if __name__ == "__main__":
    pass
    # start_time = timeit.default_timer()

    # import sys
    # from leptomts import readConfig
    # _, pdict = readConfig(sys.argv[1])

    # L=LeptoCalc()
    # L.setParams(pdict)
    # print("Eta_B full =",np.real(L.EtaB))
    # L=LeptoCalc(approx=True)
    # L.setParams(pdict)
    # print("Eta_B appx =",np.real(L.EtaB))
    # L=LeptoCalc(nds=2, controlplots=False)#, plotprefix=sys.argv[2])
    # L.setParams(pdict)
    # print("Eta_B full 2DS =",np.real(L.EtaB))
    # L=LeptoCalc(nds=2, approx=True)
    # L.setParams(pdict)
    # print("Eta_B appx 2DS =",np.real(L.EtaB))
    # L=LeptoCalc(nds=3, controlplots=False)
    # L.setParams(pdict)
    # print("Eta_B full 3DS =",np.real(L.EtaB))
    # L=LeptoCalc(nds=5, controlplots=False)
    # L.setParams(pdict)
    # print("Eta_B scattering ('nds' =5)",np.real(L.EtaB))
    # from IPython import embed
    # embed()

