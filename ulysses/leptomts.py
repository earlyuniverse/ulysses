from scipy.special import kn
import numpy as np
from odeintw import odeintw
import timeit
from numpy import linalg


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



class LeptoOld(object):
    def __init__(self, debug=False, nds=1, approx=False, xmin=1e-1, xmax=100, xsteps=1000, controlplots=False, plotprefix=""):
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
        self.debug=debug
        self.nds=nds
        self.approx=approx
        self._xmin = xmin
        self._xmax = xmax
        self._xsteps=xsteps
        self.controlplots=controlplots
        self.xs=None
        self.setXS()
        self.plotprefix=plotprefix
        self.sphalfact = 0.01


    def setXMin(self, x):
        self._xmin=x
        self.setXS()
    def setXMax(self, x):
        self._xmax=x
        self.setXS()
    def setXSteps(self, x):
        self._xsteps=x
        self.setXS()

    def setXS(self):
        self.xs = np.geomspace(self.xmin, self.xmax, self.xsteps)
        if self.debug:
            print("Integration range:",self.xs.min(),self.xs.max())


    @property
    def xmin(self):
        return self._xmin
    @property
    def xmax(self):
        return self._xmax
    @property
    def xsteps(self):
        return self._xsteps

    def setParams(self, pdict):
        """
        This set the model parameters. pdict is expected to be a dictionary
        """
        self.delta    = pdict['delta']/180*np.pi
        self.a        = pdict['a']/180*np.pi
        self.b        = pdict['b']/180*np.pi
        self.theta12  = pdict['theta12']/180*np.pi
        self.theta23  = pdict['theta23']/180*np.pi
        self.theta13  = pdict['theta13']/180*np.pi
        self.x1       = pdict['x1']/180*np.pi
        self.y1       = pdict['y1']/180*np.pi
        self.x2       = pdict['x2']/180*np.pi
        self.y2       = pdict['y2']/180*np.pi
        self.x3       = pdict['x3']/180*np.pi
        self.y3       = pdict['y3']/180*np.pi
        self.m1       = 10**pdict['m1'] * 1e-9 # NOTE input is in log10(m1) in eV --- we convert here to the real value in GeV
        self.M1       = 10**pdict['M1']  #
        self.M2       = 10**pdict['M2']  #
        self.M3       = 10**pdict['M3']  #
        self.ordering = pdict['ordering']

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

    @property
    def isTreeDominant(self):
        min_tree = abs(self.h_tree).min()
        max_loop = abs(self.h_loop).max()
        print("Want (Loop, tree)", max_loop, "to be <", min_tree)
        return max_loop < .1 * min_tree
    
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

    @property
    def h_loop(self):
        return self.h +  self.h_tree
   
    @property
    def h_tree(self):
        """
        Yukawa matrix, tree-level
        """
        return (1./self.v)*(self.U @ self.SqrtDm @ np.transpose(self.R) @ self.SqrtDM)

    def FTmeasure(self, debug=False):
        # e_loop, v_loop = linalg.eig(self.m_loop)
        # e_tot,   v_tot = linalg.eig(self.m_tot)
        U_tree, S_tree, V_tree = linalg.svd(self.m_tree)
        U_loop, S_loop, V_loop = linalg.svd(self.m_loop)
        U_tot,  S_tot,  V_tot  = linalg.svd(self.m_tot)
        U_2loop,S_2loop,V_2loop= linalg.svd(self.m_2loop)
        
        meas = sum([abs(x) for x in S_loop])/sum([abs(x) for x in S_tot])
    
        if debug:
            #print("Eigenvalues loop:", e_loop)
            #print("Eigenvalues tot:", e_tot)
            print("Total of singular values tree:", np.sum(S_tree))
            print("Total of singular values loop:", np.sum(S_loop))
            print("Total of singular values tot:" , np.sum(S_tot))
            print("Total of singular values 2loop:" , np.sum(S_2loop))
            print("measure:", meas)

        return meas


    @property
    def h(self):
        """
        Yukawa matrix (LOOP + Tree)
        """
        #return (1./self.v)*(self.U @ self.SqrtDm @ np.transpose(self.R) @ self.SqrtDM)
        #return (i1./self.v)*(self.U @ self.SqrtDm @ np.transpose(self.R) @ self.fMR)
        #import cmath
        #j=complex(0,1) 
        return (1./self.v)*(self.U @ self.SqrtDm @ np.transpose(self.R) @ self.fMR)
        #return j*(1./self.v)*(self.U.conjugate() @ self.SqrtDm @ np.transpose(self.R) @ self.fMR)

    ##########################################
    #Define functions for Boltzmann equations#
    ##########################################
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
        # print(a,b,prefactor)
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

    ################################################
    #RHS of ODE for derivative of N1, Ntau, Nmu, Ne#
    ################################################
    def RHS_1DS_DM(self, y0,z,epstt,epsmm,epsee,epstm,epste,epsme,c1t,c1m,c1e,k):
        N1      = y0[0]
        Ntt     = y0[1]
        Nmm     = y0[2]
        Nee     = y0[3]
        Ntm     = y0[4]
        Nte     = y0[5]
        Nme     = y0[6]

        d       = np.real(self.D1(k,z))
        w1      = np.real(self.W1(k,z))
        n1eq    = self.N1Eq(z)

        c1tc    = np.conjugate(c1t)
        c1mc    = np.conjugate(c1m)
        c1ec    = np.conjugate(c1e)

        #widtht  = 485e-10*self.MP/self.M1
        #widthm  = 1.7e-10*self.MP/self.M1
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

    def RHS_1DS_DM_ZeroWidth(self, y0,z,epstt,epsmm,epsee,epstm,epste,epsme,c1t,c1m,c1e,k):
        N1      = y0[0]
        Ntt     = y0[1]
        Nmm     = y0[2]
        Nee     = y0[3]
        Ntm     = y0[4]
        Nte     = y0[5]
        Nme     = y0[6]

        d       = np.real(self.D1(k,z))
        w1      = np.real(self.W1(k,z))
        n1eq    = self.N1Eq(z)

        c1tc    = np.conjugate(c1t)
        c1mc    = np.conjugate(c1m)
        c1ec    = np.conjugate(c1e)

        widtht  = 485e-10*self.MP/self.M1
        widthm  = 1.7e-10*self.MP/self.M1

        #define the different RHSs for each equation
        rhs1 =      -d*(N1-n1eq)

        rhs2 = epstt*d*(N1-n1eq)-0.5*w1*(2*c1t*c1tc*Ntt + c1m*c1tc*Ntm + c1e*c1tc*Nte + np.conjugate(c1m*c1tc*Ntm+c1e*c1tc*Nte)                  )
        rhs3 = epsmm*d*(N1-n1eq)-0.5*w1*(2*c1m*c1mc*Nmm + c1m*c1tc*Ntm + c1e*c1mc*Nme + np.conjugate(c1m*c1tc*Ntm+c1e*c1mc*Nme)                  )
        rhs4 = epsee*d*(N1-n1eq)-0.5*w1*(2*c1e*c1ec*Nee + c1e*c1mc*Nme + c1e*c1tc*Nte + np.conjugate(c1e*c1mc*Nme+c1e*c1tc*Nte)                  )
        rhs5 = epstm*d*(N1-n1eq)-0.5*w1*(  c1t*c1mc*Nmm + c1e*c1mc*Nte + c1m*c1mc*Ntm + c1mc*c1t*Ntt + c1t*c1tc*Ntm + c1t*c1ec*np.conjugate(Nme) )
        rhs6 = epste*d*(N1-n1eq)-0.5*w1*(  c1t*c1ec*Nee + c1e*c1ec*Nte + c1m*c1ec*Ntm + c1t*c1ec*Ntt + c1t*c1mc*Nme + c1t*c1tc*Nte               )
        rhs7 = epsme*d*(N1-n1eq)-0.5*w1*(  c1m*c1ec*Nee + c1e*c1ec*Nme + c1m*c1ec*Nmm + c1t*c1ec*np.conjugate(Ntm)  + c1m*c1mc*Nme + c1m*c1tc*Nte)

        return [rhs1, rhs2, rhs3, rhs4, rhs5, rhs6, rhs7]

    def RHS_1DS_DM_ZeroWidthm(self, y0,z,epstt,epsmm,epsee,epstm,epste,epsme,c1t,c1m,c1e,k):
        N1      = y0[0]
        Ntt     = y0[1]
        Nmm     = y0[2]
        Nee     = y0[3]
        Ntm     = y0[4]
        Nte     = y0[5]
        Nme     = y0[6]

        d       = np.real(self.D1(k,z))
        w1      = np.real(self.W1(k,z))
        n1eq    = self.N1Eq(z)

        c1tc    = np.conjugate(c1t)
        c1mc    = np.conjugate(c1m)
        c1ec    = np.conjugate(c1e)

        widtht  = 485e-10*self.MP/self.M1
        widthm  = 1.7e-10*self.MP/self.M1

        #define the different RHSs for each equation
        rhs1 =      -d*(N1-n1eq)

        rhs2 = epstt*d*(N1-n1eq)-0.5*w1*(2*c1t*c1tc*Ntt + c1m*c1tc*Ntm + c1e*c1tc*Nte + np.conjugate(c1m*c1tc*Ntm+c1e*c1tc*Nte)                  )
        rhs3 = epsmm*d*(N1-n1eq)-0.5*w1*(2*c1m*c1mc*Nmm + c1m*c1tc*Ntm + c1e*c1mc*Nme + np.conjugate(c1m*c1tc*Ntm+c1e*c1mc*Nme)                  )
        rhs4 = epsee*d*(N1-n1eq)-0.5*w1*(2*c1e*c1ec*Nee + c1e*c1mc*Nme + c1e*c1tc*Nte + np.conjugate(c1e*c1mc*Nme+c1e*c1tc*Nte)                  )
        rhs5 = epstm*d*(N1-n1eq)-0.5*w1*(  c1t*c1mc*Nmm + c1e*c1mc*Nte + c1m*c1mc*Ntm + c1mc*c1t*Ntt + c1t*c1tc*Ntm + c1t*c1ec*np.conjugate(Nme) ) - widtht*Ntm
        rhs6 = epste*d*(N1-n1eq)-0.5*w1*(  c1t*c1ec*Nee + c1e*c1ec*Nte + c1m*c1ec*Ntm + c1t*c1ec*Ntt + c1t*c1mc*Nme + c1t*c1tc*Nte               ) - widtht*Nte
        rhs7 = epsme*d*(N1-n1eq)-0.5*w1*(  c1m*c1ec*Nee + c1e*c1ec*Nme + c1m*c1ec*Nmm + c1t*c1ec*np.conjugate(Ntm)  + c1m*c1mc*Nme + c1m*c1tc*Nte)

        return [rhs1, rhs2, rhs3, rhs4, rhs5, rhs6, rhs7]

    def RHS_1DS_Approx(self, y0,z,epstt,epsmm,epsee,c1t,c1m,c1e,k):
        N1      = y0[0]
        Ntt     = y0[1]
        Nmm     = y0[2]
        Nee     = y0[3]

        d 	= np.real(self.D1(k,z))
        w1 	= np.real(self.W1(k,z))
        n1eq 	= self.N1Eq(z)

        c1tc    = np.conjugate(c1t)
        c1mc    = np.conjugate(c1m)
        c1ec    = np.conjugate(c1e)


        #define the different RHSs for each equation
        rhs1 =      -d*(N1-n1eq)

        rhs2 = epstt*d*(N1-n1eq)-0.5*w1*(2*c1t*c1tc*Ntt)
        rhs3 = epsmm*d*(N1-n1eq)-0.5*w1*(2*c1m*c1mc*Nmm)
        rhs4 = epsee*d*(N1-n1eq)-0.5*w1*(2*c1e*c1ec*Nee)

        return [rhs1, rhs2, rhs3, rhs4]

    def RHS_2DS_Approx(self, y0, z, ETA, C, K):
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

    def RHS_2DS_DM(self, y0, zzz, ETA, C, K, W):
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

    def RHS_3DS_DM(self, y0, zzz, ETA, C, K, W):
        N1, N2, N3, Ntt, Nmm, Nee, Ntm, Nte, Nme = y0
        (eps1tt,eps1mm,eps1ee,eps1tm,eps1te,eps1me,eps2tt,eps2mm,eps2ee,
eps2tm,eps2te,eps2me,eps3tt,eps3mm,eps3ee,eps3tm,eps3te,eps3me) = ETA
        c1t,c1m,c1e,c2t,c2m,c2e,c3t,c3m,c3e = C
        k1term,k2term,k3term = K
        widtht,widthm = W
        d1      = np.real(self.D1(k1term, zzz))
        w1      = np.real(self.W1(k1term, zzz))
        d2      = np.real(self.D2(k2term, zzz))
        w2      = np.real(self.W2(k2term, zzz))
        d3      = np.real(self.D3(k3term, zzz))
        w3      = np.real(self.W3(k3term, zzz))
        n1eq    = self.N1Eq(zzz)
        n2eq    = self.N2Eq(zzz)
        n3eq    = self.N3Eq(zzz)

        c1tc    = np.conjugate(c1t)
        c1mc    = np.conjugate(c1m)
        c1ec    = np.conjugate(c1e)

        c2tc    = np.conjugate(c2t)
        c2mc    = np.conjugate(c2m)
        c2ec    = np.conjugate(c2e)

        c3tc    = np.conjugate(c3t)
        c3mc    = np.conjugate(c3m)
        c3ec    = np.conjugate(c3e)

        #define the different RHSs for each equation
        rhs1    =      - d1 * (N1-n1eq)

        rhs2    =      - d2 * (N2-n2eq) 

        rhs3    =      - d3 * (N3-n3eq)       

        rhs4    = (eps1tt * d1 * (N1-n1eq) + eps2tt * d2 * (N2-n2eq) + eps3tt * d3 * (N3-n3eq)
		  - 0.5 * w1 * (2 * c1t * c1tc * Ntt + c1m * c1tc * Ntm + c1e * c1tc * Nte + np.conjugate(c1m * c1tc * Ntm + c1e * c1tc * Nte))
		  - 0.5 * w2 * (2 * c2t * c2tc * Ntt + c2m * c2tc * Ntm + c2e * c2tc * Nte + np.conjugate(c2m * c2tc * Ntm + c2e * c2tc * Nte))
		  - 0.5 * w3 * (2 * c3t * c3tc * Ntt + c3m * c3tc * Ntm + c3e * c3tc * Nte + np.conjugate(c3m * c3tc * Ntm + c3e * c3tc * Nte)))
    
        rhs5    = (eps1mm * d1 * (N1-n1eq) + eps2mm * d2 * (N2-n2eq) + eps3mm * d3 * (N3-n3eq)
		  - 0.5 * w1 * (2 * c1m * c1mc * Nmm + c1m * c1tc * Ntm + c1e * c1mc * Nme + np.conjugate(c1m * c1tc * Ntm + c1e * c1mc * Nme))
		  - 0.5 * w2 * (2 * c2m * c2mc * Nmm + c2m * c2tc * Ntm + c2e * c2mc * Nme + np.conjugate(c2m * c2tc * Ntm + c2e * c2mc * Nme))
		  - 0.5 * w3 * (2 * c3m * c3mc * Nmm + c3m * c3tc * Ntm + c3e * c3mc * Nme + np.conjugate(c3m * c3tc * Ntm + c3e * c3mc * Nme)))

        rhs6    = (eps1ee * d1 * (N1-n1eq) + eps2ee * d2 * (N2-n2eq) + eps3ee * d3 * (N3-n3eq)
		  - 0.5 * w1 * (2 * c1e * c1ec * Nee + c1e * c1mc * Nme + c1e * c1tc * Nte + np.conjugate(c1e * c1mc * Nme + c1e * c1tc * Nte))
		  - 0.5 * w2 * (2 * c2e * c2ec * Nee + c2e * c2mc * Nme + c2e * c2tc * Nte + np.conjugate(c2e * c2mc * Nme + c2e * c2tc * Nte))
		  - 0.5 * w3 * (2 * c3e * c3ec * Nee + c3e * c3mc * Nme + c3e * c3tc * Nte + np.conjugate(c3e * c3mc * Nme + c3e * c3tc * Nte)))
    
        rhs7    = (eps1tm * d1 * (N1-n1eq) + eps2tm * d2 * (N2-n2eq) + eps3tm * d3 * (N3-n3eq)
		  - 0.5 * ((w1 * c1t * c1mc + w2 * c2t * c2mc + w3 * c3t * c3mc) * Nmm
		  +        (w1 * c1e * c1mc + w2 * c2e * c2mc + w3 * c3e * c3mc) * Nte
		  +        (w1 * c1m * c1mc + w2 * c2m * c2mc + w3 * c3m * c3mc) * Ntm
		  +        (w1 * c1mc * c1t + w2 * c2mc * c2t + w3 * c3mc * c3t) * Ntt
		  +        (w1 * c1t * c1tc + w2 * c2t * c2tc + w3 * c3t * c3tc + 2 * widtht + 2 * widthm) * Ntm
		  +        (w1 * c1t * c1ec + w2 * c2t * c2ec + w3 * c3t * c3ec) * np.conjugate(Nme)))

        rhs8    = (eps1te * d1 * (N1-n1eq) + eps2te * d2 * (N2-n2eq) + eps3te * d3 * (N3-n3eq)
		  - 0.5 * ((w1 * c1t * c1ec + w2 * c2t * c2ec + w3 * c3t * c3ec) * Nee
		  +        (w1 * c1e * c1ec + w2 * c2e * c2ec + w3 * c3e * c3ec) * Nte
		  +        (w1 * c1m * c1ec + w2 * c2m * c2ec + w3 * c3m * c3ec) * Ntm
		  +        (w1 * c1t * c1ec + w2 * c2t * c2ec + w3 * c3t * c3ec) * Ntt
		  +        (w1 * c1t * c1mc + w2 * c2t * c2mc + w3 * c3t * c3mc) * Nme
		  +        (w1 * c1t * c1tc + w2 * c2t * c2tc + w3 * c3t * c3tc + 2 * widtht) * Nte))
    
        rhs9    = (eps1me * d1 * (N1-n1eq) + eps2me * d2 * (N2-n2eq) + eps3me * d3 * (N3-n3eq)
		  - 0.5 * ((w1 * c1m * c1ec + w2 * c2m * c2ec + w3 * c3m * c3ec) * Nee
		  +        (w1 * c1e * c1ec + w2 * c2e * c2ec + w3 * c3e * c3ec + 2 * widthm) * Nme
		  +        (w1 * c1m * c1ec + w2 * c2m * c2ec + w3 * c3m * c3ec) * Nmm
		  +        (w1 * c1t * c1ec + w2 * c2t * c2ec + w3 * c3t * c3ec) * np.conjugate(Ntm)
		  +        (w1 * c1m * c1mc + w2 * c2m * c2mc + w3 * c3m * c3mc) * Nme
		  +        (w1 * c1m * c1tc + w2 * c2m * c2tc + w3 * c3m * c3tc) * Nte))
    
        RHStemp = [rhs1, rhs2, rhs3, rhs4, rhs5, rhs6, rhs7, rhs8, rhs9]
        return RHStemp

    def RHS_3DS_DM_scattering(self, y0, zzz, ETA, C, K, W):
        N1, N2, N3, Ntt, Nmm, Nee, Ntm, Nte, Nme = y0
        (eps1tt,eps1mm,eps1ee,eps1tm,eps1te,eps1me,eps2tt,eps2mm,eps2ee,
eps2tm,eps2te,eps2me,eps3tt,eps3mm,eps3ee,eps3tm,eps3te,eps3me) = ETA
        c1t,c1m,c1e,c2t,c2m,c2e,c3t,c3m,c3e = C
        k1term,k2term,k3term = K
        widtht,widthm = W
        d1      = np.real(self.DS(k1term, zzz))
        w1      = self.j(zzz)*np.real(self.W1(k1term, zzz))
        d2      = np.real(self.D2(k2term, zzz))
        w2      = np.real(self.W2(k2term, zzz))
        d3      = np.real(self.D3(k3term, zzz))
        w3      = np.real(self.W3(k3term, zzz))
        n1eq    = self.N1Eq(zzz)
        n2eq    = self.N2Eq(zzz)
        n3eq    = self.N3Eq(zzz)

        c1tc    = np.conjugate(c1t)
        c1mc    = np.conjugate(c1m)
        c1ec    = np.conjugate(c1e)

        c2tc    = np.conjugate(c2t)
        c2mc    = np.conjugate(c2m)
        c2ec    = np.conjugate(c2e)

        c3tc    = np.conjugate(c3t)
        c3mc    = np.conjugate(c3m)
        c3ec    = np.conjugate(c3e)

        #define the different RHSs for each equation
        rhs1    =      - d1 * (N1-n1eq)

        rhs2    =      - d2 * (N2-n2eq) 

        rhs3    =      - d3 * (N3-n3eq)       

        rhs4    = (eps1tt * d1 * (N1-n1eq) + eps2tt * d2 * (N2-n2eq) + eps3tt * d3 * (N3-n3eq)
		  - 0.5 * w1 * (2 * c1t * c1tc * Ntt + c1m * c1tc * Ntm + c1e * c1tc * Nte + np.conjugate(c1m * c1tc * Ntm + c1e * c1tc * Nte))
		  - 0.5 * w2 * (2 * c2t * c2tc * Ntt + c2m * c2tc * Ntm + c2e * c2tc * Nte + np.conjugate(c2m * c2tc * Ntm + c2e * c2tc * Nte))
		  - 0.5 * w3 * (2 * c3t * c3tc * Ntt + c3m * c3tc * Ntm + c3e * c3tc * Nte + np.conjugate(c3m * c3tc * Ntm + c3e * c3tc * Nte)))
    
        rhs5    = (eps1mm * d1 * (N1-n1eq) + eps2mm * d2 * (N2-n2eq) + eps3mm * d3 * (N3-n3eq)
		  - 0.5 * w1 * (2 * c1m * c1mc * Nmm + c1m * c1tc * Ntm + c1e * c1mc * Nme + np.conjugate(c1m * c1tc * Ntm + c1e * c1mc * Nme))
		  - 0.5 * w2 * (2 * c2m * c2mc * Nmm + c2m * c2tc * Ntm + c2e * c2mc * Nme + np.conjugate(c2m * c2tc * Ntm + c2e * c2mc * Nme))
		  - 0.5 * w3 * (2 * c3m * c3mc * Nmm + c3m * c3tc * Ntm + c3e * c3mc * Nme + np.conjugate(c3m * c3tc * Ntm + c3e * c3mc * Nme)))

        rhs6    = (eps1ee * d1 * (N1-n1eq) + eps2ee * d2 * (N2-n2eq) + eps3ee * d3 * (N3-n3eq)
		  - 0.5 * w1 * (2 * c1e * c1ec * Nee + c1e * c1mc * Nme + c1e * c1tc * Nte + np.conjugate(c1e * c1mc * Nme + c1e * c1tc * Nte))
		  - 0.5 * w2 * (2 * c2e * c2ec * Nee + c2e * c2mc * Nme + c2e * c2tc * Nte + np.conjugate(c2e * c2mc * Nme + c2e * c2tc * Nte))
		  - 0.5 * w3 * (2 * c3e * c3ec * Nee + c3e * c3mc * Nme + c3e * c3tc * Nte + np.conjugate(c3e * c3mc * Nme + c3e * c3tc * Nte)))
    
        rhs7    = (eps1tm * d1 * (N1-n1eq) + eps2tm * d2 * (N2-n2eq) + eps3tm * d3 * (N3-n3eq)
		  - 0.5 * ((w1 * c1t * c1mc + w2 * c2t * c2mc + w3 * c3t * c3mc) * Nmm
		  +        (w1 * c1e * c1mc + w2 * c2e * c2mc + w3 * c3e * c3mc) * Nte
		  +        (w1 * c1m * c1mc + w2 * c2m * c2mc + w3 * c3m * c3mc) * Ntm
		  +        (w1 * c1mc * c1t + w2 * c2mc * c2t + w3 * c3mc * c3t) * Ntt
		  +        (w1 * c1t * c1tc + w2 * c2t * c2tc + w3 * c3t * c3tc + 2 * widtht + 2 * widthm) * Ntm
		  +        (w1 * c1t * c1ec + w2 * c2t * c2ec + w3 * c3t * c3ec) * np.conjugate(Nme)))

        rhs8    = (eps1te * d1 * (N1-n1eq) + eps2te * d2 * (N2-n2eq) + eps3te * d3 * (N3-n3eq)
		  - 0.5 * ((w1 * c1t * c1ec + w2 * c2t * c2ec + w3 * c3t * c3ec) * Nee
		  +        (w1 * c1e * c1ec + w2 * c2e * c2ec + w3 * c3e * c3ec) * Nte
		  +        (w1 * c1m * c1ec + w2 * c2m * c2ec + w3 * c3m * c3ec) * Ntm
		  +        (w1 * c1t * c1ec + w2 * c2t * c2ec + w3 * c3t * c3ec) * Ntt
		  +        (w1 * c1t * c1mc + w2 * c2t * c2mc + w3 * c3t * c3mc) * Nme
		  +        (w1 * c1t * c1tc + w2 * c2t * c2tc + w3 * c3t * c3tc + 2 * widtht) * Nte))
    
        rhs9    = (eps1me * d1 * (N1-n1eq) + eps2me * d2 * (N2-n2eq) + eps3me * d3 * (N3-n3eq)
		  - 0.5 * ((w1 * c1m * c1ec + w2 * c2m * c2ec + w3 * c3m * c3ec) * Nee
		  +        (w1 * c1e * c1ec + w2 * c2e * c2ec + w3 * c3e * c3ec + 2 * widthm) * Nme
		  +        (w1 * c1m * c1ec + w2 * c2m * c2ec + w3 * c3m * c3ec) * Nmm
		  +        (w1 * c1t * c1ec + w2 * c2t * c2ec + w3 * c3t * c3ec) * np.conjugate(Ntm)
		  +        (w1 * c1m * c1mc + w2 * c2m * c2mc + w3 * c3m * c3mc) * Nme
		  +        (w1 * c1m * c1tc + w2 * c2m * c2tc + w3 * c3m * c3tc) * Nte))
    
        RHStemp = [rhs1, rhs2, rhs3, rhs4, rhs5, rhs6, rhs7, rhs8, rhs9]
        return RHStemp

    def RHS_3DS_DM_OOEtauR(self, y0, zzz, ETA, C, K, W):
        N1, N2, N3, Ntt, Nmm, Nee, Ntm, Nte, Nme, Ntr = y0
        (eps1tt,eps1mm,eps1ee,eps1tm,eps1te,eps1me,eps2tt,eps2mm,eps2ee,
eps2tm,eps2te,eps2me,eps3tt,eps3mm,eps3ee,eps3tm,eps3te,eps3me) = ETA
        c1t,c1m,c1e,c2t,c2m,c2e,c3t,c3m,c3e = C
        k1term,k2term,k3term = K
        widtht,widthm = W
        d1      = np.real(self.D1(k1term, zzz))
        w1      = np.real(self.W1(k1term, zzz))
        d2      = np.real(self.D2(k2term, zzz))
        w2      = np.real(self.W2(k2term, zzz))
        d3      = np.real(self.D3(k3term, zzz))
        w3      = np.real(self.W3(k3term, zzz))
        n1eq    = self.N1Eq(zzz)
        n2eq    = self.N2Eq(zzz)
        n3eq    = self.N3Eq(zzz)

        n1eq    = self.N1Eq(zzz)
        n2eq    = self.N2Eq(zzz)
        n3eq    = self.N3Eq(zzz)

        c1tc    = np.conjugate(c1t)
        c1mc    = np.conjugate(c1m)
        c1ec    = np.conjugate(c1e)

        c2tc    = np.conjugate(c2t)
        c2mc    = np.conjugate(c2m)
        c2ec    = np.conjugate(c2e)

        c3tc    = np.conjugate(c3t)
        c3mc    = np.conjugate(c3m)
        c3ec    = np.conjugate(c3e)

        #define the different RHSs for each equation
        rhs1    =      - d1 * (N1-n1eq)

        rhs2    =      - d2 * (N2-n2eq) 

        rhs3    =      - d3 * (N3-n3eq)       

        rhs4    = (eps1tt * d1 * (N1-n1eq) + eps2tt * d2 * (N2-n2eq) + eps3tt * d3 * (N3-n3eq)
		  - 0.5 * w1 * (2 * c1t * c1tc * Ntt + c1m * c1tc * Ntm + c1e * c1tc * Nte + np.conjugate(c1m * c1tc * Ntm + c1e * c1tc * Nte))
		  - 0.5 * w2 * (2 * c2t * c2tc * Ntt + c2m * c2tc * Ntm + c2e * c2tc * Nte + np.conjugate(c2m * c2tc * Ntm + c2e * c2tc * Nte))
		  - 0.5 * w3 * (2 * c3t * c3tc * Ntt + c3m * c3tc * Ntm + c3e * c3tc * Nte + np.conjugate(c3m * c3tc * Ntm + c3e * c3tc * Nte))
                  - 2 * widtht * Ntt + 4 * widtht * Ntr)
    
        rhs5    = (eps1mm * d1 * (N1-n1eq) + eps2mm * d2 * (N2-n2eq) + eps3mm * d3 * (N3-n3eq)
		  - 0.5 * w1 * (2 * c1m * c1mc * Nmm + c1m * c1tc * Ntm + c1e * c1mc * Nme + np.conjugate(c1m * c1tc * Ntm + c1e * c1mc * Nme))
		  - 0.5 * w2 * (2 * c2m * c2mc * Nmm + c2m * c2tc * Ntm + c2e * c2mc * Nme + np.conjugate(c2m * c2tc * Ntm + c2e * c2mc * Nme))
		  - 0.5 * w3 * (2 * c3m * c3mc * Nmm + c3m * c3tc * Ntm + c3e * c3mc * Nme + np.conjugate(c3m * c3tc * Ntm + c3e * c3mc * Nme)))

        rhs6    = (eps1ee * d1 * (N1-n1eq) + eps2ee * d2 * (N2-n2eq) + eps3ee * d3 * (N3-n3eq)
		  - 0.5 * w1 * (2 * c1e * c1ec * Nee + c1e * c1mc * Nme + c1e * c1tc * Nte + np.conjugate(c1e * c1mc * Nme + c1e * c1tc * Nte))
		  - 0.5 * w2 * (2 * c2e * c2ec * Nee + c2e * c2mc * Nme + c2e * c2tc * Nte + np.conjugate(c2e * c2mc * Nme + c2e * c2tc * Nte))
		  - 0.5 * w3 * (2 * c3e * c3ec * Nee + c3e * c3mc * Nme + c3e * c3tc * Nte + np.conjugate(c3e * c3mc * Nme + c3e * c3tc * Nte)))
    
        rhs7    = (eps1tm * d1 * (N1-n1eq) + eps2tm * d2 * (N2-n2eq) + eps3tm * d3 * (N3-n3eq)
		  - 0.5 * ((w1 * c1t * c1mc + w2 * c2t * c2mc + w3 * c3t * c3mc) * Nmm
		  +        (w1 * c1e * c1mc + w2 * c2e * c2mc + w3 * c3e * c3mc) * Nte
		  +        (w1 * c1m * c1mc + w2 * c2m * c2mc + w3 * c3m * c3mc) * Ntm
		  +        (w1 * c1mc * c1t + w2 * c2mc * c2t + w3 * c3mc * c3t) * Ntt
		  +        (w1 * c1t * c1tc + w2 * c2t * c2tc + w3 * c3t * c3tc + 2 * widtht + 2 * widthm) * Ntm
		  +        (w1 * c1t * c1ec + w2 * c2t * c2ec + w3 * c3t * c3ec) * np.conjugate(Nme)))

        rhs8    = (eps1te * d1 * (N1-n1eq) + eps2te * d2 * (N2-n2eq) + eps3te * d3 * (N3-n3eq)
		  - 0.5 * ((w1 * c1t * c1ec + w2 * c2t * c2ec + w3 * c3t * c3ec) * Nee
		  +        (w1 * c1e * c1ec + w2 * c2e * c2ec + w3 * c3e * c3ec) * Nte
		  +        (w1 * c1m * c1ec + w2 * c2m * c2ec + w3 * c3m * c3ec) * Ntm
		  +        (w1 * c1t * c1ec + w2 * c2t * c2ec + w3 * c3t * c3ec) * Ntt
		  +        (w1 * c1t * c1mc + w2 * c2t * c2mc + w3 * c3t * c3mc) * Nme
		  +        (w1 * c1t * c1tc + w2 * c2t * c2tc + w3 * c3t * c3tc + 2 * widtht) * Nte))
    
        rhs9    = (eps1me * d1 * (N1-n1eq) + eps2me * d2 * (N2-n2eq) + eps3me * d3 * (N3-n3eq)
		  - 0.5 * ((w1 * c1m * c1ec + w2 * c2m * c2ec + w3 * c3m * c3ec) * Nee
		  +        (w1 * c1e * c1ec + w2 * c2e * c2ec + w3 * c3e * c3ec + 2 * widthm) * Nme
		  +        (w1 * c1m * c1ec + w2 * c2m * c2ec + w3 * c3m * c3ec) * Nmm
		  +        (w1 * c1t * c1ec + w2 * c2t * c2ec + w3 * c3t * c3ec) * np.conjugate(Ntm)
		  +        (w1 * c1m * c1mc + w2 * c2m * c2mc + w3 * c3m * c3mc) * Nme
		  +        (w1 * c1m * c1tc + w2 * c2m * c2tc + w3 * c3m * c3tc) * Nte))
    
        rhs10   = 2 * widtht * (Ntt - 2 * Ntr)

        RHStemp = [rhs1, rhs2, rhs3, rhs4, rhs5, rhs6, rhs7, rhs8, rhs9, rhs10]
        return RHStemp

    def RHS_3DS_DM_scattering_OOEtauR(self, y0, zzz, ETA, C, K, W):
        N1, N2, N3, Ntt, Nmm, Nee, Ntm, Nte, Nme, Ntr = y0
        (eps1tt,eps1mm,eps1ee,eps1tm,eps1te,eps1me,eps2tt,eps2mm,eps2ee,
eps2tm,eps2te,eps2me,eps3tt,eps3mm,eps3ee,eps3tm,eps3te,eps3me) = ETA
        c1t,c1m,c1e,c2t,c2m,c2e,c3t,c3m,c3e = C
        k1term,k2term,k3term = K
        widtht,widthm = W
        d1      = np.real(self.DS(k1term, zzz))
        w1      = self.j(zzz)*np.real(self.W1(k1term, zzz))
        d2      = np.real(self.D2(k2term, zzz))
        w2      = np.real(self.W2(k2term, zzz))
        d3      = np.real(self.D3(k3term, zzz))
        w3      = np.real(self.W3(k3term, zzz))
        n1eq    = self.N1Eq(zzz)
        n2eq    = self.N2Eq(zzz)
        n3eq    = self.N3Eq(zzz)

        c1tc    = np.conjugate(c1t)
        c1mc    = np.conjugate(c1m)
        c1ec    = np.conjugate(c1e)

        c2tc    = np.conjugate(c2t)
        c2mc    = np.conjugate(c2m)
        c2ec    = np.conjugate(c2e)

        c3tc    = np.conjugate(c3t)
        c3mc    = np.conjugate(c3m)
        c3ec    = np.conjugate(c3e)

        #define the different RHSs for each equation
        rhs1    =      - d1 * (N1-n1eq)

        rhs2    =      - d2 * (N2-n2eq) 

        rhs3    =      - d3 * (N3-n3eq)       

        rhs4    = (eps1tt * d1 * (N1-n1eq) + eps2tt * d2 * (N2-n2eq) + eps3tt * d3 * (N3-n3eq)
		  - 0.5 * w1 * (2 * c1t * c1tc * Ntt + c1m * c1tc * Ntm + c1e * c1tc * Nte + np.conjugate(c1m * c1tc * Ntm + c1e * c1tc * Nte))
		  - 0.5 * w2 * (2 * c2t * c2tc * Ntt + c2m * c2tc * Ntm + c2e * c2tc * Nte + np.conjugate(c2m * c2tc * Ntm + c2e * c2tc * Nte))
		  - 0.5 * w3 * (2 * c3t * c3tc * Ntt + c3m * c3tc * Ntm + c3e * c3tc * Nte + np.conjugate(c3m * c3tc * Ntm + c3e * c3tc * Nte))
                  - 2 * widtht * Ntt + 4 * widtht * Ntr)
    
        rhs5    = (eps1mm * d1 * (N1-n1eq) + eps2mm * d2 * (N2-n2eq) + eps3mm * d3 * (N3-n3eq)
		  - 0.5 * w1 * (2 * c1m * c1mc * Nmm + c1m * c1tc * Ntm + c1e * c1mc * Nme + np.conjugate(c1m * c1tc * Ntm + c1e * c1mc * Nme))
		  - 0.5 * w2 * (2 * c2m * c2mc * Nmm + c2m * c2tc * Ntm + c2e * c2mc * Nme + np.conjugate(c2m * c2tc * Ntm + c2e * c2mc * Nme))
		  - 0.5 * w3 * (2 * c3m * c3mc * Nmm + c3m * c3tc * Ntm + c3e * c3mc * Nme + np.conjugate(c3m * c3tc * Ntm + c3e * c3mc * Nme)))

        rhs6    = (eps1ee * d1 * (N1-n1eq) + eps2ee * d2 * (N2-n2eq) + eps3ee * d3 * (N3-n3eq)
		  - 0.5 * w1 * (2 * c1e * c1ec * Nee + c1e * c1mc * Nme + c1e * c1tc * Nte + np.conjugate(c1e * c1mc * Nme + c1e * c1tc * Nte))
		  - 0.5 * w2 * (2 * c2e * c2ec * Nee + c2e * c2mc * Nme + c2e * c2tc * Nte + np.conjugate(c2e * c2mc * Nme + c2e * c2tc * Nte))
		  - 0.5 * w3 * (2 * c3e * c3ec * Nee + c3e * c3mc * Nme + c3e * c3tc * Nte + np.conjugate(c3e * c3mc * Nme + c3e * c3tc * Nte)))
    
        rhs7    = (eps1tm * d1 * (N1-n1eq) + eps2tm * d2 * (N2-n2eq) + eps3tm * d3 * (N3-n3eq)
		  - 0.5 * ((w1 * c1t * c1mc + w2 * c2t * c2mc + w3 * c3t * c3mc) * Nmm
		  +        (w1 * c1e * c1mc + w2 * c2e * c2mc + w3 * c3e * c3mc) * Nte
		  +        (w1 * c1m * c1mc + w2 * c2m * c2mc + w3 * c3m * c3mc) * Ntm
		  +        (w1 * c1mc * c1t + w2 * c2mc * c2t + w3 * c3mc * c3t) * Ntt
		  +        (w1 * c1t * c1tc + w2 * c2t * c2tc + w3 * c3t * c3tc + 2 * widtht + 2 * widthm) * Ntm
		  +        (w1 * c1t * c1ec + w2 * c2t * c2ec + w3 * c3t * c3ec) * np.conjugate(Nme)))

        rhs8    = (eps1te * d1 * (N1-n1eq) + eps2te * d2 * (N2-n2eq) + eps3te * d3 * (N3-n3eq)
		  - 0.5 * ((w1 * c1t * c1ec + w2 * c2t * c2ec + w3 * c3t * c3ec) * Nee
		  +        (w1 * c1e * c1ec + w2 * c2e * c2ec + w3 * c3e * c3ec) * Nte
		  +        (w1 * c1m * c1ec + w2 * c2m * c2ec + w3 * c3m * c3ec) * Ntm
		  +        (w1 * c1t * c1ec + w2 * c2t * c2ec + w3 * c3t * c3ec) * Ntt
		  +        (w1 * c1t * c1mc + w2 * c2t * c2mc + w3 * c3t * c3mc) * Nme
		  +        (w1 * c1t * c1tc + w2 * c2t * c2tc + w3 * c3t * c3tc + 2 * widtht) * Nte))
    
        rhs9    = (eps1me * d1 * (N1-n1eq) + eps2me * d2 * (N2-n2eq) + eps3me * d3 * (N3-n3eq)
		  - 0.5 * ((w1 * c1m * c1ec + w2 * c2m * c2ec + w3 * c3m * c3ec) * Nee
		  +        (w1 * c1e * c1ec + w2 * c2e * c2ec + w3 * c3e * c3ec + 2 * widthm) * Nme
		  +        (w1 * c1m * c1ec + w2 * c2m * c2ec + w3 * c3m * c3ec) * Nmm
		  +        (w1 * c1t * c1ec + w2 * c2t * c2ec + w3 * c3t * c3ec) * np.conjugate(Ntm)
		  +        (w1 * c1m * c1mc + w2 * c2m * c2mc + w3 * c3m * c3mc) * Nme
		  +        (w1 * c1m * c1tc + w2 * c2m * c2tc + w3 * c3m * c3tc) * Nte))
    
        rhs10   = 2 * widtht * (Ntt - 2 * Ntr)

        RHStemp = [rhs1, rhs2, rhs3, rhs4, rhs5, rhs6, rhs7, rhs8, rhs9, rhs10]
        return RHStemp

    def RHS_1DS_DM_scattering_OOEtauR(self, y0, zzz, ETA, C, K, W):
        N1, Ntt, Nmm, Nee, Ntm, Nte, Nme, Ntr = y0
        (eps1tt,eps1mm,eps1ee,eps1tm,eps1te,eps1me,eps2tt,eps2mm,eps2ee,
eps2tm,eps2te,eps2me,eps3tt,eps3mm,eps3ee,eps3tm,eps3te,eps3me) = ETA
        c1t,c1m,c1e,c2t,c2m,c2e,c3t,c3m,c3e = C
        k1term,k2term,k3term = K
        widtht,widthm = W
        d1      = np.real(self.DS(k1term, zzz))
        w1      = self.j(zzz)*np.real(self.W1(k1term, zzz))

        n1eq    = self.N1Eq(zzz)

        c1tc    = np.conjugate(c1t)
        c1mc    = np.conjugate(c1m)
        c1ec    = np.conjugate(c1e)


        #define the different RHSs for each equation
        rhs1    =      - d1 * (N1-n1eq) 

        rhs2    = (eps1tt * d1 * (N1-n1eq)
		  - 0.5 * w1 * (2 * c1t * c1tc * Ntt + c1m * c1tc * Ntm + c1e * c1tc * Nte + np.conjugate(c1m * c1tc * Ntm + c1e * c1tc * Nte))
                  - 2 * widtht * Ntt + 4 * widtht * Ntr)
    
        rhs3    = (eps1mm * d1 * (N1-n1eq)
		  - 0.5 * w1 * (2 * c1m * c1mc * Nmm + c1m * c1tc * Ntm + c1e * c1mc * Nme + np.conjugate(c1m * c1tc * Ntm + c1e * c1mc * Nme)))

        rhs4    = (eps1ee * d1 * (N1-n1eq)
		  - 0.5 * w1 * (2 * c1e * c1ec * Nee + c1e * c1mc * Nme + c1e * c1tc * Nte + np.conjugate(c1e * c1mc * Nme + c1e * c1tc * Nte)))
    
        rhs5    = (eps1tm * d1 * (N1-n1eq)
		  - 0.5 * ((w1 * c1t * c1mc) * Nmm
		  +        (w1 * c1e * c1mc) * Nte
		  +        (w1 * c1m * c1mc) * Ntm
		  +        (w1 * c1mc * c1t) * Ntt
		  +        (w1 * c1t * c1tc + 2 * widtht + 2 * widthm) * Ntm
		  +        (w1 * c1t * c1ec) * np.conjugate(Nme)))

        rhs6    = (eps1te * d1 * (N1-n1eq)
		  - 0.5 * ((w1 * c1t * c1ec) * Nee
		  +        (w1 * c1e * c1ec) * Nte
		  +        (w1 * c1m * c1ec) * Ntm
		  +        (w1 * c1t * c1ec) * Ntt
		  +        (w1 * c1t * c1mc) * Nme
		  +        (w1 * c1t * c1tc + 2 * widtht) * Nte))
    
        rhs7    = (eps1me * d1 * (N1-n1eq)
		  - 0.5 * ((w1 * c1m * c1ec) * Nee
		  +        (w1 * c1e * c1ec + 2 * widthm) * Nme
		  +        (w1 * c1m * c1ec) * Nmm
		  +        (w1 * c1t * c1ec) * np.conjugate(Ntm)
		  +        (w1 * c1m * c1mc) * Nme
		  +        (w1 * c1m * c1tc) * Nte))
    
        rhs8   = 2 * widtht * (Ntt - 2 * Ntr)

        RHStemp = [rhs1, rhs2, rhs3, rhs4, rhs5, rhs6, rhs7, rhs8]
        return RHStemp

    def RHS_3DS_DM_Blanchett(self, y0, zzz, ETA, C, K, W):
        N1, N2, N3, Ntt, Nmm, Nee, Ntm, Nte, Nme, Ntr, ptm, pte, pme = y0
        (eps1tt,eps1mm,eps1ee,eps1tm,eps1te,eps1me,eps2tt,eps2mm,eps2ee,
eps2tm,eps2te,eps2me,eps3tt,eps3mm,eps3ee,eps3tm,eps3te,eps3me) = ETA
        c1t,c1m,c1e,c2t,c2m,c2e,c3t,c3m,c3e = C
        k1term,k2term,k3term = K
        widtht,widthm,widthRt = W

        d1      = np.real(self.DS(k1term, zzz))
        w1      = self.j(zzz)*np.real(self.W1(k1term, zzz))
        d2      = np.real(self.D2(k2term, zzz))
        w2      = np.real(self.W2(k2term, zzz))
        d3      = np.real(self.D3(k3term, zzz))
        w3      = np.real(self.W3(k3term, zzz))
        n1eq    = self.N1Eq(zzz)
        n2eq    = self.N2Eq(zzz)
        n3eq    = self.N3Eq(zzz)

        c1tc    = np.conjugate(c1t)
        c1mc    = np.conjugate(c1m)
        c1ec    = np.conjugate(c1e)

        c2tc    = np.conjugate(c2t)
        c2mc    = np.conjugate(c2m)
        c2ec    = np.conjugate(c2e)

        c3tc    = np.conjugate(c3t)
        c3mc    = np.conjugate(c3m)
        c3ec    = np.conjugate(c3e)

        #define the different RHSs for each equation
        rhs1    =      - d1 * (N1-n1eq)

        rhs2    =      - d2 * (N2-n2eq) 

        rhs3    =      - d3 * (N3-n3eq)       

        rhs4    = (eps1tt * d1 * (N1-n1eq) + eps2tt * d2 * (N2-n2eq) + eps3tt * d3 * (N3-n3eq)
		  - 0.5 * w1 * (2 * c1t * c1tc * Ntt + c1m * c1tc * Ntm + c1e * c1tc * Nte + np.conjugate(c1m * c1tc * Ntm + c1e * c1tc * Nte))
		  - 0.5 * w2 * (2 * c2t * c2tc * Ntt + c2m * c2tc * Ntm + c2e * c2tc * Nte + np.conjugate(c2m * c2tc * Ntm + c2e * c2tc * Nte))
		  - 0.5 * w3 * (2 * c3t * c3tc * Ntt + c3m * c3tc * Ntm + c3e * c3tc * Nte + np.conjugate(c3m * c3tc * Ntm + c3e * c3tc * Nte))
                  - 2 * widtht * Ntt + 4 * widtht * Ntr)
    
        rhs5    = (eps1mm * d1 * (N1-n1eq) + eps2mm * d2 * (N2-n2eq) + eps3mm * d3 * (N3-n3eq)
		  - 0.5 * w1 * (2 * c1m * c1mc * Nmm + c1m * c1tc * Ntm + c1e * c1mc * Nme + np.conjugate(c1m * c1tc * Ntm + c1e * c1mc * Nme))
		  - 0.5 * w2 * (2 * c2m * c2mc * Nmm + c2m * c2tc * Ntm + c2e * c2mc * Nme + np.conjugate(c2m * c2tc * Ntm + c2e * c2mc * Nme))
		  - 0.5 * w3 * (2 * c3m * c3mc * Nmm + c3m * c3tc * Ntm + c3e * c3mc * Nme + np.conjugate(c3m * c3tc * Ntm + c3e * c3mc * Nme)))

        rhs6    = (eps1ee * d1 * (N1-n1eq) + eps2ee * d2 * (N2-n2eq) + eps3ee * d3 * (N3-n3eq)
		  - 0.5 * w1 * (2 * c1e * c1ec * Nee + c1e * c1mc * Nme + c1e * c1tc * Nte + np.conjugate(c1e * c1mc * Nme + c1e * c1tc * Nte))
		  - 0.5 * w2 * (2 * c2e * c2ec * Nee + c2e * c2mc * Nme + c2e * c2tc * Nte + np.conjugate(c2e * c2mc * Nme + c2e * c2tc * Nte))
		  - 0.5 * w3 * (2 * c3e * c3ec * Nee + c3e * c3mc * Nme + c3e * c3tc * Nte + np.conjugate(c3e * c3mc * Nme + c3e * c3tc * Nte)))
    
        rhs7    = (eps1tm * d1 * (N1-n1eq) + eps2tm * d2 * (N2-n2eq) + eps3tm * d3 * (N3-n3eq)
		  - 0.5 * ((w1 * c1t * c1mc + w2 * c2t * c2mc + w3 * c3t * c3mc) * Nmm
		  +        (w1 * c1e * c1mc + w2 * c2e * c2mc + w3 * c3e * c3mc) * Nte
		  +        (w1 * c1m * c1mc + w2 * c2m * c2mc + w3 * c3m * c3mc) * Ntm
		  +        (w1 * c1mc * c1t + w2 * c2mc * c2t + w3 * c3mc * c3t) * Ntt
		  +        (w1 * c1t * c1tc + w2 * c2t * c2tc + w3 * c3t * c3tc + 2 * widtht + 2 * widthm) * Ntm + 1j*widthRt*ptm
		  +        (w1 * c1t * c1ec + w2 * c2t * c2ec + w3 * c3t * c3ec) * np.conjugate(Nme)))

        rhs8    = (eps1te * d1 * (N1-n1eq) + eps2te * d2 * (N2-n2eq) + eps3te * d3 * (N3-n3eq)
		  - 0.5 * ((w1 * c1t * c1ec + w2 * c2t * c2ec + w3 * c3t * c3ec) * Nee
		  +        (w1 * c1e * c1ec + w2 * c2e * c2ec + w3 * c3e * c3ec) * Nte
		  +        (w1 * c1m * c1ec + w2 * c2m * c2ec + w3 * c3m * c3ec) * Ntm
		  +        (w1 * c1t * c1ec + w2 * c2t * c2ec + w3 * c3t * c3ec) * Ntt
		  +        (w1 * c1t * c1mc + w2 * c2t * c2mc + w3 * c3t * c3mc) * Nme
		  +        (w1 * c1t * c1tc + w2 * c2t * c2tc + w3 * c3t * c3tc + 2 * widtht) * Nte)) + 1j*widthRt*pte
    
        rhs9    = (eps1me * d1 * (N1-n1eq) + eps2me * d2 * (N2-n2eq) + eps3me * d3 * (N3-n3eq)
		  - 0.5 * ((w1 * c1m * c1ec + w2 * c2m * c2ec + w3 * c3m * c3ec) * Nee
		  +        (w1 * c1e * c1ec + w2 * c2e * c2ec + w3 * c3e * c3ec + 2 * widthm) * Nme
		  +        (w1 * c1m * c1ec + w2 * c2m * c2ec + w3 * c3m * c3ec) * Nmm
		  +        (w1 * c1t * c1ec + w2 * c2t * c2ec + w3 * c3t * c3ec) * np.conjugate(Ntm)
		  +        (w1 * c1m * c1mc + w2 * c2m * c2mc + w3 * c3m * c3mc) * Nme
		  +        (w1 * c1m * c1tc + w2 * c2m * c2tc + w3 * c3m * c3tc) * Nte)) + 1j*widthRt*pme
    
        rhs10   = 2 * widtht * (Ntt - 2 * Ntr)

        rhs11   = -1j*widthRt*Ntm

        rhs12   = -1j*widthRt*Nte

        rhs13   = -1j*widthRt*Nme

        RHStemp = [rhs1, rhs2, rhs3, rhs4, rhs5, rhs6, rhs7, rhs8, rhs9, rhs10, rhs11, rhs12, rhs13]
        return RHStemp

    @property
    def getEtaB_3DS_DM(self):

        #Define fixed quantities for BEs
        _ETA = [
            np.real(self.epsilon1ab(2,2)),
            np.real(self.epsilon1ab(1,1)),
            np.real(self.epsilon1ab(0,0)),
                    self.epsilon1ab(2,1) ,
                    self.epsilon1ab(2,0) ,
                    self.epsilon1ab(1,0) ,
            np.real(self.epsilon2ab(2,2)),
            np.real(self.epsilon2ab(1,1)),
            np.real(self.epsilon2ab(0,0)),
                    self.epsilon2ab(2,1) ,
                    self.epsilon2ab(2,0) ,
                    self.epsilon2ab(1,0) ,
            np.real(self.epsilon3ab(2,2)),
            np.real(self.epsilon3ab(1,1)),
            np.real(self.epsilon3ab(0,0)),
                    self.epsilon3ab(2,1) ,
                    self.epsilon3ab(2,0) ,
                    self.epsilon3ab(1,0)]

        _C =   [self.c1a(2), self.c1a(1), self.c1a(0),
                self.c2a(2), self.c2a(1), self.c2a(0),
                self.c3a(2), self.c3a(1), self.c3a(0)]

        _K      = [np.real(self.k1), np.real(self.k2), np.real(self.k3)]
        _W      = [ 485e-10*self.MP/self.M1, 1.7e-10*self.MP/self.M1]

        y0      = np.array([0+0j,0+0j,0+0j,0+0j,0+0j,0+0j,0+0j,0+0j,0+0j], dtype=np.complex128)

        zcrit   = 1e100
        ys, _   = odeintw(self.RHS_3DS_DM, y0, self.xs, args = tuple([_ETA, _C , _K, _W]), full_output=1)
        nb      = np.real(self.sphalfact*(ys[-1,3]+ys[-1,4]+ys[-1,5]))

        return nb

    @property
    def getEtaB_3DS_DM_scattering(self):

        #Define fixed quantities for BEs
        _ETA = [
            np.real(self.epsilon1ab(2,2)),
            np.real(self.epsilon1ab(1,1)),
            np.real(self.epsilon1ab(0,0)),
                    self.epsilon1ab(2,1) ,
                    self.epsilon1ab(2,0) ,
                    self.epsilon1ab(1,0) ,
            np.real(self.epsilon2ab(2,2)),
            np.real(self.epsilon2ab(1,1)),
            np.real(self.epsilon2ab(0,0)),
                    self.epsilon2ab(2,1) ,
                    self.epsilon2ab(2,0) ,
                    self.epsilon2ab(1,0) ,
            np.real(self.epsilon3ab(2,2)),
            np.real(self.epsilon3ab(1,1)),
            np.real(self.epsilon3ab(0,0)),
                    self.epsilon3ab(2,1) ,
                    self.epsilon3ab(2,0) ,
                    self.epsilon3ab(1,0)]

        _C =   [self.c1a(2), self.c1a(1), self.c1a(0),
                self.c2a(2), self.c2a(1), self.c2a(0),
                self.c3a(2), self.c3a(1), self.c3a(0)]

        _K      = [np.real(self.k1), np.real(self.k2), np.real(self.k3)]
        _W      = [ 485e-10*self.MP/self.M1, 1.7e-10*self.MP/self.M1]

        y0      = np.array([0+0j,0+0j,0+0j,0+0j,0+0j,0+0j,0+0j,0+0j,0+0j], dtype=np.complex128)

        zcrit   = 1e100
        ys, _   = odeintw(self.RHS_3DS_DM_scattering, y0, self.xs, args = tuple([_ETA, _C , _K, _W]), full_output=1)
        nb      = np.real(self.sphalfact*(ys[-1,3]+ys[-1,4]+ys[-1,5]))

        return nb

    @property
    def getEtaB_3DS_DM_scattering_OOEtauR(self):

        #Define fixed quantities for BEs
        _ETA = [
            np.real(self.epsilon1ab(2,2)),
            np.real(self.epsilon1ab(1,1)),
            np.real(self.epsilon1ab(0,0)),
                    self.epsilon1ab(2,1) ,
                    self.epsilon1ab(2,0) ,
                    self.epsilon1ab(1,0) ,
            np.real(self.epsilon2ab(2,2)),
            np.real(self.epsilon2ab(1,1)),
            np.real(self.epsilon2ab(0,0)),
                    self.epsilon2ab(2,1) ,
                    self.epsilon2ab(2,0) ,
                    self.epsilon2ab(1,0) ,
            np.real(self.epsilon3ab(2,2)),
            np.real(self.epsilon3ab(1,1)),
            np.real(self.epsilon3ab(0,0)),
                    self.epsilon3ab(2,1) ,
                    self.epsilon3ab(2,0) ,
                    self.epsilon3ab(1,0)]

        _C =   [self.c1a(2), self.c1a(1), self.c1a(0),
                self.c2a(2), self.c2a(1), self.c2a(0),
                self.c3a(2), self.c3a(1), self.c3a(0)]

        _K      = [np.real(self.k1), np.real(self.k2), np.real(self.k3)]
        _W      = [ 485e-10*self.MP/self.M1, 1.7e-10*self.MP/self.M1]

        y0      = np.array([0+0j,0+0j,0+0j,0+0j,0+0j,0+0j,0+0j,0+0j,0+0j,0+0j], dtype=np.complex128)

        zcrit   = 1e100
        ys, _   = odeintw(self.RHS_3DS_DM_scattering_OOEtauR, y0, self.xs, args = tuple([_ETA, _C , _K, _W]), full_output=1)
        nb      = np.real(self.sphalfact*(ys[-1,3]+ys[-1,4]+ys[-1,5]))

        return nb

    @property
    def getEtaB_3DS_DM_OOEtauR(self):

        #Define fixed quantities for BEs
        _ETA = [
            np.real(self.epsilon1ab(2,2)),
            np.real(self.epsilon1ab(1,1)),
            np.real(self.epsilon1ab(0,0)),
                    self.epsilon1ab(2,1) ,
                    self.epsilon1ab(2,0) ,
                    self.epsilon1ab(1,0) ,
            np.real(self.epsilon2ab(2,2)),
            np.real(self.epsilon2ab(1,1)),
            np.real(self.epsilon2ab(0,0)),
                    self.epsilon2ab(2,1) ,
                    self.epsilon2ab(2,0) ,
                    self.epsilon2ab(1,0) ,
            np.real(self.epsilon3ab(2,2)),
            np.real(self.epsilon3ab(1,1)),
            np.real(self.epsilon3ab(0,0)),
                    self.epsilon3ab(2,1) ,
                    self.epsilon3ab(2,0) ,
                    self.epsilon3ab(1,0)]

        _C =   [self.c1a(2), self.c1a(1), self.c1a(0),
                self.c2a(2), self.c2a(1), self.c2a(0),
                self.c3a(2), self.c3a(1), self.c3a(0)]

        _K      = [np.real(self.k1), np.real(self.k2), np.real(self.k3)]
        _W      = [ 485e-10*self.MP/self.M1, 1.7e-10*self.MP/self.M1]

        y0      = np.array([0+0j,0+0j,0+0j,0+0j,0+0j,0+0j,0+0j,0+0j,0+0j,0+0j], dtype=np.complex128)

        zcrit   = 1e100
        ys, _   = odeintw(self.RHS_3DS_DM_OOEtauR, y0, self.xs, args = tuple([_ETA, _C , _K, _W]), full_output=1)
        nb      = np.real(self.sphalfact*(ys[-1,3]+ys[-1,4]+ys[-1,5]))

        return nb

    @property
    def getEtaB_3DS_DM_Blanchett(self):

        #Define fixed quantities for BEs
        _ETA = [
            np.real(self.epsilon1ab(2,2)),
            np.real(self.epsilon1ab(1,1)),
            np.real(self.epsilon1ab(0,0)),
                    self.epsilon1ab(2,1) ,
                    self.epsilon1ab(2,0) ,
                    self.epsilon1ab(1,0) ,
            np.real(self.epsilon2ab(2,2)),
            np.real(self.epsilon2ab(1,1)),
            np.real(self.epsilon2ab(0,0)),
                    self.epsilon2ab(2,1) ,
                    self.epsilon2ab(2,0) ,
                    self.epsilon2ab(1,0) ,
            np.real(self.epsilon3ab(2,2)),
            np.real(self.epsilon3ab(1,1)),
            np.real(self.epsilon3ab(0,0)),
                    self.epsilon3ab(2,1) ,
                    self.epsilon3ab(2,0) ,
                    self.epsilon3ab(1,0)]

        _C =   [self.c1a(2), self.c1a(1), self.c1a(0),
                self.c2a(2), self.c2a(1), self.c2a(0),
                self.c3a(2), self.c3a(1), self.c3a(0)]

        _K      = [np.real(self.k1), np.real(self.k2), np.real(self.k3)]
        _W      = [ 485e-10*self.MP/self.M1, 1.7e-10*self.MP/self.M1, 970e-10*self.MP/self.M1]

        y0      = np.array([0+0j,0+0j,0+0j,0+0j,0+0j,0+0j,0+0j,0+0j,0+0j,0+0j,0+0j,0+0j,0+0j], dtype=np.complex128)

        zcrit   = 1e100
        ys, _   = odeintw(self.RHS_3DS_DM_Blanchett, y0, self.xs, args = tuple([_ETA, _C , _K, _W]), full_output=1)
        nb      = np.real(self.sphalfact*(ys[-1,3]+ys[-1,4]+ys[-1,5]))

        return nb

    @property
    def getEtaB_2DS_Approx(self):
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

        ys      = odeintw(self.RHS_2DS_Approx, y0, self.xs, args = tuple([_ETA, _C, _K]))
        # ys      = odeintw(self.RHS_2DS_Approx, y0, xs, args = tuple([_ETA, _HT, _K]))
        # ys      = odeintw(self.RHS_2DS_Approx, y0, xs, args = tuple([_ETA, _HT, _K]))
        nb      = self.sphalfact*(ys[-1,2]+ys[-1,3]+ys[-1,4])

        return nb

    @property
    def getEtaB_2DS_DM(self):

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
        ys, _      = odeintw(self.RHS_2DS_DM, y0, self.xs, args = tuple([_ETA, _C , _K, _W]), full_output=1)
        nb      = np.real(self.sphalfact*(ys[-1,2]+ys[-1,3]+ys[-1,4]))

        return nb

    @property
    def getEtaB_1DS_DM(self):
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

        ys, _      = odeintw(self.RHS_1DS_DM, y0, self.xs, args = tuple(params), full_output=True)
        nb      = self.sphalfact*(ys[-1,1]+ys[-1,2]+ys[-1,3])

        return nb

    @property
    def getEtaB_1DS_DM_scattering_OOEtauR(self):

        #Define fixed quantities for BEs
        _ETA = [
            np.real(self.epsilon1ab(2,2)),
            np.real(self.epsilon1ab(1,1)),
            np.real(self.epsilon1ab(0,0)),
                    self.epsilon1ab(2,1) ,
                    self.epsilon1ab(2,0) ,
                    self.epsilon1ab(1,0) ,
            np.real(self.epsilon2ab(2,2)),
            np.real(self.epsilon2ab(1,1)),
            np.real(self.epsilon2ab(0,0)),
                    self.epsilon2ab(2,1) ,
                    self.epsilon2ab(2,0) ,
                    self.epsilon2ab(1,0) ,
            np.real(self.epsilon3ab(2,2)),
            np.real(self.epsilon3ab(1,1)),
            np.real(self.epsilon3ab(0,0)),
                    self.epsilon3ab(2,1) ,
                    self.epsilon3ab(2,0) ,
                    self.epsilon3ab(1,0)]

        _C =   [self.c1a(2), self.c1a(1), self.c1a(0),
                self.c2a(2), self.c2a(1), self.c2a(0),
                self.c3a(2), self.c3a(1), self.c3a(0)]

        _K      = [np.real(self.k1), np.real(self.k2), np.real(self.k3)]
        _W      = [ 485e-10*self.MP/self.M1, 1.7e-10*self.MP/self.M1]

        y0      = np.array([0+0j,0+0j,0+0j,0+0j,0+0j,0+0j,0+0j,0+0j], dtype=np.complex128)

        zcrit   = 1e100
        ys, _   = odeintw(self.RHS_1DS_DM_scattering_OOEtauR, y0, self.xs, args = tuple([_ETA, _C , _K, _W]), full_output=1)
        nb      = np.real(self.sphalfact*(ys[-1,1]+ys[-1,2]+ys[-1,3]))

        return nb

    @property
    def getEtaB_1DS_DM_allQuantities(self):
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

        xs    = np.linspace(self.xmin, self.xmax, self.xsteps)
        k     = np.real(self.k1)
        y0    = np.array([0+0j,0+0j,0+0j,0+0j,0+0j,0+0j,0+0j], dtype=np.complex128)

        params= np.array([epstt,epsmm,epsee,epstm,epste,epsme,c1t,c1m,c1e,k], dtype=np.complex128)

        ys, _ = odeintw(self.RHS_1DS_DM, y0, self.xs, args = tuple(params), full_output=True)

        return self.sphalfact*(np.array([ys[-1,1],ys[-1,2],ys[-1,3],ys[-1,4],ys[-1,5],ys[-1,6]]))

    @property
    def getEtaB_1DS_DM_ZeroWidthm(self):
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

        ys, _      = odeintw(self.RHS_1DS_DM_ZeroWidthm, y0, self.xs, args = tuple(params), full_output=True)
        nb      = self.sphalfact*(ys[-1,1]+ys[-1,2]+ys[-1,3])

        return nb

    @property
    def getEtaB_1DS_DM_ZeroWidth(self):
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

        ys, _      = odeintw(self.RHS_1DS_DM_ZeroWidth, y0, self.xs, args = tuple(params), full_output=True)
        nb      = self.sphalfact*(ys[-1,1]+ys[-1,2]+ys[-1,3])

        return nb

    @property
    def getEtaB_1DS_Approx(self):
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

        ys      = odeintw(self.RHS_1DS_Approx, y0, self.xs, args = tuple(params))
        nb      = self.sphalfact*(ys[-1,1]+ys[-1,2]+ys[-1,3])

        return nb

      #####   HERE we insert the code for resonant leptogenesis #######
      #CP asymmetry parameter for flavoured resonant leptogenesis from 0705.2183 averaged for large       #washout
    def epsilonaaRES(self,a):
        l         = self.h
        ldag      = np.conjugate(np.transpose(l))
        lcon      = np.conjugate(l)
        M         = self.DM
        lsquare   = np.dot(ldag,l)
        gamma1  = (lsquare[0,0]/(8*np.pi))*M[0,0]
        gamma2  = (lsquare[1,1]/(8*np.pi))*M[1,1]
        DeltaM  = M[1,1]-M[0,0]

        sum1  = np.imag(ldag[0,a]*ldag[0,0]*l[0,1]*l[a,1])
        sum2  = np.imag(ldag[0,a]*ldag[0,1]*l[1,1]*l[a,1])
        sum3  = np.imag(ldag[0,a]*ldag[0,2]*l[2,1]*l[a,1])

        epsbar1  = (sum1+sum2+sum3)/(lsquare[0,0]*lsquare[1,1])
        epsbar2  = (DeltaM/gamma2)/(1+(DeltaM/gamma2)**2)
        epsbar  = - epsbar1 * epsbar2

        return epsbar

    def RHS_2DS_Resonant(self, y0, zzz, ETA, C, K):
        N1, Ntt, Nmm, Nee = y0
        epstt,epsmm,epsee = ETA
        c1t,c1m,c1e = C
        k1term,k2term = K

        d1            = np.real(self.D1(k1term, zzz))
        w1            = np.real(self.W1(k1term, zzz))
        n1eq          = self.N1Eq(zzz)

        c1tc          = np.conjugate(c1t)
        c1mc          = np.conjugate(c1m)
        c1ec          = np.conjugate(c1e)

        #define the different RHSs for each equation
        rhs1           =      -d1*(N1-n1eq)

        rhs2           = (  2 * epstt * d1 * (N1-n1eq)
                                     -  2 * w1 * (2 * c1t * c1tc * Ntt))

        rhs3           = (  2 * epsmm * d1 * (N1-n1eq)
                                     -  2 * w1 * (2 * c1m * c1mc * Nmm))

        rhs4           = (  2 * epsee * d1 * (N1-n1eq)
                                     -  2 * w1 * (2 * c1e * c1ec * Nee))

        RHStemp = [rhs1, rhs2, rhs3, rhs4]
        return RHStemp

    @property
    def getEtaB_2DS_Resonant(self):

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

        ys      = odeintw(self.RHS_2DS_Resonant, y0, self.xs, args = tuple([_ETA, _C, _K]))
        nb      = self.sphalfact*(ys[-1,1]+ys[-1,2]+ys[-1,3])

        return nb


    @property
    def EtaB(self):
        if self.nds==1:
            if self.approx:
                return np.real(self.getEtaB_1DS_Approx)
            else:
                return np.real(self.getEtaB_1DS_DM)
        if self.nds==2:
            if self.approx:
                return np.real(self.getEtaB_2DS_Approx)
            else:
                return np.real(self.getEtaB_2DS_DM)
        if self.nds==3:
            if self.approx:
                return print("We have nocode for 3-flavoured BE")
            else:
                return np.real(self.getEtaB_3DS_DM)
        if self.nds==4:
             if self.approx:
                return 0
             else:
                return np.real(self.getEtaB_3DS_DM_scattering_OOEtauR)
        if self.nds==5:
             if self.approx:
                return 0
             else:
                return np.real(self.getEtaB_1DS_DM_scattering_OOEtauR)
        if self.nds==6:
             if self.approx:
                return 0
             else:
                return np.real(self.getEtaB_3DS_DM_OOEtauR)
        if self.nds==7:
             if self.approx:
                return 0
             else:
                 return np.real(self.getEtaB_2DS_Resonant)
#             else:
#                return np.real(self.getEtaB_1DS_DM_ZeroWidth)

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
    start_time = timeit.default_timer()

    import sys
    from ulysses import readConfig
    _, pdict = readConfig(sys.argv[1])

    L=ULSBase()
    L.setParams(pdict)
    print("Eta_B full =",np.real(L.EtaB))
    L=ULSBase(approx=True)
    L.setParams(pdict)
    print("Eta_B appx =",np.real(L.EtaB))
    L=ULSBase(nds=2, controlplots=False)#, plotprefix=sys.argv[2])
    L.setParams(pdict)
    print("Eta_B full 2DS =",np.real(L.EtaB))
    L=ULSBase(nds=2, approx=True)
    L.setParams(pdict)
    print("Eta_B appx 2DS =",np.real(L.EtaB))
    L=ULSBase(nds=3, controlplots=False)
    L.setParams(pdict)
    print("Eta_B full 3DS =",np.real(L.EtaB))
    L=ULSBase(nds=5, controlplots=False)
    L.setParams(pdict)
    print("Eta_B scattering ('nds' =5)",np.real(L.EtaB))
    from IPython import embed
    embed()

