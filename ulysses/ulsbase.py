from scipy.special import kn
import numpy as np
from odeintw import odeintw
import timeit
from numpy import linalg
import cmath

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
class ULSBase(object):
    def __init__(self, *args, **kwargs):
        r"""
        Base class for all EtaB calculators.

        Set global constants here.

        :Keyword Arguments:
            * *vev* (``float``) --
              Higgs VEV in GeV
            * *mhiggs* (``float``) --
              Higgs mass in GeV
            * *mz* (``float``) --
              Z-boson mass in GeV
            * *gstar* (``float``) --
              Relativistic degrees of freedom at high temperature
            * *mplanck* (``float``) --
              Planck mass in GeV
            * *mstar* (``float``) --
              Neutrino cosmological mass in GeV


        """
        #Higgs vev, mass and Z-mass in GeV
        self.v =  kwargs["vev"]                        if kwargs.get("vev")      is not None else 174.
        self.MH = kwargs["mhiggs"]                     if kwargs.get("mhiggs")   is not None else 125.35
        self.MZ = kwargs["mz"]                         if kwargs.get("mz")       is not None else 91.1876
        #relativistic degrees of freedom at high temperature
        self.gstar = kwargs["gstar"]                   if kwargs.get("gstar")    is not None else 106.75
        #Planck mass in GeV
        self.MP    = kwargs["mplanck"]                 if kwargs.get("mplanck")  is not None else 1.22e+19
        #neutrino cosmological mass in GeV
        self.mstar = kwargs["mstar"]                   if kwargs.get("mstar")    is not None else 1.0e-12
        # Mass-splittings, all in GeV^2
        self.msplit2_solar       = kwargs["m2solar"]   if kwargs.get("m2solar")  is not None else 7.420e-5*1e-18 # 2018
        self.msplit2_athm_normal = kwargs["m2atm"]     if kwargs.get("m2atm")    is not None else 2.515e-3*1e-18 # Values
        self.msplit2_athm_invert = kwargs["m2atminv"]  if kwargs.get("m2atminv") is not None else 2.498e-3*1e-18 # from nu-fit 5.1 WITHOUT SK atmospheric data

        # Flags
        self.debug   = kwargs["debug"]                 if kwargs.get("debug")    is not None else False

        # Parameters of the solver
        self._zmin   = kwargs["zmin"]                  if kwargs.get("zmin")     is not None else 0.1
        self._zmax   = kwargs["zmax"]                  if kwargs.get("zmax")     is not None else 1000
        self._zsteps = kwargs["zsteps"]                if kwargs.get("zsteps")   is not None else 1000
        self._currz  = self.zmin

        # Model switches
        self.ordering = kwargs["ordering"]             if kwargs.get("ordering") is not None else 0
        self.loop     = kwargs["loop"]                 if kwargs.get("loop")     is not None else False

        self._zcut   = kwargs["zcut"]                  if kwargs.get("zcut")     is not None else 1.0 # zcut value for ARS model

        self.zs=None
        self.ys=None
        self.setZS()
        self.normfact = kwargs["normfact"] if kwargs.get("normfact") is not None else 0.013

        self.isCasasIbarrra = True
        self._manualh = np.zeros((3,3), dtype=np.complex128)
        self.pnames = ['m', 'M1', 'M2', 'M3', 'delta', 'a21', 'a31', 'x1', 'x2', 'x3', 'y1', 'y2', 'y3', 't12', 't13', 't23']
        self.evolname="z"


    def shortname(self):
        return ""

    @property
    def flavourindices(self): return None

    @property
    def flavourlabels(self): return None

    @property
    def constants(self):
        s=" Global constants:"
        s+= "\n\t Higgs VEV {} GeV ['vev']".format(self.v)
        s+= "\n\t Higgs mass {} GeV ['mhiggs']".format(self.MH)
        s+= "\n\t Z-boson mass {} GeV ['mz']".format(self.MZ)
        s+= "\n\t Planck mass {} GeV ['mplanck']".format(self.MP)
        s+= "\n\t Neutrino cosmological mass {} GeV ['mstar']".format(self.mstar)
        s+= "\n\t Relativistic degrees of freedom at high temperature {} GeV ['gstar']".format(self.gstar)
        s+= "\n\t Solar mass square splitting {} GeV^2 ['m2solar']".format(self.msplit2_solar)
        s+= "\n\t Atmospheric mass square splitting, normal ordering {} GeV^2 ['m2atm']".format(self.msplit2_athm_normal)
        s+= "\n\t Atmospheric mass square splitting, inverted ordering {} GeV^2 ['m2atminv']".format(self.msplit2_athm_invert)
        return s

    def __call__(self, x):
        r"""
        Operator that returns EtaB for a given parameter point.

        :Arguments:
            * *x* (``dict``) --
              parameter dictionary

        NOTE --- this operator is intended to be used with derived classes where EtaB is implemented
        """
        self.setParams(x)
        return self.EtaB

    def __str__(self):
        s="Model:\n{}".format(self.__doc__)
        s+= "\nNormal ordering\n" if self.ordering==0 else "\nInverted ordering\n"
        s+= "Loop-corrected Yukawa\n" if self.loop else "Tree-level Yukawa\n"
        s+="Integration in [{}, {}] in {} steps\n".format(self._zmin, self._zmax, self._zsteps)
        if self.debug: s+=self.constants
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
    def zcut(self): return self._zcut

    def setEvolData(self, ys):
        self.ys = np.empty((len(self.zs), max(self.flavourindices()) + 2))
        self.ys[:,0] = self.zs
        self.ys[:, self.flavourindices()] = ys[:, self.flavourindices()].real
        self.ys[:,-1] = self.normfact*np.sum(self.ys[:,self.flavourindices()], axis=1)

    def setEvolDataARS(self, ys):
        self.ys = np.empty((len(self.zs), max(self.flavourindices()) + 2))
        self.ys[:,0] = self.zs
        self.ys[:, self.flavourindices()] = ys[:, self.flavourindices()].real
        self.ys[:,-1] = self.normfact*np.sum(self.ys[:,self.flavourindices()], axis=1)


    def setEvolDataPBH(self, ys):
        self.ys = ys
        #self.ys = np.empty((len(self.zs), max(self.flavourindices()) + 2))
        #self.ys[:,0] = self.zs
        #self.ys[:, self.flavourindices()] = ys[:, self.flavourindices()].real
        #self.ys[:,-1] = self.normfact*np.sum(self.ys[:,self.flavourindices()], axis=1)

        
        
    @property
    def evolData(self):
        r"""

        :getter: Return an N-D array of the evolution data.

        The first column is the evolution variable
        The second column corresponds to Ntautau, the
        third to Nmumu and the last columnd to Nee

        """
        if self.flavourindices is not None:
            return self.ys
        else: # FIXME this is only for some compatibility
            pd = np.empty((self.zsteps, 4))
            pd[:,      0] = self.zs
            pd[:,[1,2,3]] = self.ys
        return pd

    def setParams(self, pdict):
        """
        This set the model parameters. pdict is expected to be a dictionary
        """
        if self.isCasasIbarrra:
            self.delta    = pdict['delta']/180*np.pi
            self.a21      = pdict['a21']/180*np.pi
            self.a31      = pdict['a31']/180*np.pi
            self.t12      = pdict['t12']/180*np.pi
            self.t23      = pdict['t23']/180*np.pi
            self.t13      = pdict['t13']/180*np.pi
            self.x1       = pdict['x1']/180*np.pi
            self.y1       = pdict['y1']/180*np.pi
            self.x2       = pdict['x2']/180*np.pi
            self.y2       = pdict['y2']/180*np.pi
            self.x3       = pdict['x3']/180*np.pi
            self.y3       = pdict['y3']/180*np.pi
            self.m        = 10**pdict['m'] * 1e-9 # NOTE input is in log10(m1) in eV --- we convert here to the real value in GeV
            self.M1       = 10**pdict['M1']
            self.M2       = 10**pdict['M2']
            self.M3       = 10**pdict['M3']
        else:
            self.M1       = 10**pdict['M1']
            self.M2       = 10**pdict['M2']
            self.M3       = 10**pdict['M3']
            # Explicit setting of yukawa matix entries
            self._manualh[0][0] = cmath.rect(pdict["Y11_mag"], pdict["Y11_phs"])
            self._manualh[0][1] = cmath.rect(pdict["Y12_mag"], pdict["Y12_phs"])
            self._manualh[0][2] = cmath.rect(pdict["Y13_mag"], pdict["Y13_phs"])
            self._manualh[1][0] = cmath.rect(pdict["Y21_mag"], pdict["Y21_phs"])
            self._manualh[1][1] = cmath.rect(pdict["Y22_mag"], pdict["Y22_phs"])
            self._manualh[1][2] = cmath.rect(pdict["Y23_mag"], pdict["Y23_phs"])
            self._manualh[2][0] = cmath.rect(pdict["Y31_mag"], pdict["Y31_phs"])
            self._manualh[2][1] = cmath.rect(pdict["Y32_mag"], pdict["Y32_phs"])
            self._manualh[2][2] = cmath.rect(pdict["Y33_mag"], pdict["Y33_phs"])




    def printParams(self):
        """
        Print current parameters.
        """
        K = ('delta','a21','a31','t12','t23','t13','x1','y1','x2','y2','x3','y3','m','M1','M2','M3')
        V = (self.delta/np.pi*180, self.a21/np.pi*180, self.a31/np.pi*180, self.t12/np.pi*180, self.t23/np.pi*180, self.t13/np.pi*180, self.x1/np.pi*180, self.y1/np.pi*180, self.x2/np.pi*180, self.y2/np.pi*180, self.x3/np.pi*180, self.y3/np.pi*180,
                np.log10(self.m/1e-9), np.log10(self.M1), np.log10(self.M2), np.log10(self.M3))
        for k, v in zip(K,V):
            print(k,v)


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
                       [0., 0., 1.]], dtype=np.complex128)

        return R1 @ R2 @ R3

    @property
    def SqrtDM(self):
        """
        Square root of diagonal heavy mass matrix.
        TODO: is this np.sqrt(self.DM) ???
        """
        return np.array([[np.sqrt(self.M1), 0., 0.],
                         [0., np.sqrt(self.M2), 0.],
                         [0., 0., np.sqrt(self.M3)]], dtype=np.complex128)

    @property
    def DM(self):
        """
        Diagonal heavy mass matrix.
        """
        return np.array([[self.M1, 0., 0.],
                         [0., self.M2, 0.],
                         [0., 0., self.M3]], dtype=np.complex128)

    @property
    def SqrtDm(self):
        """
        Square root of diagonal light mass matrix.
        Everything is in GeV.
        """

        if self.ordering==0:
            m11 = np.sqrt(self.m)
            m22 = np.sqrt(np.sqrt(self.msplit2_solar       + self.m*self.m))
            m33 = np.sqrt(np.sqrt(self.msplit2_athm_normal + self.m*self.m))
        elif self.ordering==1:
            m11 = np.sqrt(np.sqrt(self.msplit2_athm_invert + self.m*self.m - self.msplit2_solar))
            m22 = np.sqrt(np.sqrt(self.msplit2_athm_invert + self.m*self.m))
            m33 = np.sqrt(self.m)
        else:
            raise Exception("ordering %i not implemented"%self.ordering)

        return np.array([ [m11,  0.,  0.],
                          [ 0., m22,  0.],
                          [ 0.,  0., m33] ], dtype=np.complex128)

    @property
    def U(self):
        """
        PMNS matrix using the PDG parametrisation convention.
        """
        s12     = np.sin(self.t12)
        s23     = np.sin(self.t23)
        s13     = np.sin(self.t13)
        c12     = (1-s12*s12)**0.5
        c23     = (1-s23*s23)**0.5
        c13     = (1-s13*s13)**0.5
        return np.array([ [c12*c13,c13*s12*np.exp(self.a21*1j/2.), s13*np.exp(self.a31*1j/2-self.delta*1j)],
                           [-c23*s12 - c12*np.exp(self.delta*1j)*s13*s23,np.exp((self.a21*1j)/2.)*(c12*c23 - np.exp(self.delta*1j)*s12*s13*s23) , c13*np.exp((self.a31*1j)/2.)*s23],
                           [-c12*c23*np.exp(self.delta*1j)*s13 + s12*s23,np.exp((self.a21*1j)/2.)*(-c23*np.exp(self.delta*1j)*s12*s13 - c12*s23) ,c13*c23*np.exp((self.a31*1j)/2.)]], dtype=np.complex128)


    @property
    def fMR(self):
        """
        This function returns the diagonl heavy neutrino mass matrix taking
        radiative corrections into account. (see h_loop) I.e. equivalent to self.SqrtDM in loop case
        """
        a = 1./self.M1
        b = 1./self.M2
        c = 1./self.M3

        prefactor = -1./(32*np.pi**2*self.v**2)
        d = self.fMLoop(self.M1)
        e = self.fMLoop(self.M2)
        f = self.fMLoop(self.M3)
        A = np.diag([a,b,c])
        B = prefactor*np.diag([d,e,f])

        return np.sqrt(linalg.inv(A+B))


    @property
    def fMLoopHelper(self):
        """
        Helper function.
        """
        prefactor = 1./(32*np.pi**2*self.v**2)
        d = self.fMLoop(self.M1)
        e = self.fMLoop(self.M2)
        f = self.fMLoop(self.M3)
        B = prefactor*np.diag([d,e,f])
        return B

    def fMLoop(self, x):
        """
        The loop function.
        """
        rH2 = (x/self.MH)**2
        rZ2 = (x/self.MZ)**2
        return x*(np.log(rH2)/(rH2-1) + 3*np.log(rZ2)/(rZ2-1) )

    @property
    def m_tree(self):
        """
        Tree-level mass matrix.
        """
        return self.v**2 * self.h @ np.linalg.inv(self.DM) @ np.transpose(self.h)

    @property
    def m_loop(self):
        """
        One-loop-level mass matrix.
        """
        return -1 * self.v**2 * self.h @ self.fMLoopHelper @ np.transpose(self.h)

    @property
    def h(self):
        """
        YUKAWA matrix.
        """
        if self.isCasasIbarrra:
            return self.h_loop if self.loop else self.h_tree
        else:
            return self._manualh

    @property
    def h_loop(self):
        """
        Yukawa matrix (LOOP + Tree).
        """
        return (1./self.v)*(self.U @ self.SqrtDm @ np.transpose(self.R) @ self.fMR)

    @property
    def h_tree(self):
        """
        Yukawa matrix, tree-level.
        """
        return (1./self.v)*(self.U @ self.SqrtDm @ np.transpose(self.R) @ self.SqrtDM)

    @property
    def meff1(self):
        """
        Effective mass 1 used for decay ans washout.
        """
        return np.dot(np.conjugate(np.transpose(self.h)),self.h)[0,0]*(self.v**2)/self.M1

    @property
    def meff2(self,):
        """
        Effective mass 2 used for decay ans washout.
        """
        return np.dot(np.conjugate(np.transpose(self.h)),self.h)[1,1]*(self.v**2)/self.M2

    @property
    def meff3(self,):
        """
        Effective mass 3 used for decay and washout.
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

    def scat(self, z):
        """
        Function that multiplies washouts to incorporate scatterings.
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

    def D2(self, k,z):
        """
        Decay term for Boltzmann equation with two decaying steriles.
        """
        r =self.M2/self.M1
        a = r*r
        x=np.real(r*z)
        b = my_kn1(x)
        c = my_kn2(x)
        return k*z*a*b/c

    def D3(self, k,z):
        """
        Decay term for Boltzmann equation with three decaying steriles.
        """
        r =self.M3/self.M1
        a = r*r
        x=np.real(r*z)
        b = my_kn1(x)
        c = my_kn2(x)
        return k*z*a*b/c

    def NDW1(self, k, z):
        mk1 = my_kn1(z)
        mk2 = my_kn2(z)

        return [3./8*z*z*mk2, k*z*mk1/mk2, 1./4*z*z*z*k*mk1]

    def N1Eq(self, z):
        """
        Equilibrium number density with one decaying sterile.
        """
        n1 = 3./8.*(z**2)*my_kn2(z)
        return n1

    def N2Eq(self, z):
        """
        Equilibrium number density with two decaying steriles.
        """
        r = self.M2/self.M1
        n2 = 3./8.*np.power(r*z,2)*my_kn2(r*z)
        return n2

    def N3Eq(self, z):
        """
        Equilibrium number density with three decaying steriles.
        """
        r = self.M3/self.M1
        n3 = 3./8.*np.power(r*z,2)*my_kn2(r*z)
        return n3

    def W1(self, k1, z):
        """
        Washout parameter with one decaying sterile.
        """
        w1 = 1./4*(z**3)*k1*my_kn1(z)
        return w1

    def W2(self, k, z):
        """
        Washout parameter with two decaying steriles.
        """
        r = self.M2/self.M1
        w2 = k*r/4*np.power(r*z,3) * my_kn1(r*z)
        return w2

    def W3(self, k, z):
        """
        Washout parameter with three decaying steriles.
        """
        r = self.M3/self.M1
        w3 = k*r/4*np.power(r*z,3) * my_kn1(r*z)
        return w3

    def hterm(self, a, b):
        """
        Projection probability projecting onto certain directions
        of flavour space indicated by the indices a and b.

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
        Probability coefficient for 2 a
        """
        norm          = np.sqrt(1./((np.dot(np.conjugate(np.transpose(self.h)), self.h))[1,1]))
        return norm*(self.h[a,1])

    def c3a(self, a):
        """
        Probability coefficient for 3 a
        """
        norm          = np.sqrt(1./((np.dot(np.conjugate(np.transpose(self.h)), self.h))[2,2]))
        return norm*(self.h[a,2])

    def f1(self, x):
        """
        Loop function in epsilon.
        """
        r2=np.power(x,2)

        f1temp = 2./3.*r2*( (1.+r2) * np.log( (1.+r2) / r2 ) - (2.-r2)/(1.-r2) )

        return f1temp if x<10000 else 1

    def f2(self, x):
        """
        Loop function in epsilon.
        """
        return (2./3.)*(1/(np.power(x,2)-1.))

    def epsilon(self, i, j, k, m):
        """
        CP asymmetry parameter.
        i,j,k,m denote indices in the heavy neutrino mass matrix.
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

    def epsilon1ab(self,a,b):
        """
        Off-diagonal CP asymmetry parameter for decays of N1. a and b denote lepton flavour.
        """
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

    def epsilon2ab(self,a,b):
        """
        Off-diagonal CP asymmetry parameter for decays of N2. a and b denote lepton flavour.
        """
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

    def epsilon3ab(self,a,b):
        """
        Off-diagonal CP asymmetry parameter for decays of N3. a and b denote lepton flavour.
        """
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
    def eps100(self): return self.epsilon1ab(0,0)

    @property
    def eps111(self): return self.epsilon1ab(1,1)

    @property
    def eps122(self): return self.epsilon1ab(2,2)

    @property
    def eps200(self): return self.epsilon2ab(0,0)

    @property
    def eps211(self): return self.epsilon2ab(1,1)

    @property
    def eps222(self): return self.epsilon2ab(2,2)

    @property
    def eps300(self): return self.epsilon3ab(0,0)

    @property
    def eps311(self): return self.epsilon3ab(1,1)

    @property
    def eps322(self): return self.epsilon3ab(2,2)

    def epsilonaaRES(self,a):
        """
        CP asymmetry for resonant Leptogenesis.
        """
        l         = self.h
        ldag      = np.conjugate(np.transpose(l))
        lcon      = np.conjugate(l)
        M         = self.DM
        lsquare   = np.dot(ldag,l)
        gamma1    = (lsquare[0,0]/(8*np.pi))*M[0,0]
        gamma2    = (lsquare[1,1]/(8*np.pi))*M[1,1]
        DeltaM    = M[1,1]-M[0,0]
        sum1      = np.imag(ldag[0,a]*l[a,1]*lsquare[0,1]) + (M[0,0]/M[1,1]) *  np.imag(ldag[0,a]*l[a,1]*lsquare[1,0])
        epsbar1   = sum1/(lsquare[0,0]*lsquare[1,1])
        fmix      =  -2*(DeltaM/gamma2)/(1+(2*DeltaM/gamma2)**2)
        fosc      = -0.5*(DeltaM/gamma2)/(1+(DeltaM/gamma2)**2)
        epsbar    = -0.5*epsbar1 * (fosc + fmix)
        return epsbar

    def epsiloniaaRES(self,a,i,j):
        """
        CP asymmetry for resonant leptogenesis in terms of i and j.
        """
        l         = self.h
        ldag      = np.conjugate(np.transpose(l))
        lcon      = np.conjugate(l)
        M         = self.DM
        lsquare   = np.dot(ldag,l)[0:2,0:2]
        detpiece  = np.linalg.det(np.real(lsquare))/(lsquare[0,0]*lsquare[1,1])
        gammai    = (lsquare[i,i]/(8*np.pi))*M[i,i]
        gammaj    = (lsquare[j,j]/(8*np.pi))*M[j,j]
        DeltaM    = M[j,j]-M[i,i]
        sum1      = np.imag(ldag[i,a]*l[a,j]*lsquare[i,j]) + (M[i,i]/M[j,j]) * np.imag(ldag[i,a]*l[a,j]*lsquare[j,i])
        epsbar1   = sum1/(lsquare[0,0]*lsquare[1,1])
        fmix      = ((M[i,i]**2 - M[j,j]**2) * M[i,i] * gammaj)/((M[i,i]**2 - M[j,j]**2)**2 + M[i,i]**2 * gammaj**2)#-2*(DeltaM/gammaj)/(1+(2*DeltaM/gammaj)**2)
        fosc      = ((M[i,i]**2 - M[j,j]**2) * M[i,i] * gammaj) /((M[i,i]**2-M[j,j]**2)**2 +  detpiece * (M[i,i] * gammai + M[j,j] * gammaj)**2)
        epsbar    = -epsbar1 * (fosc + fmix)
        return epsbar
    
    def epsiloniaaRESmix(self,a,i,j):
        """
        CP asymmetry for resonant leptogenesis in terms of i and j.
        """
        l         = self.h
        ldag      = np.conjugate(np.transpose(l))
        lcon      = np.conjugate(l)
        M         = self.DM
        lsquare   = np.dot(ldag,l)[0:2,0:2]
        detpiece  = np.linalg.det(np.real(lsquare))/(lsquare[0,0]*lsquare[1,1])
        gammai    = (lsquare[i,i]/(8*np.pi))*M[i,i]
        gammaj    = (lsquare[j,j]/(8*np.pi))*M[j,j]
        DeltaM    = M[j,j]-M[i,i]
        sum1      = np.imag(ldag[i,a]*l[a,j]*lsquare[i,j]) + (M[i,i]/M[j,j]) * np.imag(ldag[i,a]*l[a,j]*lsquare[j,i])
        epsbar1   = sum1/(lsquare[0,0]*lsquare[1,1])
        fmix      = ((M[i,i]**2 - M[j,j]**2) * M[i,i] * gammaj)/((M[i,i]**2 - M[j,j]**2)**2 + M[i,i]**2 * gammaj**2)#-2*(DeltaM/gammaj)/(1+(2*DeltaM/gammaj)**2)
        epsbar    = -epsbar1 * fmix
        return epsbar

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
    def Gamma1(self):
        """
        Decay rate of N1.
        """
        l         = self.h
        ldag      = np.conjugate(np.transpose(l))
        lcon      = np.conjugate(l)
        M         = self.DM
        lsquare   = np.dot(ldag,l)

        return (M[0,0]/(8*np.pi))*lsquare[0,0]

    @property
    def Gamma2(self):
        """
        Decay rate of N2.
        """
        l         = self.h
        ldag      = np.conjugate(np.transpose(l))
        lcon      = np.conjugate(l)
        M         = self.DM
        lsquare   = np.dot(ldag,l)

        return (M[1,1]/(8*np.pi))*lsquare[1,1]

    @property
    def Gamma3(self):
        """
        Decay rate of N3.
        """
        l         = self.h
        ldag      = np.conjugate(np.transpose(l))
        lcon      = np.conjugate(l)
        M         = self.DM
        lsquare   = np.dot(ldag,l)

        return (M[2,2]/(8*np.pi))*lsquare[2,2]

    @property
    def isPerturbative(self):
        """
        Check for perturbative nature of Yukawas
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


if __name__ == "__main__":
    L=ULSBase()
    print(L)
    L=ULSBase(debug=True)
    print(L)
