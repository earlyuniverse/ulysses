from scipy.special import kn
import numpy as np
from odeintw import odeintw
import timeit
from numpy import linalg
import cmath

def fast_isPerturbative(y):
    """
    Check perturbativity of Yukawa
    """
    # Limit of perturbativity for y
    limit = np.power(4 * np.pi, 0.5)
    
    # Check if any element of column 1, 2 or 3 is larger than limit
    col1 = (y[0, 0] < limit) * (y[1, 0] < limit) * (y[2, 0] < limit)
    col2 = (y[0, 1] < limit) * (y[1, 1] < limit) * (y[2, 1] < limit)
    col3 = (y[0, 2] < limit) * (y[1, 2] < limit) * (y[2, 2] < limit)
    
    return col1 * col2 * col3

def my_kn1(x):
    """
    Convenience wrapper for kn(1, x)
    """
    return kn(1, x) if x <= 600 else 1e-100


def my_kn2(x):
    """
    Convenience wrapper for kn(2, x)
    """
    return kn(2, x) if x <= 600 else 1e-100


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
        # Higgs vev, mass and Z-mass in GeV
        self.v = kwargs.get("vev", 174.0)
        self.MH = kwargs.get("mhiggs", 125.35)
        self.MZ = kwargs.get("mz", 91.1876)
        
        # Relativistic degrees of freedom at high temperature
        self.gstar = kwargs.get("gstar", 106.75)
        
        # Planck mass in GeV
        self.MP = kwargs.get("mplanck", 1.22e+19)
        
        # Neutrino cosmological mass in GeV
        self.mstar = kwargs.get("mstar", 1.0e-12)
        
        # Mass-splittings, all in GeV^2 -- from NuFit 6.1 (2025), IC23 without SK atm data, best fit
        self.msplit2_solar = kwargs.get("m2solar", 7.537e-5 * 1e-18)
        self.msplit2_athm_normal = kwargs.get("m2atm", 2.521e-3 * 1e-18)
        self.msplit2_athm_invert = kwargs.get("m2atminv", 2.500e-3 * 1e-18)

        # Flags
        self.debug     = kwargs.get("debug", False)
        self.use_hind  = kwargs.get("use_hind", False)

        # Parameters of the z-solver (standard leptogenesis)
        self._zmin = kwargs.get("zmin", 0.001)
        self._zmax = kwargs.get("zmax", 1000)
        self._zsteps = kwargs.get("zsteps", 1000)
        self._currz = self.zmin

        # Parameters of the x-solver (ARS); xmax=None means auto from M1
        self._xmin = kwargs.get("xmin", None)
        self._xmax = kwargs.get("xmax", None)
        self._xsteps = kwargs.get("xsteps", 500)

        # Model switches
        self.ordering = kwargs.get("ordering", 0)
        self.loop = kwargs.get("loop", False)
        self._zcut = kwargs.get("zcut", 1.0)  # zcut value for ARS model
        self.extended_mode = kwargs.get("extended_mode", False)

        # model output path
        self.path = kwargs.get("path", "./")

        self.zs = None
        self.ys = None
        self.setZS()
        self.normfact = kwargs.get("normfact", 0.013)
        #ARS helpers
        self.normfact_ars = (28/79) * np.pi**2 / (27 * 6 * 1.2020569)
        self.plot = None

        # INCLUDEd IN ULYSSESv3
        # Choose the parameterisation of the Yukawa matrix

        # Default parameterisation is CasasIbarra with Euler Angles
        self._manualh = np.zeros((3, 3), dtype=np.complex128)
        self._which_param = None  # must exist before the setter fires
        self.which_param = kwargs.get("which_param", 'euler')  # setter sets pnames
        self.evolname = "z"

        #switch for the initial abundance of RHNs
        self.initial_abundance = kwargs.get("initial_abundance", 0)  #1 for thermal initial abundance, 0 for vanishing initial abundance

        #ARS parameter Lambda which controls the scale at which the fast modes in linear term are forced to enter quasi-static regime
        self.Lambda = kwargs.get("Lambda", None)

        # Storage for the full parameter dict passed to setParams (includes extended params)
        self.pdict = {}


    def shortname(self):
        return ""

    @property
    def which_param(self):
        return self._which_param

    @which_param.setter
    def which_param(self, value):
        if value not in {'euler', 'single_imaginary', 'manual'}:
            raise ValueError("Invalid parameterisation '{}'. Choose from 'euler', 'single_imaginary', or 'manual'.".format(value))
        self._which_param = value
        if value == "euler":
            self.pnames = ['m', 'M1', 'M2', 'M3', 'delta', 'a21', 'a31', 'x1', 'x2', 'x3', 'y1', 'y2', 'y3', 't12', 't13', 't23']
        elif value == "single_imaginary":
            self.pnames = ['m', 'M1', 'M2', 'M3', 'delta', 'a21', 'a31', 'xN1', 'xN2', 'xnu1', 'xnu2', 'x', 'y', 't12', 't13', 't23']
        else:
            self.pnames = ['m', 'M1', 'M2', 'M3',
                           'Y11_mag', 'Y12_mag', 'Y13_mag',
                           'Y21_mag', 'Y22_mag', 'Y23_mag',
                           'Y31_mag', 'Y32_mag', 'Y33_mag',
                           'Y11_phs', 'Y12_phs', 'Y13_phs',
                           'Y21_phs', 'Y22_phs', 'Y23_phs',
                           'Y31_phs', 'Y32_phs', 'Y33_phs']

    def flavourindices(self): return []

    def flavourlabels(self): return []

    def extendedindices(self): return []

    def extendedlabels(self): return []

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

    @property
    def xmin(self): return self._xmin

    @property
    def xmax(self): return self._xmax

    @property
    def xsteps(self): return self._xsteps

    def setEvolData(self, ys):
        fi = self.flavourindices()
        ext = self.extendedindices()
        if not fi:
            return
        self.ys = np.empty((len(self.zs), max(fi) + len(ext) + 3))
        self.ys[:,0] = self.zs
        self.ys[:, fi] = ys[:, fi].real
        if ext:
            self.ys[:, ext] = ys[:, ext].real
        self.ys[:,-2] = ys[:, 0].real  # second to last column: RHN abundance
        self.ys[:,-1] = self.normfact*np.sum(self.ys[:, fi], axis=1)

    def setEvolDataARS(self, ys):
        fi = self.flavourindices()
        ext = self.extendedindices()
        if not fi:
            return
        self.ys = np.empty((len(ys[:,0]), max(fi) + len(ext) + 2))
        self.ys[:,0] = ys[:,0].real
        if ext:
            self.ys[:, ext] = ys[:, ext].real
        self.ys[:, fi] = ys[:, fi].real
        self.ys[:,-1] = self.normfact_ars*np.sum(self.ys[:, fi], axis=1)


    def setEvolDataPBH(self, ys):
        fi  = self.flavourindices()
        ext = self.extendedindices()
        if not fi:
            return
        self.ys = ys


    #CHANGES TO BE INCLUDED in ULYSSESv3
    @property
    def evolData(self):
        r"""

        :getter: Return an N-D array of the evolution data.

        The first column is the evolution variable
        The second column corresponds to Ntautau, the
        third to Nmumu and the last columnd to Nee

        """
        if self.flavourindices():
            return self.ys
        else: # fallback for compatibility
            pd = np.empty((self.zsteps, 4))
            pd[:,      0] = self.zs
            pd[:,[1,2,3]] = self.ys
        return pd

    def setParams(self, pdict):
        """
        Set the model parameters from a dictionary.

        Args:
            pdict (dict): Parameter dictionary containing model parameters.
                In extended_mode the dict may contain extra model-specific keys
                beyond pnames; these are stored in self.pdict for access by the
                subclass EtaB implementation.

        Note:
            Angular parameters are expected in degrees and converted to radians.
            Mass parameters are in specific units (see inline comments).
        """
        self.pdict = pdict
        if "Lambda" in pdict:
            self.Lambda = pdict["Lambda"]
        # Optional mass splittings in eV^2 — converted to GeV^2 internally
        if "m2solar" in pdict:
            self.msplit2_solar = pdict["m2solar"] * 1e-18
        if "m2atm" in pdict:
            self.msplit2_athm_normal = pdict["m2atm"] * 1e-18
        if "m2atminv" in pdict:
            self.msplit2_athm_invert = pdict["m2atminv"] * 1e-18
        # Helper function to convert degrees to radians
        def deg_to_rad(key):
            return pdict[key] / 180 * np.pi
        
        # Helper function to convert log10 masses
        def log10_mass(key):
            return float(10**pdict[key])
        
        if self.which_param == 'euler':
            # CasasIbarra with Euler Angles parameterization
            self._set_common_angles(pdict, deg_to_rad)
            self._set_casas_ibarra_angles(pdict, deg_to_rad)
            self._set_masses(pdict, log10_mass)
            
        elif self.which_param == 'single_imaginary':
            # single_imaginary Casas Ibarra parameterization
            self._set_common_angles(pdict, deg_to_rad)
            self._set_single_imaginary_angles(pdict, deg_to_rad)
            self._set_masses(pdict, log10_mass)
            
        else:
            # Manual Yukawa matrix entries
            self._set_masses(pdict, log10_mass)
            self._set_manual_yukawa_matrix(pdict)
    
    # NuFIT best-fit defaults for PMNS mixing angles (degrees), normal / inverted ordering
    _DEFAULT_T12_NO = 33.76
    _DEFAULT_T13_NO =  8.62
    _DEFAULT_T23_NO = 43.27
    _DEFAULT_T12_IO = 33.76
    _DEFAULT_T13_IO =  8.65
    _DEFAULT_T23_IO = 48.15

    # Parameters recognised in param files but not required (silently defaulted)
    _OPTIONAL_PARAMS = {'t12', 't13', 't23', 'm2solar', 'm2atm', 'm2atminv'}

    def _set_common_angles(self, pdict, deg_to_rad):
        """Set angles common to all parameterisations."""
        self.delta = deg_to_rad('delta')
        self.a21 = deg_to_rad('a21')
        self.a31 = deg_to_rad('a31')
        if self.ordering == 0:
            self.t12 = pdict.get('t12', self._DEFAULT_T12_NO) / 180 * np.pi
            self.t23 = pdict.get('t23', self._DEFAULT_T23_NO) / 180 * np.pi
            self.t13 = pdict.get('t13', self._DEFAULT_T13_NO) / 180 * np.pi
        else:
            self.t12 = pdict.get('t12', self._DEFAULT_T12_IO) / 180 * np.pi
            self.t23 = pdict.get('t23', self._DEFAULT_T23_IO) / 180 * np.pi
            self.t13 = pdict.get('t13', self._DEFAULT_T13_IO) / 180 * np.pi
    
    def _set_casas_ibarra_angles(self, pdict, deg_to_rad):
        """Set angles specific to standard Casas-Ibarra parameterisation."""
        self.x1 = deg_to_rad('x1')
        self.y1 = deg_to_rad('y1')
        self.x2 = deg_to_rad('x2')
        self.y2 = deg_to_rad('y2')
        self.x3 = deg_to_rad('x3')
        self.y3 = deg_to_rad('y3')
    
    def _set_single_imaginary_angles(self, pdict, deg_to_rad):
        """Set angles specific to single_imaginary Casas-Ibarra parameterisation."""
        self.xnu2 = deg_to_rad('xnu2')
        self.xnu1 = deg_to_rad('xnu1')
        self.xN2 = deg_to_rad('xN2')
        self.xN1 = deg_to_rad('xN1')
        self.x = deg_to_rad('x')
        self.y = deg_to_rad('y')
    
    def _set_masses(self, pdict, log10_mass):
        """Set mass parameters."""
        # NOTE: m input is in log10(m) in eV --- we convert here to the real value in GeV
        self.m = log10_mass('m') * 1e-9
        self.M1 = log10_mass('M1')
        self.M2 = log10_mass('M2')
        self.M3 = log10_mass('M3')
    
    def _set_manual_yukawa_matrix(self, pdict):
        """Set Yukawa matrix entries manually using magnitude and phase."""
        yukawa_params = [
            ("Y11_mag", "Y11_phs", 0, 0),
            ("Y12_mag", "Y12_phs", 0, 1),
            ("Y13_mag", "Y13_phs", 0, 2),
            ("Y21_mag", "Y21_phs", 1, 0),
            ("Y22_mag", "Y22_phs", 1, 1),
            ("Y23_mag", "Y23_phs", 1, 2),
            ("Y31_mag", "Y31_phs", 2, 0),
            ("Y32_mag", "Y32_phs", 2, 1),
            ("Y33_mag", "Y33_phs", 2, 2),
        ]
        
        for mag_key, phs_key, i, j in yukawa_params:
            self._manualh[i, j] = cmath.rect(pdict[mag_key], pdict[phs_key])


    # ---------------------------------------------------------------------- #
    #  DM yield utilities — available to any subclass that sets extended_mode #
    # ---------------------------------------------------------------------- #

    # Cosmological constants for DM yield / relic density calculations
    # Entropy density today s0 [cm^-3], rho_c/h^2 [GeV cm^-3] (PDG 2020)
    _s0    = 2891.2
    _rho_c = 1.05371e-5

    @property
    def _ngamma_over_s(self):
        """n_gamma / s = 45 zeta(3) / (pi^4 g*s,0),  g*s,0 = 43/11."""
        return 45.0 * 1.20206 / (np.pi**4 * (43.0 / 11.0))

    def Y_DM(self, N_DM_final):
        """Convert final comoving DM number N_DM to yield Y_DM = n_DM/s."""
        return self._ngamma_over_s * np.real(N_DM_final) / 27.0

    def OmegaDMh2(self, Y_DM, m_dm_GeV):
        """Compute Omega_DM h^2 = m_DM * s0 * Y_DM / (rho_c/h^2)."""
        return m_dm_GeV * self._s0 * Y_DM / self._rho_c

    @property
    def dmYield(self):
        """DM yield Y_DM = n_DM/s after EtaB has been evaluated."""
        return getattr(self, "_Y_DM", None)

    @property
    def OmDMh2(self):
        """Omega_DM h^2 after EtaB has been evaluated."""
        return getattr(self, "_OmDMh2", None)

    def extended_summary(self):
        """Optional extra output printed by uls-calc after eta_B.

        Override in a subclass to return a formatted string of any additional
        quantities computed inside EtaB (e.g. DM yield, relic abundance).
        Return None to suppress extra output (default).
        """
        return None

    # ---------------------------------------------------------------------- #

    def printParams(self):
        """
        Print current parameters.
        """
        K = ('delta','a21','a31','t12','t23','t13','x1','y1','x2','y2','x3','y3','m','M1','M2','M3')
        V = (self.delta/np.pi*180, self.a21/np.pi*180, self.a31/np.pi*180, self.t12/np.pi*180, self.t23/np.pi*180, self.t13/np.pi*180, self.x1/np.pi*180, self.y1/np.pi*180, self.x2/np.pi*180, self.y2/np.pi*180, self.x3/np.pi*180, self.y3/np.pi*180,
                np.log10(self.m/1e-9), np.log10(self.M1), np.log10(self.M2), np.log10(self.M3))
        for k, v in zip(K,V):
            print(k,v)


    @property
    def R_ord(self):
        """
        Orthogonal matrix R = R1.R2.R3
        """

        R23 = np.array([[1., 0., 0.],
                       [0.,  np.cos(self.x1+self.y1*1j), np.sin(self.x1+self.y1*1j)],
                       [0., -np.sin(self.x1+self.y1*1j), np.cos(self.x1+self.y1*1j)]], dtype=np.complex128)

        R13 = np.array([[ np.cos(self.x2+self.y2*1j), 0., np.sin(self.x2+self.y2*1j)],
                       [0., 1. , 0.],
                       [-np.sin(self.x2+self.y2*1j), 0., np.cos(self.x2+self.y2*1j)]], dtype=np.complex128)

        R12 = np.array([[ np.cos(self.x3+self.y3*1j), np.sin(self.x3+self.y3*1j), 0.],
                       [-np.sin(self.x3+self.y3*1j), np.cos(self.x3+self.y3*1j), 0.],
                       [0., 0., 1.]], dtype=np.complex128)

        return R23 @ R13 @ R12

    #single_imaginary R matrix parametrisation as given in arXiv:2106.16226 -- include in ULYSSESv3
    @property
    def R_alt(self):
        """
        Orthogonal matrix R = Onu13.Onu23.RC12.ON23.ON13
        """

        Onu23 = np.array([[1., 0., 0.],
                          [0.,  np.cos(self.xnu1), np.sin(self.xnu1)],
                          [0., -np.sin(self.xnu1), np.cos(self.xnu1)]], dtype=np.complex128)

        Onu13 = np.array([[ np.cos(self.xnu2), 0., np.sin(self.xnu2)],
                          [0., 1. , 0.],
                          [-np.sin(self.xnu2), 0., np.cos(self.xnu2)]], dtype=np.complex128)
        
        ON23 = np.array([[1., 0., 0.],
                         [0.,  np.cos(self.xN1), np.sin(self.xN1)],
                         [0., -np.sin(self.xN1), np.cos(self.xN1)]], dtype=np.complex128)

        ON13 = np.array([[ np.cos(self.xN2), 0., np.sin(self.xN2)],
                         [0., 1. , 0.],
                         [-np.sin(self.xN2), 0., np.cos(self.xN2)]], dtype=np.complex128)

        RC12 = np.array([[ np.cos(self.x+self.y*1j), np.sin(self.x+self.y*1j), 0.],
                         [-np.sin(self.x+self.y*1j), np.cos(self.x+self.y*1j), 0.],
                         [0., 0., 1.]], dtype=np.complex128)
        

        return np.transpose(Onu13 @ Onu23 @ RC12 @ ON23 @ ON13)
    
    
    @property
    def R(self):
        if self.which_param == 'single_imaginary':
            return self.R_alt
        elif self.which_param == 'CP_conserving':
            return self.R_CP_cons
        else:
            return self.R_ord

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
    def Theta(self):
        """
        Mixing of the heavy Majorana neutrinos with SM charged leptons, aka mixing parameter
        """
        Theta_mat = np.zeros((3,3), dtype=np.complex128)
        for a in range(3):
            for j in range(3):
                Theta_mat[a,j] = (self.v)*self.h[a, j]/self.DM[j,j]
        return Theta_mat

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
        return -1 * self.v**2 * self.h @ np.linalg.inv(self.DM) @ np.transpose(self.h)

    @property
    def m_loop(self):
        """
        One-loop-level mass matrix.
        
        """
        print(self.fMLoopHelper)
        return -1 * self.v**2 * self.h @ self.fMLoopHelper @ np.transpose(self.h)

    @property
    def h(self):
        """
        YUKAWA matrix.
        """
        if self.which_param in {'euler', 'single_imaginary', 'CP_conserving'}:
            #Some test Yukawas
            h_test = -1j*np.array([[ 3.10409870e-04-2.88390858e-06j, -3.15972857e-04+2.93988151e-06j,-4.11823008e-06-4.42936908e-04j],
                                   [ 5.63377207e-04-1.13668872e-06j, -5.73902621e-04+1.15657100e-06j,-1.62164031e-06-8.04212584e-04j],
                                   [-4.11904691e-07-9.81246627e-07j, -1.60266127e-07+9.98410006e-07j,-1.39988112e-06+1.74359234e-07j]])
            #BMCI from 1810.12463, Inverted Ordering
            h_BMCI = np.array([[-2.0e-5 - 1j*7.9e-5, 7.9e-5 - 1j*2.0e-5, 1.8e-8 - 1j*9.5e-8],
                       [2.7e-5 - 1j*1.3e-5, 1.3e-5 + 1j*2.7e-5, 4.6e-8 - 1j*2.8e-8],
                       [-2.9e-5 - 1j*0.4e-5, 0.4e-5 - 1j*2.9e-5, -4.0e-8 + 1j*1.0e-8]])
            ##BMCII from 1810.12463, Normal Ordering
            h_BMCII = np.array([[2.8e-5 - 1j*0.4e-5, 0.4e-5 + 1j*2.8e-5, 0.4e-8 - 1j*1.3e-8],
                        [-3.0e-7 - 1j*0.8e-7, 0.8e-7 - 1j*3.0e-7, -3.9e-8 + 1j*5.7e-8],
                        [-4.9e-5 + 1j*0.4e-5, -0.4e-5 - 1j*4.9e-5, -3.4e-8 + 1j*5.0e-8]])
            #BMCIII from 1810.12463, Normal Ordering
            h_BMCIII = np.array([[-3.2e-8 - 1j*4.5e-8, -1.1e-7 - 1j*1.7e-7, -2.4e-7 + 1j*1.6e-7],
                        [1.7e-7 - 1j*5.9e-7, 0.6e-6 - 1j*2.1e-6, -2.9e-6 - 1j*0.8e-6],
                        [4.4e-7 - 1j*3.0e-7, 1.5e-6 - 1j*1.1e-6, -1.5e-6 - 1j*2.1e-6]])
            
            return self.h_loop if self.loop else self.h_tree
        else:
            return self._manualh

    @property
    def h_loop(self):
        """
        Yukawa matrix (LOOP + Tree).
        """
        return (1j/self.v)*(self.U @ self.SqrtDm @ np.transpose(self.R) @ self.fMR)

    @property
    def h_tree(self):
        """
        Yukawa matrix, tree-level.
        """
        return (1j/self.v)*(self.U @ self.SqrtDm @ np.transpose(self.R) @ self.SqrtDM)

    #Include in ULYSSESv3, the heavy Majorana neutrino mixing
    @property
    def RV(self):
        """
        Heavy Majorana neutrino mixing
        """
        return 1j*self.U @ self.SqrtDm @ np.transpose(self.R) @ linalg.inv(self.SqrtDM)
    
    
    @property
    def eta(self):
        """
        Heavy Majorana neutrino mixing
        """
        return -(1/2)*self.RV@np.transpose(np.conjugate(self.RV))
    
    @property
    def h_non_unitary(self):
        """
        Yukawa matrix, tree-level.
        """
        return self.h_tree + (1./self.v)*(self.eta@self.U @ self.SqrtDm @ np.transpose(self.R) @ self.SqrtDM)
    
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
    def zosc_mat(self):
        Tew    = 131.7
        M0     = 7.112582895088419e+17
        zosc_mat = np.zeros((3,3), dtype=np.complex128)
        deltaM2_mat = np.zeros((3,3), dtype=np.complex128)
        for i in range(3):
            for j in range(3):
                deltaM2_mat[i,j] = self.DM[j,j]**2 - self.DM[i,i]**2
                zosc_mat[i,j] = np.cbrt(np.real(12.*Tew**3/(deltaM2_mat[i,j] * M0)))
        return np.real(zosc_mat)

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
        DStemp    = 0.11399*k1*(1+t*(z**2)*np.log(1+8.77298/(t*z)))
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
    
    def epsilon1ab_reg(self,a,b):
        """
        Off-diagonal CP asymmetry parameter for decays of N1. a and b denote lepton flavour. It is regularised to
        include the resonant effects between N1 and N2.
        """
        l         = self.h
        ldag      = np.conjugate(np.transpose(l))
        lcon      = np.conjugate(l)
        M         = self.DM
        lsquare   = np.dot(ldag,l)
        
        DeltaM21 = M[1,1]-M[0,0]
        DeltaM31 = M[2,2]-M[0,0]
        x21 = (M[1,1]/M[0,0])**2
        x31 = (M[2,2]/M[0,0])**2
        reg21 = 1/(1 + self.Gamma2**2 * (1/(DeltaM21 * (1 + x21))**2))
        reg31 = 1/(1 + self.Gamma3**2 * (1/(DeltaM31 * (1 + x31))**2))
        if x21<10000:
            f21 = 2./3.* x21 *( (1.+x21) * np.log( (1.+x21) / x21 ) - 1 + reg21 * (1.)/(x21 - 1.) )
        else:
            f21 = 1
        if x31<10000:
            f31 = 2./3.* x31 *( (1.+x31) * np.log( (1.+x31) / x31 ) - 1 + reg31 * (1.)/(x31 - 1.) )
        else:
            f31 = 1

        prefactor   = (3/(32*np.pi))*(1/(lsquare[0,0]))
        first       = 1j*(lsquare[1,0]*l[a,0]*lcon[b,1]-lsquare[0,1]*l[a,1]*lcon[b,0])*(M[0,0]/M[1,1])*f21
        second      = 1j*(2./3.)*(reg21/(np.power(M[1,1]/M[0,0],2)-1.))*(l[a,0]*lcon[b,1]*lsquare[0,1]-lcon[b,0]*l[a,1]*lsquare[1,0])

        third       = 1j*(lsquare[2,0]*l[a,0]*lcon[b,2]-lsquare[0,2]*l[a,2]*lcon[b,0])*(M[0,0]/M[2,2])*f31
        fourth      = 1j*(2./3.)*(reg31/(np.power(M[2,2]/M[0,0],2)-1.))*(l[a,0]*lcon[b,2]*lsquare[0,2]-lcon[b,0]*l[a,2]*lsquare[2,0])
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
    
    
    def epsilon2ab_reg(self,a,b):
        """
        Off-diagonal CP asymmetry parameter for decays of N2. a and b denote lepton flavour. It is regularised to
        include the resonant effects between N2 and N3
        """
        l         = self.h
        ldag      = np.conjugate(np.transpose(l))
        lcon      = np.conjugate(l)
        M         = self.DM
        lsquare   = np.dot(ldag,l)
        
        DeltaM12 = M[0,0]-M[1,1]
        DeltaM32 = M[2,2]-M[1,1]
        x12 = (M[0,0]/M[1,1])**2
        x32 = (M[2,2]/M[1,1])**2
        reg12 = 1/(1 + self.Gamma1**2 * (1/(DeltaM12 * (1 + x12))**2))
        reg32 = 1/(1 + self.Gamma3**2 * (1/(DeltaM32 * (1 + x32))**2))
        #Bad behaviour in the limit of large xk/xj
        if x12 < 10000:
            f12 = 2./3.* x12 *( (1.+x12) * np.log( (1.+x12) / x12 ) - 1 + reg12 * (1.)/(x12 - 1.) )
        else:
            f12 = 1
        if x32 < 10000:
            f32 = 2./3.* x32 *( (1.+x32) * np.log( (1.+x32) / x32 ) - 1 + reg32 * (1.)/(x32 - 1.) )
        else:
            f32 = 1
    

        #define terms of epsilon: prefactor and first term (first), second term (second) etc.
        prefactor   = (3/(32*np.pi))*(1/(lsquare[1,1]))
        first       = 1j*(lsquare[0,1]*l[a,1]*lcon[b,0]-lsquare[1,0]*l[a,0]*lcon[b,1])*(M[1,1]/M[0,0])*f12
        second      = 1j*(2./3.)*(reg12/(np.power(M[0,0]/M[1,1],2)-1.))*(l[a,1]*lcon[b,0]*lsquare[1,0]-lcon[b,1]*l[a,0]*lsquare[0,1])
        
        #regularised
        third       = 1j*(lsquare[2,1]*l[a,1]*lcon[b,2]-lsquare[1,2]*l[a,2]*lcon[b,1])*(M[1,1]/M[2,2])*f32
        fourth      = 1j*(2./3.)*(reg32/(np.power(M[2,2]/M[1,1],2)-1.))*(l[a,1]*lcon[b,2]*lsquare[1,2]-lcon[b,1]*l[a,2]*lsquare[2,1])
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
    
    
    def epsilon3ab_reg(self,a,b):
        """
        Off-diagonal CP asymmetry parameter for decays of N3. a and b denote lepton flavour. It is regularised to
        include the resonant effects in the case of M1 < M2 sim M3
        """
        l         = self.h
        ldag      = np.conjugate(np.transpose(l))
        lcon      = np.conjugate(l)
        M         = self.DM
        lsquare   = np.dot(ldag,l)

        DeltaM13 = M[0,0]-M[2,2]        
        DeltaM23 = M[1,1]-M[2,2]
        x13 = (M[0,0]/M[2,2])**2
        x23 = (M[1,1]/M[2,2])**2
        reg13 = 1/(1 + self.Gamma1**2 * (1/(DeltaM13 * (1 + x13))**2))
        reg23 = 1/(1 + self.Gamma2**2 * (1/(DeltaM23 * (1 + x23))**2))
        if x13 < 10000:
            f13 = 2./3.* x13 *( (1.+x13) * np.log( (1.+x13) / x13 ) - 1 + reg13 * (1.)/(x13 - 1.) )
        else:
            f13 = 1
        if x23 < 10000:
            f23 = 2./3.* x23 *( (1.+x23) * np.log( (1.+x23) / x23 ) - 1 + reg23 * (1.)/(x23 - 1.) )
        else:
            f23 = 1

        #define terms of epsilon: prefactor and first term (first), second term (second) etc.
        prefactor   = (3/(32*np.pi))*(1/(lsquare[2,2]))
        first       = 1j*(lsquare[0,2]*l[a,2]*lcon[b,0]-lsquare[2,0]*l[a,0]*lcon[b,2])*(M[2,2]/M[0,0])*f13
        second      = 1j*(2./3.)*(reg13/(np.power(M[0,0]/M[2,2],2)-1.))*(l[a,2]*lcon[b,0]*lsquare[2,0]-lcon[b,2]*l[a,0]*lsquare[0,2])
        
        #regularised
        third       = 1j*(lsquare[1,2]*l[a,2]*lcon[b,1]-lsquare[2,1]*l[a,1]*lcon[b,2])*(M[2,2]/M[1,1])*f23
        fourth      = 1j*(2./3.)*(reg23/(np.power(M[1,1]/M[2,2],2)-1.))*(l[a,2]*lcon[b,1]*lsquare[2,1]-lcon[b,2]*l[a,1]*lsquare[1,2])
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
        return Gamma/(self.M2-self.M1)

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
