"""
etab1BE1F_DM_FreezeIn.py — Freeze-in dark matter + leptogenesis.

Extends the standard 1BE1F leptogenesis model with a third Boltzmann equation tracking a dark matter particle X produced
via freeze-in from N1 decays through a dark sector coupling lambda.


    dN_DM/dz   =  Br_dark * D(k,z) * (N1 - N1eq)               [new]

where Br_dark = (Gamma_dark / Gamma_total) is the branching ratio of N1
into the dark sector.  In the FIMP limit the DM is never in equilibrium and
its equilibrium abundance is effectively zero, so no washout term appears.

The dark partial width is parametrised as

    Gamma_dark = lambda^2 * M1 / (8*pi)

Runcard additions
-----------------
Two new parameters are read from the parameter card alongside the standard
1BE1F parameters:

    lam     <value>    # dark sector coupling  (dimensionless, e.g. 1e-8)
    m_dm    <value>    # DM mass in GeV        (e.g. 1e3)

If absent, defaults of lam=0 and m_dm=M1/10 are used (safe fallback).

Outputs
-------
N_B-L  — standard B-L number density 
EtaB   — standard baryon asymmetry eta_B   
N_DM   — Dark Matter number density
"""

import ulysses
import numpy as np
from odeintw import odeintw
from numba import jit
from scipy.special import kn


# ---------------------------------------------------------------------- #
#  JIT-compiled RHS — three coupled equations                            #
# ---------------------------------------------------------------------- #

@jit(nopython=True)
def fast_RHS_FreezeIn(y0, d_sm, d_dark, w1, n1eq, epstt, epsmm, epsee):
    """
    Right-hand side for freeze-in leptogenesis + DM system.
    """
    N1   = y0[0]
    NBL  = y0[1]
    N_DM = y0[2]

    dN1   = -(d_sm + d_dark) * (N1 - n1eq)   # total width depletes N1
    dNBL  = (epstt + epsmm + epsee) * d_sm * (N1 - n1eq) - w1 * NBL
    dNDM  =  d_dark * N1   # FIMP: only forward decays, no inverse process

    return [dN1, dNBL, dNDM]


# ---------------------------------------------------------------------- #
#  Model class                                                            #
# ---------------------------------------------------------------------- #

class EtaB_1BE1F_DM_FreezeIn(ulysses.ULSBase):

    def shortname(self): return "1BE1F_DM_FreezeIn"

    def flavourindices(self): return [1]

    def flavourlabels(self): return ["$N_{B-L}$"]

    def extendedindices(self): return [2, 3]

    def extendedlabels(self): return [r"$N_{DM}$", r"$N_{DM}^{\rm eq}$"]

    def extended_summary(self):
        if self.dmYield is None:
            return None
        return "{:<14}{}\n{:<14}{}".format(
            "Y_DM",         self.dmYield,
            "Omega_DM h^2", self.OmDMh2)

    # ------------------------------------------------------------------ #
    #  Dark sector decay parameter                                         #
    # ------------------------------------------------------------------ #

    def _Br_dark(self):
        """
        Branching ratio of N1 into the dark sector.

        Gamma_dark = lambda^2 * M1 / (8*pi)
        Gamma_SM   is encoded via the decay parameter k1 = Gamma_SM / H(z=1)

        Br_dark = (lam^2 / (8*pi)) / (tilde_m1 / v^2 + lam^2 / (8*pi))
        k_dark = lambda^2 * M1 / (8*pi*H*)

        so  Br_dark = k_dark / (k1 + k_dark).
        """
        lam   = float(self.pdict.get("lam", 1.0))
        M1    = self.M1   # GeV
        # H* = H(T = M1) in the radiation dominated era
        # H* = 1.66 * sqrt(g*) * M1^2 / M_Pl  (standard formula)
        gstar = 106.75
        M_Pl  = 1.22e19  # GeV
        Hstar = 1.66 * np.sqrt(gstar) * M1**2 / M_Pl

        k_dark = lam**2 * M1 / (8.0 * np.pi * Hstar)
        k1     = np.real(self.k1)

        denom = k1 + k_dark
        if denom == 0.0:
            return 0.0
        return k_dark / denom

    def RHS(self, y0, z, epstt, epsmm, epsee, k, br_dark):
        """Evaluate the RHS at redshift parameter z."""
        if z != self._currz or z == self.zmin:
            self._d    = np.real(self.D1(k, z))
            self._w1   = np.real(self.W1(k, z))
            self._n1eq = self.N1Eq(z)
            self._currz = z

        # Split total decay parameter into SM and dark contributions.
        # D is proportional to Gamma_total; we need to rescale.
        # d_sm   drives leptogenesis  (fraction 1 - Br_dark)
        # d_dark drives DM production (fraction Br_dark)
        d_total = self._d
        d_sm    = d_total * (1.0 - br_dark)
        d_dark  = d_total * br_dark

        return fast_RHS_FreezeIn(
            y0, d_sm, d_dark, self._w1, self._n1eq,
            epstt, epsmm, epsee
        )

    @property
    def EtaB(self):
        """
        Solve the coupled BE system and return the baryon asymmetry eta_B.

        Also stores self._Y_DM and self._OmDMh2 for retrieval after the call.
        """
        # Standard CP asymmetry parameters (same as 1BE1F)
        epstt = np.real(self.epsilon1ab(2, 2))
        epsmm = np.real(self.epsilon1ab(1, 1))
        epsee = np.real(self.epsilon1ab(0, 0))

        k      = np.real(self.k1)
        br_dark = self._Br_dark()

        # Initial conditions: N1 follows initial_abundance, asymmetry and N_DM start at zero
        y0 = np.array([self.initial_abundance * self.N1Eq(self.zmin) + 0j, 0. + 0j, 0. + 0j], dtype=np.complex128)

        params = np.array([epstt, epsmm, epsee, k, br_dark],
                          dtype=np.complex128)

        ys = odeintw(
            self.RHS, y0, self.zs, args=tuple(params),
            atol=1e-14, rtol=1e-10
        )

        # DM mass — read from runcard, fallback to M1/10
        m_dm_GeV = float(self.pdict.get("m_dm", self.M1 / 10.0))

        # Compute equilibrium number density for DM at each z
        # N_DM_eq(z) = (3/8) * (m_dm/M1 * z)^2 * K2(m_dm/M1 * z)
        # same normalization as N1Eq; in FIMP regime N_DM << N_DM_eq throughout
        r = m_dm_GeV / self.M1
        x_vals = r * self.zs
        x_safe = np.where(x_vals < 1e-10, 1e-10, x_vals)
        N_DM_eq_vals = np.where(
            x_vals < 1e-10,
            0.75,  # relativistic limit: (3/8) * x^2 * (2/x^2) = 3/4
            (3.0 / 8.0) * x_safe**2 * kn(2, x_safe)
        )

        ys_ext = np.column_stack([np.real(ys), N_DM_eq_vals])
        self.setEvolData(ys_ext)

        # DM outputs (stored for external access)
        N_DM_final   = ys[-1][2]
        self._Y_DM   = self.Y_DM(N_DM_final)
        self._OmDMh2 = self.OmegaDMh2(self._Y_DM, m_dm_GeV)

        # Return standard baryon asymmetry (last column of self.ys)
        return self.ys[-1][-1]

