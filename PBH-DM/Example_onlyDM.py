###################################################################################################
#                                                                                                 #
#                               Primordial Black Hole Evaporation                                 #
#                            Purely Gravitational Interacting Dark Matter.                        #
#                                                                                                 #
#         Authors: Andrew Cheek, Lucien Heurtier, Yuber F. Perez-Gonzalez, Jessica Turner         #
#                                   Based on: arXiv:2107.xxxxx                                    #
#                                                                                                 #
###################################################################################################

#======================================================#
#                                                      #
#                     Example script                   #  
#                                                      #
#======================================================#

import sys
import numpy as np
from odeintw import odeintw
import pandas as pd
from scipy import interpolate
import matplotlib.pyplot as plt
import scipy.integrate as integrate
from scipy.integrate import quad, ode, solve_ivp, odeint
from scipy.optimize import root
from scipy.special import zeta
from scipy.special import kn
import mpmath
from mpmath import polylog

import BHProp as bh

from Omega_h2_onlyDM import FrInPBH

Mi   = 2.      # Log10@ Initial BH mass in grams
ai   = 0.      # Initial a* value, a* = 0. -> Schwarzschild, a* > 0. -> Kerr.
bi   = -17.5  # Log10@beta^\prime 
mDM  = 7.      # Log10 @ DM Mass in GeV
g_DM = 2       # DM degrees of freedom

Oh2 = FrInPBH(Mi, ai, bi, mDM, g_DM) # Omega * h^2

Z = Oh2.Omegah2()

print('{:.7E}'.format(Mi), '{:.7E}'.format(bi), '{:.5E}'.format(ai), '{:.7E}'.format(mDM), Z, sep='\t')

#++++++++++++++++++++++++++++++++++++++++++#
#                                          #
#        Example of scan in DM Mass        #
#                                          #
#++++++++++++++++++++++++++++++++++++++++++#
#   Caveat: Need to comment plotting setup #
#           in the main code               #
#++++++++++++++++++++++++++++++++++++++++++#

# ns = 50

# MDM = np.linspace(-2., 16, num = ns, endpoint=True)

# Z1 = np.zeros((ns))

# for i in range(ns):
    
#     Oh2 = FrInPBH(Mi, ai, bi, MDM[i], g_DM)
    
#     Z1[i] = Oh2.Omegah2()
    
#     print(Mi, MDM[i], Z1[i])

# plt.plot(MDM, Z1, 'o-', markersize=2, label = r"$s_i = 1/2$")
# plt.xlabel(r"$m_{\rm DM}$")
# plt.ylabel(r"$\Omega h^2$")
# plt.legend(loc="upper right")
# plt.title(r"PBHs + DM")
# plt.axhline(y=0.1193, alpha=0.5, color = '#4E2A84', linestyle='--')
# plt.yscale('log')
# plt.show()
