r"""
This script can be used to compute the tables of the rates that enter the denisty matrix equations 
relevant for low-scale leptogenesis via oscillations. 

The columns of the table correspond to the following quantities:
T       <\gamma_LNC^(0)>/T      <S_LNV^(0)>/T      <\gamma_LNC^(1)>/T      <S_LNV^(2)>/T

The results are stored in a different file named rt_#M#units.txt, where #M is the averaged mass of the heavy neutrinos written with three digits while #units are either MeV, GeV or TeV.
E.g. the rates for M = 1GeV are stored in the file rt_001GeV.txt.
"""


import os
import numpy as np
import ProdRates as PD


HEADER = 'T_GeV <gamma_LNC^(0)>/T <S_LNV^(0)>/T <gamma_LNC^(1)>/T <S_LNV^(2)>/T'


def format_rate_filename(mass_gev):
    if mass_gev < 1:
        mass_value = mass_gev * 1e3
        unit = 'MeV'
    elif mass_gev < 1e3:
        mass_value = mass_gev
        unit = 'GeV'
    else:
        mass_value = mass_gev / 1e3
        unit = 'TeV'

    rounded_mass_value = round(mass_value)
    if not np.isclose(mass_value, rounded_mass_value):
        raise ValueError('Mtest must correspond to an integer number of MeV, GeV, or TeV.')

    return f'rt_{rounded_mass_value:03d}{unit}.txt'


def build_mass_list():
    masses_gev = []

    for mass in [1, 2, 3, 5, 7]:
        masses_gev.append(float(mass))
        masses_gev.append(float(mass) * 1e3)

    for mass in [100, 140, 200, 300, 500, 700]:
        masses_gev.append(float(mass) * 1e-3)
        masses_gev.append(float(mass))
        masses_gev.append(float(mass) * 1e3)

    for mass in [10, 20, 30, 50, 70]:
        masses_gev.append(float(mass))
        masses_gev.append(float(mass) * 1e3)

    return sorted(set(masses_gev))


def compute_rate_table(mass_gev, temperatures_gev):
    gam_lnc0 = np.array([PD.averaged_gamma_p(mass_gev, temperature) for temperature in temperatures_gev])
    s_lnv0   = np.array([PD.averaged_gamma_m(mass_gev, temperature) * temperature**2 / mass_gev**2 for temperature in temperatures_gev])
    gam_lnc1 = np.array([PD.averaged_gamma1_p(mass_gev, temperature, b = 0) for temperature in temperatures_gev])
    s_lnv2   = np.array([PD.averaged_gamma1_m(mass_gev, temperature, b = 0) * temperature**2 / mass_gev**2 for temperature in temperatures_gev])

    return np.column_stack((temperatures_gev, gam_lnc0, s_lnv0, gam_lnc1, s_lnv2))


def save_rate_table(path, mass_gev, temperatures_gev):
    filename   = format_rate_filename(mass_gev)
    rate_table = compute_rate_table(mass_gev, temperatures_gev)
    
    np.savetxt(os.path.join(path, filename), rate_table, header = HEADER, comments = '')


#CHANGE THIS ACCORDINGLY
PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'ARS_rates_ULYSSES')


Npoints = 50
Ts    = np.logspace(7, 2, Npoints)#GeV
MASSES = build_mass_list()

if __name__ == '__main__':
    for mass in [1]:#MASSES:
        save_rate_table(PATH, mass, Ts)
