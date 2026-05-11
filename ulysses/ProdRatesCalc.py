r"""
This script can be used to compute the tables of the rates that enter the denisty matrix equations 
relevant for low-scale leptogenesis via oscillations. 

Each output table has three columns:
M_GeV   T_GeV   rate

One file is produced per rate quantity:
rates_averaged_G0.txt
rates_averaged_S0.txt
rates_averaged_G1.txt
rates_averaged_S1.txt
rates_averaged_hind_p.txt
rates_averaged_hind_m.txt
"""


import os
import multiprocessing as mp
import numpy as np
import ProdRates as PD


def rate_averaged_G0(mass_gev, temperature_gev):
    return PD.averaged_gamma_p(mass_gev, temperature_gev)

def rate_averaged_S0(mass_gev, temperature_gev):
    return PD.averaged_gamma_m(mass_gev, temperature_gev)

def rate_averaged_G1(mass_gev, temperature_gev):
    return PD.averaged_gamma1_p(mass_gev, temperature_gev, b = 0)

def rate_averaged_S1(mass_gev, temperature_gev):
    return PD.averaged_gamma1_m(mass_gev, temperature_gev, b = 0)

def rate_averaged_h_ind_p(mass_gev, temperature_gev):
    return PD.averaged_h_indirect_p(mass_gev, temperature_gev)

def rate_averaged_h_ind_m(mass_gev, temperature_gev):
    return PD.averaged_h_indirect_m(mass_gev, temperature_gev)


RATE_FUNCTIONS = {
    'G0': rate_averaged_G0,
    'S0': rate_averaged_S0,
    'G1': rate_averaged_G1,
    'S1': rate_averaged_S1,
    'hindp': rate_averaged_h_ind_p,
    'hindm': rate_averaged_h_ind_m,
}


TABLE_SPECS = {
    'G0': ('rates_averaged_G0.txt', 'M_GeV T_GeV <g0>/T'),
    'S0': ('rates_averaged_S0.txt', 'M_GeV T_GeV <S0>/T'),
    'G1': ('rates_averaged_G1.txt', 'M_GeV T_GeV <g1>/T'),
    'S1': ('rates_averaged_S1.txt', 'M_GeV T_GeV <S1>/T'),
    'hindp': ('rates_averaged_hind_p.txt', 'M_GeV T_GeV <hindp>/T'),
    'hindm': ('rates_averaged_hind_m.txt', 'M_GeV T_GeV <hindm>/T'),
}


def _compute_rate_row(task):
    rate_key, mass_gev, temperature_gev = task
    value = RATE_FUNCTIONS[rate_key](mass_gev, temperature_gev)
    return rate_key, mass_gev, temperature_gev, value


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


def save_rate_tables(path, masses_gev, temperatures_gev, n_processes = None, chunksize = 20):
    os.makedirs(path, exist_ok = True)

    if n_processes is None:
        n_processes = max(1, (os.cpu_count() or 1) - 1)#use all but one CPU cores to avoid overloading the system

    tasks = [(rate_key, mass_gev, temperature_gev) for rate_key in TABLE_SPECS  for mass_gev in masses_gev   for temperature_gev in temperatures_gev]

    rows_by_rate = {rate_key: [] for rate_key in TABLE_SPECS}

    with mp.Pool(processes = n_processes) as pool:
        for rate_key, mass_gev, temperature_gev, value in pool.imap_unordered(_compute_rate_row, tasks, chunksize = chunksize):
            rows_by_rate[rate_key].append((mass_gev, temperature_gev, value))

    for rate_key, (filename, header) in TABLE_SPECS.items():
        rows_by_rate[rate_key].sort(key = lambda row: (row[0], row[1]))
        rate_table = np.array(rows_by_rate[rate_key], dtype = float)
        np.savetxt(os.path.join(path, filename), rate_table, header = header, comments = '')


#CHANGE THIS ACCORDINGLY
PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'ARS_rates_ULYSSES')

Npoints = 100
Ts      = np.append(np.logspace(7, 3, Npoints), np.logspace(3, 2, Npoints)) #GeV, more points close to the electroweak phase transition (Tc \sim 160 GeV) where the rates change rapidly
MASSES = build_mass_list()

if __name__ == '__main__':
    save_rate_tables(PATH, MASSES, Ts)
