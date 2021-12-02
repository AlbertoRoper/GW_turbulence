"""
initialize_JCAP_2021.py is a Python routine that reads the data
and stores the resulting run variables as pickle files.
The simulations are those of A. Roper Pol, S. Mandal, A. Brandenburg, and
T. Kahniashvili, "Polarization of gravitational waves from helical MHD
turbulent sources," https://arxiv.org/abs/2107.05356.

The function run() executes the code.
"""

import os
import numpy as np

# get working directory, where the runs and routines should be stored
dir0 = os.getcwd() + '/'
HOME = dir0 + '/..'
os.chdir(HOME)

import run as r
from dirs import read_dirs as rd
import spectra as sp

os.chdir(dir0)

def run():

    # read the runs
    os.chdir(HOME)
    runs = read_runs()
    # compute polarization spectra
    polarization(runs)
    # compute max, min, mean of GW and helical GW spectra
    compute_min_max_hel(runs)
    # assign parameters of the run and some characteristic values
    assign_pars(runs)
    # save variables
    save_runs(runs)
    os.chdir(dir0)

def read_runs():

    """
    Function that reads the runs from the Pencil Code simulations.

    Returns:
        runs -- dictionary that contains the run variables
    """

    # dictionary with the name identifying
    # the runs and pointing to the corresponding directory
    dirs = {}
    dirs = rd('JCAP_2021_ini', dirs)
    dirs = rd('JCAP_2021_dri', dirs)
    R = [s for s in dirs]

    # set quiet to False to see the spectra available, the runs read,
    # and some characteristic info of the run
    runs = r.initialize_runs(R, dir0, dirs, quiet=False)
    r.characterize_runs(runs, quiet=False)

    return runs

def polarization(runs):

    """
    Function that computes the polarization spectra.
    """

    for i in runs:
        run = runs.get(i)
        run.compute_pol()

def compute_min_max_hel(runs):

    """
    Function that computes the average, maximum and minimum spectra over time
    and uses the indices of the time array after which the spectra are no
    longer growing but just oscillating (the specific indices have to be
    studied separately and one by one for all the runs).
    For helical magnetic and GW spectra, it computes the maximum and minimum
    positive and negative values separately.

    Arguments:
        runs -- dictionary of the run variables obtained from the Pencil Code
                simulations

    Returns:
        runs -- dictionary of the run variables updated with the maximum and
                averaged spectra of GWs
    """

    print('Computing the maximum, minimum and average values of the spectra ',
          'at each wave number over times in the oscillatory regime.')

    indsst = []
    tss = []
    for i in runs:
        run = runs.get(i)
        t = run.spectra.get('t_EGW')
        indt = np.argmin(abs(t - max(1 + 1e-1, run.tini)))
        if i == 'f_s05': indt = np.argmin(abs(t - max(1 + 2e-1, run.tini)))
        if i == 'f_s1_neg': indt = np.argmin(abs(t - max(1 + 2e-1, run.tini)))
        if i == 'f_s07': indt = np.argmin(abs(t - max(1 + 3e-1, run.tini)))
        if i == 'f_s03': indt = np.argmin(abs(t - max(1 + 2e-1, run.tini)))
        if i == 'i_s07': indt = np.argmin(abs(t - max(1 + 5e-2, run.tini)))
        if i == 'i_s1': indt = np.argmin(abs(t - max(1 + 1e-1, run.tini)))
        if t[indt] < 1 + 1e-2: indt += 1
        indsst = np.append(indsst, indt)
        tss = np.append(tss, t[indt])
        run.min_max_stat(indt=indt, sp='EGW')
        run.min_max_stat(indt=indt, sp='helEGW', hel=True)
        run.min_max_stat(indt=indt, sp='PGW', abs_b=False)
        run.min_max_stat(indt=indt, sp='Ph', abs_b=False)
    indsst = np.array(indsst)
    tss = np.array(tss)

    print('The resulting indices for the current runs are:\n')
    j = 0
    for i in runs:
        print('time t[%i] = %.4f for run %s'%(indsst[j], tss[j], i))
        j += 1

def assign_pars(runs):

    """
    Function that assigns some parameters of the specific runs.
    """

    for i in runs:
        run = runs.get(i)
        nm = run.name_run
        if 'i' in nm: run.type = 'ini'
        if 'f' in nm: run.type = 'forc'
        if nm == 'i_s01': run.sig = '0.1'
        if '001' in nm:
            if 'neg' in nm: run.sig = '-0.01'
            else: run.sig = '0.01'
        if '03' in nm: run.sig = '0.3'
        if '05' in nm: run.sig = '0.5'
        if '07' in nm: run.sig = '0.7'
        if nm == 'i_s1': run.sig = '1'
        if nm == 'f_s1_neg': run.sig = '-1'
        t = run.spectra.get('t_mag')
        k = run.spectra.get('k')
        indt = np.argmin(abs(t - run.tini))
        EMmax = run.spectra.get('mag')[indt, :]
        EM0 = run.spectra.get('mag')[0, :]
        HM0 = run.spectra.get('helmag_comp')[0, :]
        run.EMmax = np.trapz(EMmax, k)
        run.PM = np.trapz(HM0, k)/np.trapz(EM0, k)
        EGW = run.spectra.get('EGW_stat_sp')
        run.GWstat = np.trapz(EGW, k[1:])
        HGW = run.spectra.get('helEGW_stat_sp')
        helGWstat = np.trapz(HGW, k[1:])
        run.PGW = helGWstat/run.GWstat
        run.n = 1152
        if 'i' in nm: run.eta = 5e-8
        if 'f' in nm: run.eta = 5e-7
        run.k = 600

def save_runs(runs):

    """
    Function that saves the run variables as pickle variables.

    Arguments:
        runs -- dictionary of the run variables
    """

    for i in runs:
        run = runs.get(i)
        run.save(dir0=dir0)
