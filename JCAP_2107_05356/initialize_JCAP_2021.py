"""
initialize_JCAP_2021.py is a Python routine that reads the data
and stores the resulting run variables as pickle files.
The simulations are those of A. Roper Pol, S. Mandal, A. Brandenburg, and
T. Kahniashvili, "Polarization of gravitational waves from helical MHD
turbulent sources", https://arxiv.org/abs/2107.05356.

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
    # compute max, min, mean of GW and helical GW spectra
    compute_min_max_hel(runs)
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
    dirs = rd('JCAP_2021_ini')
    dirs = rd('JCAP_2021_dri', dirs)
    R = [s for s in dirs]

    # set quiet to False to see the spectra available, the runs read,
    # and some characteristic info of the run
    runs = r.initialize_runs(R, dir0, dirs, quiet=False)
    r.characterize_runs(runs, quiet=False)

    return runs

def save_runs(runs):

    """
    Function that saves the run variables as pickle variables.

    Arguments:
        runs -- dictionary of the run variables
    """

    for i in runs:
        run = runs.get(i)
        run.save(dir0=dir0)

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
    indsst = np.array(indsst)
    tss = np.array(tss)

    print('The resulting indices for the current runs are:\n')
    j = 0
    for i in runs:
        print('time t[%i] = %.4f for run %s'%(indsst[j], tss[j], i))
        j += 1
