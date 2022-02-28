"""
initialize_memory.py is a Python routine that reads the data
and stores the resulting run variables as pickle files of Y. He,
A. Roper Pol, and A. Brandenburg, "Leading-order nonlinear gravitational
waves from reheating magnetogeneses".

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
import astropy.units as u

os.chdir(dir0)

def run():

    os.chdir(HOME)

    # read the runs
    runs = read_runs()
    # assign B0
    assign_B0(runs)
    # assign parameters of T scale, helicity and beta parameter to
    # each of the runs
    assign_pars(runs)
    # compute average spectra between t = 2 and t = 10 for toff runs
    compute_stat(runs)
    # save variables
    save_runs(runs)

    os.chdir(dir0)

def read_runs(dirs=[]):

    """
    Function that reads the runs from the Pencil Code simulations.

    Returns:
        runs -- dictionary that contains the run variables
    """

    os.chdir(HOME)

    # dictionary with the name identifying
    # the runs and pointing to the corresponding directory
    if len(dirs) == 0:
        dirs = rd('memory_nonhelical_b73')
        dirs = rd('memory_nonhelical_b27')
        dirs = rd('memory_helical_b73')
        dirs = rd('memory_helical_b27')
        dirs = rd('memory_helical_b17')
        dirs = rd('memory_nonhelical_b17')
        dirs = rd('memory_helical_toff')
        dirs = rd('memory_nonhelical_toff')
    R = [s for s in dirs]

    # set quiet to False to see the spectra available, the runs read,
    # and some characteristic info of the run
    runs = r.initialize_runs(R, dir0, dirs, quiet=False, opt=1)
    r.characterize_runs(runs, quiet=False)

    os.chdir(dir0)

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

def assign_B0(runs):

    for i in runs:
        run = runs.get(i)
        A = run.name_run
        if 'A1' in A: B0 = 3.3e-19
        if 'A2' in A: B0 = 7.5e-19
        if 'A3' in A: B0 = 2.4e-18
        if 'A4' in A: B0 = 7.5e-18
        if 'B1' in A: B0 = 2.7e-7
        if 'B2' in A: B0 = 6e-7
        if 'B3' in A: B0 = 1.9e-6
        if 'B4' in A: B0 = 6e-6
        if 'C1' in A: B0 = 1.7e-24
        if 'C2' in A: B0 = 3.9e-24
        if 'C3' in A: B0 = 1.2e-23
        if 'C4' in A: B0 = 3.9e-23
        if 'D1' in A: B0 = 2.4e-9
        if 'D2' in A: B0 = 5.4e-9
        if 'D3' in A: B0 = 1.7e-8
        if 'D4' in A: B0 = 5.4e-8
        if 'E1' in A: B0 = 4.5e-6
        if 'E2' in A: B0 = 1e-5
        if 'E3' in A: B0 = 3.2e-5
        run.B0 = B0

def assign_pars(runs):

    """
    Function that assigns the temperature, number of degrees of freedom (g),
    helicity (gamma) and beta parameter of the specific runs.
    """

    for i in runs:
        run = runs.get(i)
        pars_n = np.array(['T', 'g', 'gamma', 'beta'])
        if 'A' in run.name_run: pars = np.array([100, 106, 0, 7.3])
        if 'B' in run.name_run: pars = np.array([.15, 15, 0, 2.7])
        if 'C' in run.name_run: pars = np.array([8, 86, 1, 7.3])
        if 'D' in run.name_run: pars = np.array([.12, 20, 1, 2.7])
        if 'E' in run.name_run: pars = np.array([3e5, 106, 1, 1.7])
        run.pars = pars
        run.pars_n = pars_n

def compute_stat(runs):

    """
    Function that computes the average of the GW and helical GW spectra
    over times between 2 and 10 for the additional 't_off' runs.
    """

    for i in runs:
        run = runs.get(i)
        if 'toff' in run.name_run:
            t = run.spectra.get('t_EGW')
            good = np.where(t > 2)[0]
            k = run.spectra.get('k')
            EGW = run.spectra.get('EGW')
            XiGW = run.spectra.get('helEGW')
            EGW_stat = np.trapz(EGW[good, :], t[good], axis=0)/ \
                                    (t[-1] - t[good][0])
            XiGW_stat = np.trapz(XiGW[good, :], t[good], axis=0)/ \
                                    (t[-1] - t[good][0])
            run.spectra.update({'EGW_stat': EGW_stat})
            run.spectra.update({'helEGW_stat': XiGW_stat})
