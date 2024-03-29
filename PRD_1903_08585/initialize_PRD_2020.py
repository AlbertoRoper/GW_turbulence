"""
initialize_PRD_2020.py is a Python routine that reads the data, computes
the averaged spectra using the specific time indices for each run at which
the GW spectrum reaches its stationary regime (added by hand), and stores the
resulting run variables as pickle files.
The simulations are those of A. Roper Pol, S. Mandal, A. Brandenburg,
T. Kahniashvili, and A. Kosowsky, "Numerical simulations of gravitational
waves from early-universe turbulence," Phys. Rev. D 102, 083512 (2020),
https://arxiv.org/abs/1903.08585.

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

    os.chdir(HOME)
    # read the runs
    runs = read_runs()

    # runs ini1, ini2, ini3 had an error in the normalization of the time series
    # of the GW energy density EEGW, so we correct it here before
    # saving the new variables, such that the pickle variables already contain
    # the corrected values
    #### ONLY RUN THIS ONCE AFTER READING THE DATA FILES!!
    correct_ini1_3(runs)

    # computes the average and maximum spectra over time and
    # uses the indices of the time array after which the spectra are no
    # longer growing but just oscillating (the specific indices have to be
    # studied separately and one by one for all the runs).
    compute_aver_spec(runs)

    # we can save the variables runs now to avoid repeating the
    # computation of the averaged and maximum spectra
    save_runs(runs)

    os.chdir(dir0)

def read_runs():

    """
    Function that reads the runs from the Pencil Code simulations.

    Returns:
        runs -- dictionary that contains the run variables
    """

    os.chdir(HOME)

    # dictionary with the name identifying
    # the runs and pointing to the corresponding directory
    dirs = rd('PRD_2020_ini')
    dirs = rd('PRD_2020_hel', dirs)
    dirs = rd('PRD_2020_noh', dirs)
    dirs = rd('PRD_2020_ac', dirs)
    R = [s for s in dirs]

    # set quiet to False to see the spectra available, the runs read,
    # and some characteristic info of the run
    runs = r.initialize_runs(R, dir0, dirs, quiet=False)
    r.characterize_runs(runs, quiet=False)

    os.chdir(dir0)

    return runs

def correct_ini1_3(runs):

    """
    Function that corrects the wrong normalization factor in runs ini1, ini2,
    and ini3.
    """

    print('\n')
    print('CORRECTING NORMALIZATION OF EGW IN RUNS ini1, ini2, AND ini3',
          'BY A FACTOR 32 pi/6 \n')

    rrs = ['ini1', 'ini2', 'ini3']
    for i in rrs:
        run = runs.get(i)
        EEGW = run.ts.get('EEGW')
        EEGW *= 16*np.pi/3
        run.ts.update({'EEGW':EEGW})

def compute_aver_spec(runs):

    """
    Function that computes the average and maximum spectra over time and
    uses the indices of the time array after which the spectra are no
    longer growing but just oscillating (the specific indices have to be
    studied separately and one by one for all the runs).

    Arguments:
        runs -- dictionary of the run variables obtained from the Pencil Code
                simulations

    Returns:
        runs -- dictionary of the run variables updated with the maximum and
                averaged spectra of GWs
    """

    print('Computing the maximum value of the spectrum at each wave number',
          ' over times and the averaged spectrum in the oscillatory regime.',
          ' The resulting indices for the current runs are:\n')
    print('   runs: ini1, ini2, ini3, hel1, hel2, hel3, hel4, noh1, noh2,'
          ' ac1, ac2, ac3 \n')
    # once we know the indices corresponding to each run we can run all of them
    indsst = np.array([3, 50, 300, 10, 20, 80, 15, 10, 30, 11, 10, 20])
    print('   indt:', [s for s in indsst], ' \n')

    j = 0
    for i in runs:
        run = runs.get(i)
        t = run.spectra.get('t_GWs')
        indt = indsst[j]
        j += 1
        run.min_max_stat(indt=indt, sp='EGW')
        GWs_stat_sp = run.spectra.get('EGW_stat_sp')
        k = run.spectra.get('k')[1:]
        run.OmGWsat = np.trapz(GWs_stat_sp, k)

def save_runs(runs):

    """
    Function that saves the run variables as pickle variables.

    Arguments:
        runs -- dictionary of the run variables
    """

    for i in runs:
        run = runs.get(i)
        run.save(dir0=dir0)
