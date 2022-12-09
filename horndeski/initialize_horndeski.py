"""
initialize_horndeski.py is a Python routine that reads the data
and stores the resulting run variables as pickle files.
The simulations are those of Y. He, A. Roper Pol, A. Brandenburg,
"Modified propagation of gravitational waves from the early radiation
era," submitted to JCAP (2022).

The function run() executes the code. It is only required to be run once
and once the pickle variables have been stored one can directly use the
Jupyter notebook generate_horndeski.ipynb.

Author: Alberto Roper Pol
created: 01/12/2022
"""

import os
import numpy as np

# get working directory, where the runs and routines should be stored
dir0 = os.getcwd() + '/'
HOME = dir0 + '/..'
os.chdir(HOME)

import run as r
from dirs import read_dirs as rd

os.chdir(dir0)

def run(rsd='all'):

    os.chdir(HOME)

    # read the runs
    runs = read_runs(rsd=rsd)

    # save variables
    save_runs(runs)

    os.chdir(dir0)
    
    return runs

def read_runs(rsd='all'):

    """
    Function that reads the runs from the Pencil Code simulations.

    Returns:
        runs -- dictionary that contains the run variables
        rsd -- option that allows to read only one of the four subset of runs
               (default 'all', options 'M0', 'M1', 'M2' or 'M3')
    """

    os.chdir(HOME)

    # dictionary with the name identifying
    # the runs and pointing to the corresponding directory
    dirs = {}
    if rsd == 'all': dirs = rd('horndeski', dirs={})
    if rsd == 'M0': dirs = rd('horndeski_M0', dirs={})
    if rsd == 'M0_lowk': dirs = rd('horndeski_M0_lowk', dirs)
    if rsd == 'M1': dirs = rd('horndeski_M1', dirs={})
    if rsd == 'M2': dirs = rd('horndeski_M2', dirs={})
    if rsd == 'M3': dirs = rd('horndeski_M3', dirs={})
    if rsd == 'M3_lowk': dirs = rd('horndeski_M3_lowk', dirs)
    print('List of runs: ', dirs)
    print('\n')
    R = [s for s in dirs]

    # set quiet to False to see the spectra available, the runs read,
    # and some characteristic info of the run
    runs = r.initialize_runs(R, dir0, dirs, quiet=False, debug=False)
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
