"""
initialize_horndeski.py is a Python routine that reads the data
and stores the resulting run variables as pickle files.
The simulations are those of Y. He, A. Roper Pol, A. Brandenburg,
"Modified propagation of gravitational waves from the early radiation
era," in preparation (2022).

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
import GW_analytical as an
import reading as re
import cosmoGW

os.chdir(dir0)

def run(g=10):

    os.chdir(HOME)

    # read the runs
    runs = read_runs()

    # save variables
    save_runs(runs)

    os.chdir(dir0)
    
    return runs

def read_runs():

    """
    Function that reads the runs from the Pencil Code simulations.

    Returns:
        runs -- dictionary that contains the run variables
    """

    os.chdir(HOME)

    # dictionary with the name identifying
    # the runs and pointing to the corresponding directory
    dirs = {}
#    dirs = rd('horndeski_M0', dirs)
    dirs = rd('horndeski_M1', dirs)
#    dirs = rd('horndeski_M2', dirs)
#    dirs = rd('horndeski_M3', dirs)
#    dirs = rd('horndeski_lowk', dirs=dirs)
    print(dirs)
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
