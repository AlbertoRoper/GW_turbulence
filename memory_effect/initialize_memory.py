"""
initialize_memory.py is a Python routine that reads the data
and stores the resulting run variables as pickle files.

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
    # save variables
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
    dirs = rd('memory_nonhelical_b73')
    dirs = rd('memory_nonhelical_b27')
    dirs = rd('memory_helical_b73')
    dirs = rd('memory_helical_b27')
    R = [s for s in dirs]

    # set quiet to False to see the spectra available, the runs read,
    # and some characteristic info of the run
    runs = r.initialize_runs(R, dir0, dirs, quiet=False, opt=1)
    r.characterize_runs(runs, quiet=True)

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
