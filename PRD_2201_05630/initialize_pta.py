"""
initialize_pta.py is a Python routine that reads the data
and stores the resulting run variables as pickle files.
The simulations are those of A. Roper Pol, C. Caprini, A. Neronov, D. Semikoz,
"The gravitational wave signal from primordial magnetic fields in the Pulsar
Timing Array frequency band," https://arxiv.org/abs/2201.05630 (2022).

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

    # save tfin values obtained after analyzing the results
    assign_tfin(runs)

    # compute maximum and mean GW spectra
    compute_min_max_mean(runs)
    # compute GW spectra using analytical model of constant Pi
    compute_Pi0_model(runs, g=g)

    # select runs that can be combined
    for i in runs:
        run = runs.get(i)
        run.comb = False
    combine_runs(runs, 'A1', 'A2')
    combine_runs(runs, 'C1', 'C2')
    combine_runs(runs, 'D1', 'D2')

    # assign some run parameters
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

    os.chdir(HOME)

    # dictionary with the name identifying
    # the runs and pointing to the corresponding directory
    dirs = {}
    dirs = rd('pta', dirs)
    R = [s for s in dirs]

    # set quiet to False to see the spectra available, the runs read,
    # and some characteristic info of the run
    runs = r.initialize_runs(R, dir0, dirs, quiet=False)
    r.characterize_runs(runs, quiet=False)

    os.chdir(dir0)

    return runs

def compute_min_max_mean(runs):

    """
    Function that computes the minimum, maximum, and mean spectra over time.

    Arguments:
        runs -- dictionary of the run variables obtained from the Pencil Code
                simulations

    Returns:
        runs -- dictionary of the run variables updated with the maximum and
                averaged spectra of GWs
    """

    print('Computing the maximum and minimum values of the GW spectrum at each,'
          ' wave number over times and the averaged spectrum.')

    j = 0
    for i in runs:
        print(j, '/%i'%(len(runs) - 1))
        run = runs.get(i)
        run.min_max_stat(sp='EGW')
        EGW_max = run.spectra.get('EGW_max_sp')
        EGW_stat = run.spectra.get('EGW_stat_sp')
        k = run.spectra.get('k')[1:]
        OmGW_max = EGW_max*k
        OmGW_stat = EGW_stat*k
        run.spectra.update({'OmGW_max_sp': OmGW_max})
        run.spectra.update({'OmGW_stat_sp': OmGW_stat})
        j += 1

def compute_Pi0_model(runs, g=10):

    """
    Function that computes the GW spectra of all runs using the analytical model
    that assumes constant Pi (k, t) = Pi(k, t = tini) at initial time.
    It also computes the Pi spectrum using the model given for a Gaussian
    magnetic field, fitted to a smoothed broken power law with smooting
    parameter alpha = 2.

    Arguments:
        runs -- dictionary of the run variables obtained from the Pencil Code
                simulations
    """

    print('Computing the GW spectrum of each run using the analytical model',
          ' that assumes constant Pi equal to its value at the initial time.')

    for i in runs:
        run = runs.get(i)
        run.compute_PiM()
        run.compute_GWsmodel_Pi0()

        # compute GW numerical spectrum at the peak
        OmGW_k = max(run.spectra.get('OmGW_max_sp'))
        run.OmGW_num = cosmoGW.shift_onlyOmGW_today(OmGW_k, g)

        # compute GW analytical envelope from model Pi spectrum, which is
        # computed from the numerical magnetic spectrum, by fitting it to a
        # broken smoothed power law and assuming Gaussianity
        tini = run.tini
        tfin = run.dtfin + tini
        ks = run.spectra.get('ks_Pi_model')
        PiM = run.spectra.get('Pi_model')
        _, OmGW_tfin_k = an.OmGW_from_Pi0(PiM, ks, tfin, tfin=tfin, tini=tini)
        OmGW_tfin = cosmoGW.shift_onlyOmGW_today(OmGW_tfin_k, g)
        # compute the spectral peak
        run.OmGW_env = max(OmGW_tfin[-1, :])
        run.rat_OmGW = run.OmGW_num/run.OmGW_env

def combine_runs(runs, run1, run2):

    """
    Function that combines 2 equivalent runs with different dynamical ranges
    and combines them compensating by the amplitude of the initial magnetic
    spectrum (due to low resolution around the magnetic peak for large domain
    simulations).

    Arguments:
        runs -- dictionary of the run variables obtained from the Pencil Code
                simulations
        run1, run2 -- strings with the names of runs 1 and 2
    """

    print('Combining the results of runs ', run1, ' and ', run2)

    # read spectra of both runs
    rr1 = runs.get(run1)
    k1 = rr1.spectra.get('k')[1:]
    GW1 = rr1.spectra.get('EGW_max_sp')*k1
    mag1 = rr1.spectra.get('mag')[:, 1:]
    rr2 = runs.get(run2)
    k2 = rr2.spectra.get('k')[1:]
    mag2 = rr2.spectra.get('mag')[:, 1:]
    GW2 = rr2.spectra.get('EGW_max_sp')*k2
    # compute factor between magnetic spectra amplitudes
    factM = an.factM(k1, k2, mag1, mag2)
    rr1.spectra.update({'factM': factM})
    # compute combined k and GW spectra
    k_comb, GW_comb = sp.combine(k1, k2, GW1, GW2, factM)
    rr1.spectra.update({'k_comb':k_comb})
    rr1.spectra.update({'EGW_max_sp_comb':GW_comb/k_comb})
    rr1.spectra.update({'OmGW_max_sp_comb':GW_comb})
    rr1.comb = True

def assign_tfin(runs):

    """
    Function that saves the tfin values obtained to characterize the
    numerical simulations using the constant source model.
    This has been obtained after analyzing one by one all of the simulations
    and checking which value of tfin fits better the results.
    """

    for i in runs:
        run = runs.get(i)
        nm = run.name_run
        if 'A' in nm: run.dtfin = .6
        if 'B' in nm: run.dtfin = .6
        if 'C' in nm: run.dtfin = .75
        if 'D' in nm: run.dtfin = .86
        if 'E' in nm: run.dtfin = 2.9

def assign_pars(runs):

    """
    Function that assigns some numerical parameters of the simulations to
    the selected runs.
    """

    for i in runs:
        run = runs.get(i)
        nm = run.name_run
        run.eta = 1e-7
        if 'E' in nm:
            run.n = 512
        else:
            run.n = 768
            if nm == 'A2': run.eta = 1e-6
            if nm == 'B': run.eta = 1e-6
            if nm == 'C1': run.eta = 1e-6
        run.L = re.read_L(dir_data=run.dir_run + '/data')

def save_runs(runs):

    """
    Function that saves the run variables as pickle variables.

    Arguments:
        runs -- dictionary of the run variables
    """

    for i in runs:
        run = runs.get(i)
        run.save(dir0=dir0)
