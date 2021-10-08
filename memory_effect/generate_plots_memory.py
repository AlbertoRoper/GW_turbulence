"""
generate_plots_memory.py is a Python routine  that can be used to generate
the plots of Y. He, A. Brandenburg, and A. Roper Pol, "Leading-order nonlinear
gravitational waves from reheating magnetogeneses".

It reads the pickle run variables that can be generated by the routine
initialize_memory.py.

The function run() executes the code.
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u

# get working directory, where the runs and routines should be stored
dir0 = os.getcwd() + '/'
HOME = dir0 + '/..'
os.chdir(HOME)

from dirs import read_dirs as rd
import run as r
import plot_sets
import spectra
import cosmoGW
import interferometry as inte
import pta

os.chdir(dir0)

def run():

    os.chdir(HOME)
    # import dictionary with the names identifying
    # the runs and pointing to the corresponding directory
    dirs = rd('memory_nonhelical_b73')
    dirs = rd('memory_nonhelical_b27')
    dirs = rd('memory_helical_b73')
    dirs = rd('memory_helical_b27')
    R = [s for s in dirs]

    # read the runs stored in the pickle variables
    runs = r.load_runs(R, dir0, dirs, quiet=False)
    os.chdir(dir0)

    return runs

def select_runs(runs, A='A'):

    """
    Function that returns linear and nonlinear runs corresponding to the
    type of simulations.

    Arguments:
        runs -- variable that contains the memory project runs with the
                stored spectra
        A -- option to chose the type of runs to be plotted (default 'A',
             other options are 'B', 'C', 'D')

    Returns:
        runs_l -- array with linear run variables
        runs_nl -- array with nonlinear run variables
        col -- color corresponding to A for plots
    """

    col = 'blue'
    if A == 'A':
        run1_l = runs.get('A1_l')
        run1_nl = runs.get('A1_nl')
        run2_l = runs.get('A2_l')
        run2_nl = runs.get('A2_nl')
        run3_l = runs.get('A3_l')
        run3_nl = runs.get('A3_nl')
        run4_l = runs.get('A4_l')
        run4_nl = runs.get('A4_nl')

    if A == 'B':
        run1_l = runs.get('B1_l')
        run1_nl = runs.get('B1_nl')
        run2_l = runs.get('B2_l')
        run2_nl = runs.get('B2_nl')
        run3_l = runs.get('B3_l')
        run3_nl = runs.get('B3_nl')
        run4_l = runs.get('B4_l')
        run4_nl = runs.get('B4_nl')
        col = 'darkgreen'

    if A == 'C':
        run1_l = runs.get('C1_l')
        run1_nl = runs.get('C1_nl')
        run2_l = runs.get('C2_l')
        run2_nl = runs.get('C2_nl')
        run3_l = runs.get('C3_l')
        run3_nl = runs.get('C3_nl')
        run4_l = runs.get('C4_l')
        run4_nl = runs.get('C4_nl')
        col = 'orange'

    if A == 'D':
        run1_l = runs.get('D1_l')
        run1_nl = runs.get('D1_nl')
        run2_l = runs.get('D2_l')
        run2_nl = runs.get('D2_nl')
        run3_l = runs.get('D3_l')
        run3_nl = runs.get('D3_nl')
        run4_l = runs.get('D4_l')
        run4_nl = runs.get('D4_nl')
        col = 'red'

    runs_l = [run1_l, run2_l, run3_l, run4_l]
    runs_nl = [run1_nl, run2_nl, run3_nl, run4_nl]

    return runs_l, runs_nl, col

def plot_EGW(runs, A='A', diff=False, save=True):

    """
    Function that plots the resulting GW energy density spectrum at the end
    of inflation (reheating) and compares the result from linear theory to
    the result after adding the leading-order non-linear term (memory effect).

    It generates the plots corresponding to figure 1 of
    Y. He, A. Brandenburg, and A. Roper Pol, "Leading-order nonlinear
    gravitational waves from reheating magnetogeneses".
    It generates the left panels if diff = False and the right panels if
    diff = True

    Arguments:
        runs -- variable that contains the memory project runs with the
                stored spectra
        diff -- option to plot the EGW spectrum or the difference between
                linear and nonlinear when diff = True (default False)
        A -- option to chose the type of runs to be plotted (default 'A',
             other options are 'B', 'C', 'D')
        save -- option to save the plot in plots/EGW_k_'A'_'diff'.pdf
                (default True)
    """

    fig, ax = plt.subplots(figsize=(12, 8))
    plot_sets.axes_lines()
    if diff: ax2 = ax.twinx()
    plot_sets.axes_lines(both=False)

    # chose linear and nonlinear runs corresponding to A
    runs_l, runs_nl, col = select_runs(runs, A=A)
    EEM = [0.02, 0.1, 1, 10]

    if diff:
        for i in range(0, 4):
            run_l = runs_l[i]
            t_l = run_l.spectra.get('t_EGW')
            ind_tl = np.argmin(abs(t_l - 1.))
            if abs(t_l[ind_tl] - 1) > 1e-2:
                print('The time t = 1 is not available in the spectra of',
                      ' the run %s, so t = %.2f has been',
                      ' taken'%(run_l.name_run, t_l[ind_tl]))
            EGW_l = run_l.spectra.get('EGW')[ind_tl, 1:]
            k = run_l.spectra.get('k')[1:]
            run_nl = runs_nl[i]
            t_nl = run_nl.spectra.get('t_EGW')
            ind_tnl = np.argmin(abs(t_nl - 1.))
            if abs(t_nl[ind_tnl] - 1) > 1e-2:
                print('The time t = 1 is not available in the spectra of',
                      ' the run %s, so t = %.2f has been',
                      ' taken'%(run.name_run, t_l[ind_tnl]))
            EGW_nl = run_nl.spectra.get('EGW')[ind_tnl, 1:]
            dif = abs(EGW_nl - EGW_l)
            ax.plot(k, dif, color=col, alpha = .1 + i*.3,
                        label=r'${\cal E}_{\rm EM} = %.2f$'%EEM[i])
            good = np.where(EGW_l != 0)
            ax2.plot(k, dif/EGW_nl[good], '.', color=col,
                     alpha = .15 + i*.15)
            ax2.set_ylim(1e-5, 2.)

    else:
        j = 0
        for i in runs_l:
            t_l = i.spectra.get('t_EGW')
            ind_tl = np.argmin(abs(t_l - 1.))
            if abs(t_l[ind_tl] - 1) > 1e-2:
                print('The time t = 1 is not available in the spectra of',
                      ' the run %s, so t = %.2f has been',
                      ' taken'%(i.name_run, t_l[ind_tl]))
            EGW_l = i.spectra.get('EGW')[ind_tl, 1:]
            k = i.spectra.get('k')[1:]
            ax.plot(k, EGW_l, color=col, alpha = .1 + j*.3,
                     label=r'${\cal E}_{\rm EM} = %.2f$'%EEM[j])
            j += 1
        j = 0
        for i in runs_nl:
            t_nl = i.spectra.get('t_EGW')
            ind_tnl = np.argmin(abs(t_nl - 1.))
            if abs(t_nl[ind_tnl] - 1) > 1e-2:
                print('The time t = 1 is not available in the spectra of',
                      ' the run %s, so t = %.2f has been',
                      ' taken'%(i.name_run, t_nl[ind_tnl]))
            EGW_nl = i.spectra.get('EGW')[ind_tnl, 1:]
            k = i.spectra.get('k')[1:]
            ax.plot(k, EGW_nl, color=col, ls='--', alpha = .1 + j*.3)
            j += 1
        if A == 'A':
            xx = np.linspace(1.2, 5)
            ax.plot(xx, 1e-8*xx, color=col, ls='-.', lw=.8)
            ax.text(2, 1e-10, r'$\sim\!k$', color=col)
            xx = np.linspace(15, 60)
            ax.plot(xx, 1e-12*(xx/10)**(-32), color=col, ls='-.', lw=.8)
            ax.text(18, 1e-30, r'$\sim\!k^{-32}$', color=col)
        if A == 'B':
            xx = np.linspace(1.15, 3)
            ax.plot(xx, 1e-6*xx, color=col, ls='-.', lw=.8)
            ax.text(1.5, 3e-8, r'$\sim\!k$', color=col)
            xx = np.linspace(10, 100)
            ax.plot(xx, 2e-1*(xx/10)**(-10), color=col, ls='-.', lw=.8)
            ax.text(27, 6e-5, r'$\sim\!k^{-10}$', color=col)
        if A == 'C':
            xx = np.linspace(1.4, 20)
            ax.plot(xx, 1e-10*xx**1.5, color=col, ls='-.', lw=.8)
            ax.text(4, 1e-12, r'$\sim\!k^{3/2}$', color=col)
            xx = np.linspace(30, 100)
            ax.plot(xx, 1e10*(xx/10)**(-45), color=col, ls='-.', lw=.8)
            ax.text(30, 1e-24, r'$\sim\!k^{-45}$', color=col)
        if A == 'D':
            xx = np.linspace(1.25, 8)
            ax.plot(xx, 1e-8*xx**1.5, color=col, ls='-.', lw=.8)
            ax.text(2.5, 5e-10, r'$\sim\!k^{3/2}$', color=col)
            xx = np.linspace(11, 50)
            ax.plot(xx, 1e-7*(xx/10)**(-15), color=col, ls='-.', lw=.8)
            ax.text(15, 6e-15, r'$\sim\!k^{-15}$', color=col)

    if not diff:
        ax.legend(fontsize=18, loc='lower left', frameon=False)
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlim(1, 300)
    if diff: ax2.set_yscale('log')
    if A == 'A' or A == 'B': h = 'non-helical'
    else: h = 'helical'
    if A == 'A' or A == 'C':
        b = 7.3
        if diff:
            ax.set_ylim(1e-40, 1e2)
            ax.set_yticks(np.logspace(-46, 2, 13))
        else:
            ax.set_ylim(1e-42, 1e2)
            ax.set_yticks(np.logspace(-42, 2, 12))
    else:
        b = 2.7
        if diff:
            ax.set_ylim(1e-32, 1e2)
            ax.set_yticks(np.logspace(-34, 2, 10))
        else:
            ax.set_ylim(1e-30, 1e2)
            ax.set_yticks(np.logspace(-30, 2, 9))

    ax.set_title(r'%s runs with $\beta = %.1f$'%(h, b), pad=15)
    ax.set_xlabel('$k$')
    if diff:
        ax.set_ylabel(r'$|\Delta E_{\rm GW} (k)|$')
        ax2.set_ylabel(r'$|\Delta E_{\rm GW} (k)|$' + \
                       r'$/E_{\rm GW}^{\rm nlin} (k)$')
    else: ax.set_ylabel(r'$E_{\rm GW} (k)$')
    dff = ''
    if diff: dff = '_diff'
    if save: plt.savefig('plots/' + 'EGW_k_' + A + dff + '.pdf',
                         bbox_inches='tight')

def plot_PGW(runs, A='A', save=True):

    """
    Function that plots the resulting GW polarization spectrum at the end
    of inflation (reheating) and compares the result from linear theory to
    the result after adding the leading-order non-linear term (memory effect).

    It generates the plots corresponding to figure 2 of
    Y. He, A. Brandenburg, and A. Roper Pol, "Leading-order nonlinear
    gravitational waves from reheating magnetogeneses".

    Arguments:
        runs -- variable that contains the memory project runs with the
                stored spectra
        A -- option to chose the type of runs to be plotted (default 'A',
             other options are 'B', 'C', 'D')
        save -- option to save the plot in plots/PGW_k_'A'.pdf (default True)
    """

    plt.figure(figsize=(12, 8))

    # chose linear and nonlinear runs corresponding to A
    runs_l, runs_nl, col = select_runs(runs, A=A)
    EEM = [0.02, 0.1, 1, 10]

    j = 0
    for i in runs_l:
        t_l = i.spectra.get('t_EGW')
        ind_tl = np.argmin(abs(t_l - 1.))
        if abs(t_l[ind_tl] - 1) > 1e-2:
            print('The time t = 1 is not available in the spectra of',
                  ' the run %s, so t = %.2f has been',
                  ' taken'%(i.name_run, t_l[ind_tl]))
        EGW_l = i.spectra.get('EGW')[ind_tl, 1:]
        XiGW_l = i.spectra.get('helEGW')[ind_tl, 1:]
        k = i.spectra.get('k')[1:]
        good = np.where(EGW_l != 0)
        plt.plot(k[good], XiGW_l[good]/EGW_l[good],
                 color=col, alpha=.1 + j*.3,
                 label=r'${\cal E}_{\rm EM} = %.2f$'%EEM[j])
        j += 1

    j = 0
    for i in runs_nl:
        t_nl = i.spectra.get('t_EGW')
        ind_tnl = np.argmin(abs(t_nl - 1.))
        if abs(t_nl[ind_tnl] - 1) > 1e-2:
            print('The time t = 1 is not available in the spectra of',
                  ' the run %s, so t = %.2f has been',
                  ' taken'%(i.name_run, t_nl[ind_tnl]))
        EGW_nl = i.spectra.get('EGW')[ind_tnl, 1:]
        XiGW_nl = i.spectra.get('helEGW')[ind_tnl, 1:]
        k = i.spectra.get('k')[1:]
        good = np.where(EGW_nl != 0)
        plt.plot(k[good], XiGW_nl[good]/EGW_nl[good], '--',
                 color=col, alpha=.1 + j*.3)
        j += 1

    if A == 'A' or A == 'B': h = 'non-helical'
    else: h = 'helical'
    if A == 'A' or A == 'C': b = 7.3
    else: b = 2.7
    plt.title(r'%s runs with $\beta = %.1f$'%(h, b), pad=15)

    plt.xscale('log')
    plt.xlim(1, 300)
    if A=='C':
        plt.ylim(-1, 1.1)
        plt.legend(loc='lower right', fontsize=18, frameon=False)
    if A=='D':
        plt.ylim(-.2, 1.1)
        plt.legend(loc='lower left', fontsize=18, frameon=False)
    plot_sets.axes_lines()
    plt.xlabel('$k$')
    plt.ylabel(r'${\cal P}_{\rm GW} (k)$')
    ax = plt.gca()
    ax.tick_params(axis='x', pad=15)

    if save: plt.savefig('plots/' + 'PGW_k_' + A + '.pdf',
                         bbox_inches='tight')

def plot_OmGW_f(run_l, run_nl, T, g, col='blue'):

    t_l = run_l.spectra.get('t_EGW')
    ind_tl = np.argmin(abs(t_l - 1.))
    if abs(t_l[ind_tl] - 1) > 1e-2:
        print('The time t = 1 is not available in the spectra of',
              ' the run %s, so t = %.2f has been',
              ' taken'%(run_l.name_run, t_l[ind_tl]))
    EGW_l = run_l.spectra.get('EGW')[ind_tl, 1:]
    k_l = run_l.spectra.get('k')[1:]
    f_l, OmGW_l = cosmoGW.shift_OmGW_today(k_l, EGW_l*k_l, T, g)
    t_nl = run_nl.spectra.get('t_EGW')
    ind_tnl = np.argmin(abs(t_nl - 1.))
    if abs(t_nl[ind_tnl] - 1) > 1e-2:
        print('The time t = 1 is not available in the spectra of',
              ' the run %s, so t = %.2f has been',
              ' taken'%(run_nl.name_run, t_nl[ind_tnl]))
    EGW_nl = run_nl.spectra.get('EGW')[ind_tnl, 1:]
    k_nl = run_nl.spectra.get('k')[1:]
    f_nl, OmGW_nl = cosmoGW.shift_OmGW_today(k_nl, EGW_nl*k_nl, T, g)

    plt.plot(f_l, OmGW_nl, color=col)
    plt.plot(f_nl, abs(OmGW_nl - OmGW_l), color=col, ls='-.')

def plot_OmGW_vs_f(runs, save=True):

    """
    Function that plots the resulting GW energy density spectrum at the end
    of inflation (reheating) as an observable at the present time and compares
    it with LISA sensitivity and NANOGrav 12.5 yr results.
    It plots the leading-order nonlinear term as a GW spectrum separately
    for comparison.

    It generates the plots corresponding to figure 4 of
    Y. He, A. Brandenburg, and A. Roper Pol, "Leading-order nonlinear
    gravitational waves from reheating magnetogeneses".

    Arguments:
        runs -- variable that contains the memory project runs with the
                stored spectra
        save -- option to save the plot in plots/OmGW_f_detectors.pdf
                (default True)
    """

    plt.figure(figsize=(12, 8))

    # read LISA PLS
    CWD = os.getcwd()
    os.chdir('..')
    fs, LISA_Om, LISA_OmPLS = inte.read_sens(SNR=10, T=4)
    fs = fs*u.Hz
    os.chdir(CWD)

    # read NANOGrav data
    CWD = os.getcwd()
    os.chdir('../runs_nonhelical_ini')
    _ = pta.read_PTA_data(beta_b=False, Omega_b=False, return_all=True)
    gamma_NG_sPL_1s, A1_NG_sPL_1s, A2_NG_sPL_1s = [_[3], _[4], _[5]]
    betas = np.linspace(-2, 5, 100)
    colors = ['blue']*len(betas)
    _ = pta.CP_delay(betas, colors, obs='NANOGrav_singlePL_1s',
             plot=False)
    fNG = _[0]
    OmGW_NG_a = _[3]
    OmGW_NG_b = _[6]
    _ = pta.CP_delay(betas, colors, obs='NANOGrav_brokenPL_1s',
             plot=False)
    fNGb = _[0]
    OmGW_NGb_a = _[3]
    OmGW_NGb_b = _[6]
    os.chdir(CWD)

    plt.plot(fs, LISA_OmPLS, color='limegreen')
    plt.plot(fs, LISA_Om, color='limegreen', ls='-.')
    minOm, maxOm = pta.get_min_max(fNG, OmGW_NG_a, OmGW_NG_b)
    plt.fill_between(fNG, minOm, maxOm, color='blue', alpha=.3)
    #for i in range(0, len(betas)):
    #    plt.fill_between(fNG, OmGW_NG_a[i, :],
    #                     OmGW_NG_b[i, :], color='blue', alpha=.2)
        #plt.fill_between(fNGb, OmGW_NGb_a[i, :],
        #                 OmGW_NGb_b[i, :], color='blue', alpha=.2)
    plt.text(1e-8, 2e-5, 'NANOGrav 12.5yr', color='blue', fontsize=16)
    plt.text(1e-4, 2e-13, 'LISA PLS', color='limegreen', fontsize=16)

    # chose linear and nonlinear runs corresponding to A
    runsA_l, runsA_nl, colA = select_runs(runs, A='A')
    runsB_l, runsB_nl, colB = select_runs(runs, A='B')
    runsC_l, runsC_nl, colC = select_runs(runs, A='C')
    runsD_l, runsD_nl, colD = select_runs(runs, A='D')
    EEM = [0.02, 0.1, 1, 10]

    # select runs corresponding to EM = 0.1
    runA_l = runsA_l[1]
    runA_nl = runsA_nl[1]
    runB_l = runsB_l[1]
    runB_nl = runsB_nl[1]
    runC_l = runsC_l[1]
    runC_nl = runsC_nl[1]
    runD_l = runsD_l[1]
    runD_nl = runsD_nl[1]

    # Note that T and g are different for every run
    plot_OmGW_f(runA_l, runA_nl, 100*u.GeV, 100, col=colA)
    plot_OmGW_f(runB_l, runB_nl, 150*u.MeV, 15, col=colB)
    plot_OmGW_f(runC_l, runC_nl, 8*u.GeV, 86, col=colC)
    plot_OmGW_f(runD_l, runD_nl, 120*u.MeV, 20, col=colD)

    plot_sets.axes_lines()

    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('$f$ [Hz]')
    plt.ylabel(r'$h_0^2 \Omega_{\rm GW} (f)$')
    plt.xticks(np.logspace(-9, -2, 8))
    plt.xlim(1e-9, 1e-2)
    plt.yticks(np.logspace(-16, -3, 14))
    plt.ylim(1e-16, 1e-3)

    if save: plt.savefig('plots/' + 'OmGW_f_detectors.pdf',
                         bbox_inches='tight')


def generate_table(runs, save=True):

    """
    Function that generates the Table II of Y. He, A. Brandenburg, and
    A. Roper Pol, "Leading-order nonlinear gravitational waves from
    reheating magnetogeneses" that contains the relevant results of the runs.

    Arguments:
        runs -- variable that contains the memory project runs with the
                stored spectra
        save -- option to save the table in tableII.csv
                (default True)
    """

    import pandas as pd

    EGW_ar = []
    DEGW_ar = []
    rat_DEGW_ar = []
    hr_ar = []
    Dhr_ar = []
    rat_Dhr_ar = []
    name = []

    for i in runs:
        run = runs.get(i)
        t = run.ts.get('t')
        indt = np.argmin(abs(t - 1.))
        GW_ts = run.ts.get('EEGW')[indt]
        hr = run.ts.get('hrms')[indt]

        tk = run.spectra.get('t_EGW')
        indt = np.argmin(abs(tk - 1))
        k = run.spectra.get('k')
        GW = np.trapz(run.spectra.get('EGW')[indt, :], k)
        hc = np.trapz(run.spectra.get('GWh')[indt, :], k)

        if '_l' in run.name_run:
            EGW_ar.append(GW)
            hr_ar.append(hr)
            name.append(run.name_run)
            DEGW_ar.append(GW_nl - GW)
            rat_DEGW_ar.append((GW_nl - GW)/GW_nl)
            Dhr_ar.append(hc_nl - hc)
            rat_Dhr_ar.append((hc_nl - hc)/hc_nl)
        else:
            tk = run.spectra.get('t_EGW')
            indt = np.argmin(abs(tk - 1))
            k = run.spectra.get('k')
            GW_nl = np.trapz(run.spectra.get('EGW')[indt, :], k)
            hc_nl = np.trapz(run.spectra.get('GWh')[indt, :], k)

    name = np.array(name)
    EGW_ar = np.array(EGW_ar)
    DEGW_ar = np.array(DEGW_ar)
    rat_DEGW_ar = np.array(rat_DEGW_ar)
    hr_ar = np.array(hr_ar)
    Dhr_ar = np.array(Dhr_ar)
    rat_Dhr_ar = np.array(rat_Dhr_ar)

    df = pd.DataFrame({'name': name, 'EGW': EGW_ar, 'Del EGW': DEGW_ar,
                   'ratio Del EGW': rat_DEGW_ar, 'hrms': hr_ar,
                   'Del hrms': Dhr_ar, 'ratio Del hrms': rat_Dhr_ar})

    if save: df.to_csv('tableII.csv')

    return df
