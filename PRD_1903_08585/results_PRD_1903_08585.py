"""
results_PRD_1903_08585.py is a Python routine that reads the data, computes
the averaged spectra using the specific time indices for each run at which
the GW spectrum reaches its stationary regime (added by hand), and stores the
resulting run variables as pickle files.

It can be used to generate the plots and results of the numerical simulations
of A. Roper Pol, S. Mandal, A. Brandenburg, T. Kahniashvili, and A. Kosowsky,
"Numerical simulations of gravitational waves from early-universe turbulence,"
Phys. Rev. D 102, 083512 (2020), https://arxiv.org/abs/1903.08585.

Author: Alberto Roper Pol
Created: 01/11/2021
Updated: 05/09/2023 (new release of the cosmoGW code)
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u

# get working directory, where the runs and routines should be stored
dir0 = os.getcwd() + '/'
HOME = dir0 + '/..'

# create directory to store plots
try:
    os.mkdir('plots')
except:
    aa_tst = 1

os.chdir(HOME)
import plot_sets
import run as r
import interferometry as inte
import cosmoGW
os.chdir(dir0)

# reference values and constants
Tref = 100*u.GeV    # EWPT
gref = 100          # EWPT
SNR_PLS = 10        # LISA PLS 
T_PLS = 4           # LISA PLS

def correct_ini1_3(runs):

    """
    Function that corrects the wrong normalization factor in runs ini1, ini2,
    and ini3.
    """

    rrs = ['ini1', 'ini2', 'ini3']
    print('\n')
    print('CORRECTING NORMALIZATION OF EGW IN RUNS ini1, ini2, AND ini3',
          'BY A FACTOR 32 pi/6 \n')
    for i in rrs:
        run = runs.get(i)
        # check if this factor has already been corrected
        if 'corr' not in run.ts:
            EEGW = run.ts.get('EEGW')
            EEGW *= 16*np.pi/3
            run.ts.update({'EEGW': EEGW})
            run.ts.update({'corr': True})
        else:
            print('\n')
            print('NORMALIZATION OF EGW IN RUN ', i, ' ',
                  'HAS ALREADY BEEN CORRECTED')
            
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
    ts = np.zeros(len(indsst))
    for i in runs:
        run = runs.get(i)
        t = run.spectra.get('t_GWs')
        ts[j] = t[indsst[j]]
        indt = indsst[j]
        j += 1
        run.min_max_stat(indt=indt, sp='EGW')
        GWs_stat_sp = run.spectra.get('EGW_stat_sp')
        k = run.spectra.get('k')
        run.OmGWsat = np.trapz(GWs_stat_sp, k)
        
    print('   times:', [s for s in ts], ' \n')
    
def plot_EGW_EM_vs_k(runs, rr='ini2', save=True, show=True):

    """
    Function that generates the plot of the magnetic spectrum
    EM (k) = Omega_M(k)/k at the initial time of turbulence generation
    and the GW spectrum EGW (k) = Omega_GW(k)/k, averaged over oscillations
    in time.

    It corresponds to figure 1 of A. Roper Pol, S. Mandal, A. Brandenburg,
    T. Kahniashvili, and A. Kosowsky, "Numerical simulations of gravitational
    waves from early-universe turbulence," Phys. Rev. D 102, 083512 (2020),
    https://arxiv.org/abs/1903.08585.

    Arguments:
        runs -- dictionary that includes the run variables
        rr -- string that selects which run to plot (default 'ini2')
        save -- option to save the resulting figure as
                plots/EGW_EM_vs_k_'name_run'.pdf' (default True)
        show -- option to show the resulting figure (default True)
    """

    plt.figure(figsize=(10,6))
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(80, 7e4)
    plt.ylim(1e-19, 1e-4)
    plt.xlabel('$k$')
    plt.ylabel(r'$\Omega_{\rm GW}(k)$ and $\Omega_{\rm M}(k)$',
               fontsize=20)

    run = runs.get(rr)
    # plot the averaged over times GW spectrum
    GWs_stat_sp = run.spectra.get('EGW_stat_sp')
    k = run.spectra.get('k')
    plt.plot(k[1:], k[1:]*GWs_stat_sp[1:], color='black')
    # plot magnetic spectrum at the initial time
    mag = run.spectra.get('mag')
    plt.plot(k[1:], k[1:]*mag[0, 1:], color='black')

    # plot k^4 line
    k0 = np.logspace(np.log10(150), np.log10(500), 5)
    plt.plot(k0, k0*1e-9*(k0/100)**4, color='black', ls='-.', lw=.7)
    plt.text(300, 8e-9*300, r'$\sim\!k^5$', fontsize=20)

    # plot k^(-5/3) line
    k0 = np.logspace(np.log10(2000), np.log10(8000), 5)
    plt.plot(k0, k0*1e-5*(k0/1000)**(-5/3)/10, color='black', ls='-.', lw=.7)
    plt.text(5e3, 1.6e-6*5e3/1000, r'$\sim\!k^{-2/3}$', fontsize=20)

    # plot k^(-11/3) line
    k0 = np.logspace(np.log10(3000), np.log10(30000), 5)
    plt.plot(k0, k0*1e-12*(k0/1000)**(-11/3), color='black', ls='-.', lw=.7)
    plt.text(1e4, 5e-16*1e4, r'$\sim\!k^{-8/3}$', fontsize=20)
    
    # plot k line
    k0 = np.logspace(np.log10(300), 3, 5)
    plt.plot(k0, k0*3e-13*(k0/1000)**(1), color='black', ls='-.', lw=.7)
    plt.text(5e2, 3e-14*1e4, r'$\sim\!k$', fontsize=20)

    plt.text(1500, 1e-16*1500, r'$\Omega_{\rm GW} (k)$', fontsize=20)
    plt.text(800, 5e-8*800, r'$\Omega_{\rm M} (k)$', fontsize=20)

    ax = plt.gca()
    ax.set_xticks([100, 1000, 10000])
    ytics = 10**np.array(np.linspace(-19, -5, 7))
    ytics2 = 10**np.array(np.linspace(-19, -5, 15))
    ytics2 = 10**np.array(np.linspace(-15, -1, 8))
    yticss = ['$10^{-19}$', '', '$10^{-17}$', '', '$10^{-15}$', '',
              '$10^{-13}$', '', '$10^{-11}$', '', '$10^{-9}$', '',
              '$10^{-7}$', '', '$10^{-5}$']
    yticss = ['$10^{-15}$',
              '$10^{-13}$', '$10^{-11}$', '$10^{-9}$', 
              '$10^{-7}$', '$10^{-5}$', '$10^{-3}$', '$10^{-1}$']
    ax.set_yticks(ytics2)
    ax.set_yticklabels(yticss)
    plot_sets.axes_lines()
    ax.tick_params(pad=10)
    plt.ylim(1e-16, 1e-1)

    if save: plt.savefig('plots/EGW_EM_vs_k_' + run.name_run + '.pdf',
                         bbox_inches='tight')
    if not show: plt.close()
    
def plot_EGW_vs_kt(runs, rr='ini2', save=True, show=True, k0=100):

    """
    Function that generates the plot of the compensated GW spectrum as a
    function of k(t - tini) for the smallest wave numbers of the run.

    It corresponds to figure 3 of A. Roper Pol, S. Mandal, A. Brandenburg,
    T. Kahniashvili, and A. Kosowsky, "Numerical simulations of gravitational
    waves from early-universe turbulence," Phys. Rev. D 102, 083512 (2020),
    https://arxiv.org/abs/1903.08585.

    Arguments:
        runs -- dictionary that includes the run variables
        rr -- string that selects which run to plot (default 'ini2')
        save -- option to save the resulting figure as
                plots/EGW_vs_kt_'name_run'.pdf (default True)
        show -- option to show the resulting figure (default True)
        k0 -- lowest wave number (used to show 10 instead of 9 in the plot)
    """

    run = runs.get(rr)
    k = run.spectra.get('k')
    EGW = np.array(run.spectra.get('EGW'), dtype='float')
    t = run.spectra.get('t_EGW')

    plt.figure(figsize=(8,5))
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(5e-3, 50)
    plt.ylim(3e-7, 3e-4)
    plt.xlabel('$k (t - 1)$')
    plt.ylabel(r'$\left[k_* \Omega_{\rm GW} (k, t)/k\right]^{1/2}$')
    plot_sets.axes_lines()

    # plot for initial wave numbers
    cols = ['black', 'red', 'blue', 'green']
    lss = ['solid', '-.', 'dashed', 'dotted']
    j = [1, 2, 4, 8]
    for i in range(0, len(j)):
        plt.plot(k[j[i]]*(t - 1), np.sqrt(EGW[:, j[i]]*k0*2*np.pi),
                 color=cols[i - 1], lw=.8, ls=lss[i - 1],
                 label='$k = %i$'%(k0*j[i]))
    plt.legend(fontsize=22, loc='lower right', frameon=False)

    if save: plt.savefig('plots/EGW_vs_kt_' + run.name_run + '.pdf', bbox_inches='tight')
    if not show: plt.close()
    
def plot_OmGW_hc_vs_f_ini(runs, T=Tref, g=gref, SNR=SNR_PLS, Td=T_PLS,
                          save=True, show=True):

    """
    Function that generates the plot of the GW energy density spectrum
    of initial runs (ini1, ini2, and ini3), compared to the LISA sensitivity
    and power law sensitivity (PLS).

    It corresponds to figure 4 of A. Roper Pol, S. Mandal, A. Brandenburg,
    T. Kahniashvili, and A. Kosowsky, "Numerical simulations of gravitational
    waves from early-universe turbulence," Phys. Rev. D 102, 083512 (2020),
    https://arxiv.org/abs/1903.08585.

    Arguments:
        runs -- dictionary that includes the run variables
        T -- temperature scale (in natural units) at the time of turbulence
             generation (default 100 GeV, i.e., electroweak scale)
        g -- number of relativistic degrees of freedom at the time of
             turbulence generation (default 100, i.e., electroweak scale)
        SNR -- signal-to-noise ratio (SNR) of the resulting PLS (default 10)
        Td -- duration of the mission (in years) of the resulting PLS
             (default 4)
        save -- option to save the resulting figure as
                plots/OmGW_vs_f_ini.pdf (default True)
        show -- option to show the resulting figure (default True)
    """

    # read LISA sensitivity
    fs, LISA_Om, LISA_OmPLS = inte.read_sens(SNR=SNR, T=Td, interf='LISA')
    fs = fs*u.Hz

    # chose the runs to be shown
    rrs = ['ini1', 'ini2', 'ini3']
    # chose the colors of each run
    col = ['black', 'red', 'blue']

    plt.figure(1, figsize=(12,5))
    plt.figure(2, figsize=(12,5))

    j = 0
    for i in rrs:
        run = runs.get(i)
        k = run.spectra.get('k')[1:]
        EGW_stat = run.spectra.get('EGW_stat_sp')[1:]
        f, OmGW_stat = cosmoGW.shift_OmGW_today(k, EGW_stat*k, T=T, g=g)
        OmGW_stat = np.array(OmGW_stat, dtype='float')
        hc_stat = cosmoGW.hc_OmGW(f, OmGW_stat, h0=1.)
        plt.figure(1)
        plt.plot(f, OmGW_stat, color=col[j], lw=.8)
        if i == 'ini1': plt.text(5e-2, 1.5e-15, i, color=col[j])
        if i == 'ini2': plt.text(3e-2, 2e-17, i, color=col[j])
        if i == 'ini3': plt.text(3e-3, 4e-17, i, color=col[j])
        plt.figure(2)
        plt.plot(f, hc_stat, color=col[j], lw=.8)
        if i == 'ini1': plt.text(5e-2, 1.5e-24, i, color=col[j])
        if i == 'ini2': plt.text(1e-2, 3e-24, i, color=col[j])
        if i == 'ini3': plt.text(4e-3, 1e-24, i, color=col[j])
        j += 1

    plt.figure(1)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(1e-4, 1e-1)
    plt.ylim(1e-19, 1e-9)
    plt.xlabel('$f$ [Hz]')
    plt.ylabel(r'$h_0^2 \Omega_{\rm GW} (f)$')
    plt.plot(fs, LISA_OmPLS, color='lime', ls='dashdot')
    plt.plot(fs, LISA_Om, color='lime')

    # plot f^(-8/3) line
    fs0 = np.logspace(-2.1, -1.5, 5)
    plt.plot(fs0, 3e-15*(fs0/2e-2)**(-8/3), color='black',
             ls='dashdot', lw=.7)
    plt.text(1e-2, 5e-16, '$\sim\!f^{-8/3}$')

    # plot f line
    fs0 = np.logspace(-3.45, -2.8, 5)
    plt.plot(fs0, 2e-13*(fs0/1e-3)**(1), color='black',
             ls='dashdot', lw=.7)
    plt.text(4e-4, 3e-13, '$\sim\!f$')
    plt.yticks(np.logspace(-18, -10, 5))
    plot_sets.axes_lines()

    if save: plt.savefig('plots/OmGW_vs_f_ini.pdf', bbox_inches='tight')
    if not show: plt.close()

    plt.figure(2)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(1e-4, 1e-1)
    plt.ylim(1e-25, 1e-20)
    plt.xlabel('$f$ [Hz]')
    plt.ylabel(r'$h_{\rm c}(f)$')

    LISA_hc_PLS = cosmoGW.hc_OmGW(fs, LISA_OmPLS)
    LISA_hc = cosmoGW.hc_OmGW(fs, LISA_Om)
    plt.plot(fs, LISA_hc_PLS, color='lime', ls='dashdot')
    plt.plot(fs, LISA_hc, color='lime')

    # plot f^(-1/2) line
    fs0 = np.logspace(-3.4, -2.6, 5)
    plt.plot(fs0, 8e-22*(fs0/1e-3)**(-1/2), color='black',
             ls='dashdot', lw=.7)
    plt.text(8e-4, 1.5e-21, '$\sim\!f^{-1/2}$')

    # plot f^(-7/3) line
    fs0 = np.logspace(-2, -1.4, 5)
    plt.plot(fs0, 1e-23*(fs0/2e-2)**(-7/3), color='black',
             ls='dashdot', lw=.7)
    plt.text(2e-2, 2e-23, '$\sim\!f^{-7/3}$')

    ax = plt.gca()
    ytics2 = 10**np.array(np.linspace(-25, -20, 6))
    ax.set_yticks(ytics2)
    plot_sets.axes_lines()

    if save: plt.savefig('plots/hc_vs_f_ini.pdf', bbox_inches='tight')
    if not show: plt.close()


def plot_OmGW_hc_vs_f_driven(runs, T=Tref, g=gref, SNR=SNR_PLS, Td=T_PLS,
                             save=True, show=True):

    """
    Function that generates the plot of the GW energy density spectrum
    of some of the initially driven runs (ac1, hel1, hel2, hel3, and noh1),
    compared to the LISA sensitivity and power law sensitivity (PLS).

    It corresponds to figure 6 of A. Roper Pol, S. Mandal, A. Brandenburg,
    T. Kahniashvili, and A. Kosowsky, "Numerical simulations of gravitational
    waves from early-universe turbulence," Phys. Rev. D 102, 083512 (2020),
    https://arxiv.org/abs/1903.08585.

    Arguments:
        runs -- dictionary that includes the run variables
        T -- temperature scale (in natural units) at the time of turbulence
             generation (default 100 GeV, i.e., electroweak scale)
        g -- number of relativistic degrees of freedom at the time of
             turbulence generation (default 100, i.e., electroweak scale)
        SNR -- signal-to-noise ratio (SNR) of the resulting PLS (default 10)
        Td -- duration of the mission (in years) of the resulting PLS
             (default 4)
        save -- option to save the resulting figure as
                plots/OmGW_vs_f_driven.pdf (default True)
        show -- option to show the resulting figure (default True)
    """

    # read LISA and Taiji sensitivities
    fs, LISA_Om, LISA_OmPLS = inte.read_sens(SNR=SNR, T=Td, interf='LISA')
    fs = fs*u.Hz

    # chose the runs to be shown
    rrs = ['ac1', 'hel1', 'hel2', 'hel3', 'noh1']
    # chose the colors of each run
    col = ['black', 'red', 'blue', 'blue', 'blue']
    # chose the line style for the plots
    ls = ['solid', 'solid', 'solid', 'dotted', 'dashed']

    plt.figure(1, figsize=(12,5))
    plt.figure(2, figsize=(12,5))

    j = 0
    for i in rrs:
        run = runs.get(i)
        k = run.spectra.get('k')
        EGW_stat = run.spectra.get('EGW_stat_sp')
        f, OmGW_stat = cosmoGW.shift_OmGW_today(k, EGW_stat*k, T=T, g=g)
        OmGW_stat = np.array(OmGW_stat, dtype='float')
        hc_stat = cosmoGW.hc_OmGW(f, OmGW_stat)

        plt.figure(1)
        # omit largest frequencies where there is not enough numerical accuracy
        if i == 'hel3':
            OmGW_stat = OmGW_stat[np.where(f.value<3e-2)]
            hc_stat = hc_stat[np.where(f.value<3e-2)]
            f = f[np.where(f.value<3e-2)]
        if i=='hel3' or i=='hel2' or i=='noh1' or i=='hel1':
            plt.plot(f, OmGW_stat, color=col[j],
                     ls=ls[j], lw=.8, label=i)
        else:
            plt.plot(f, OmGW_stat, color=col[j], ls=ls[j], lw=.8)
            plt.text(5e-2, 2e-19, i, fontsize=20, color='black')

        plt.figure(2)
        if i=='hel3' or i=='hel2' or i=='noh1'or i=='hel1':
            plt.plot(f, hc_stat, color=col[j], ls=ls[j], lw=.8, label=i)
        else:
            plt.plot(f, hc_stat, color=col[j], ls=ls[j], lw=.8)
            plt.text(3e-3, 5e-22, i, color=col[j])
        j += 1

    plt.figure(1)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(1e-4, 1e-1)
    plt.ylim(1e-26, 1e-9)
    plt.xlabel('$f$ [Hz]')
    plt.ylabel(r'$h_0^2 \Omega_{\rm GW} (f)$')
    plt.legend(loc='lower left', frameon=False, fontsize=20)
    plt.plot(fs, LISA_OmPLS, color='lime', ls='dashdot')
    plt.plot(fs, LISA_Om, color='lime')

    # plot f^(-5) line
    fk0 = np.logspace(-2.2, -1.6, 5)
    plt.plot(fk0, 1e-14*(fk0/1e-2)**(-5), color='black',
             ls='dashdot', lw=.7)
    plt.text(1.3e-2, 1e-14, '$\sim\!f^{-5}$')

    # plot f line
    fk0 = np.logspace(-3.3, -2.8, 5)
    plt.plot(fk0, 2e-16*(fk0/1e-3)**(1), color='black',
             ls='dashdot', lw=.7)
    plt.text(6e-4, 1e-17, '$\sim\!f$')

    ax = plt.gca()
    ytics2 = 10**np.array(np.linspace(-25, -9, 16))
    yticss = ['$10^{-25}$', '', '$10^{-23}$', '',  '$10^{-21}$', '',
              '$10^{-19}$','', '$10^{-17}$', '','$10^{-15}$','',
              '$10^{-13}$', '','$10^{-11}$','']
    ax.set_yticks(ytics2)
    ax.set_yticklabels(yticss)
    plot_sets.axes_lines()

    if save: plt.savefig('plots/OmGW_vs_f_driven.pdf', bbox_inches='tight')
    if not show: plt.close()

    plt.figure(2)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(1e-4, 1e-1)
    plt.ylim(1e-28, 1e-20)
    plt.xlabel('$f$ [Hz]')
    plt.ylabel(r'$h_{\rm c}(f)$')
    plt.legend(loc='lower left', frameon=False, fontsize=20)

    LISA_hc_PLS = cosmoGW.hc_OmGW(fs, LISA_OmPLS)
    LISA_hc = cosmoGW.hc_OmGW(fs, LISA_Om)
    plt.plot(fs, LISA_hc_PLS, color='lime', ls='dashdot')
    plt.plot(fs, LISA_hc, color='lime')

    # plot f^(-7/2) line
    fk0 = np.logspace(-2.2, -1.6, 5)
    plt.plot(fk0, 1e-23*(fk0/1e-2)**(-7/2), color='black',
             ls='dashdot', lw=.7)
    plt.text(1.3e-2, 1e-23, '$\sim\!f^{-7/2}$')

    # plot f^(-1/2) line
    fk0 = np.logspace(-3.3, -2.8, 5)
    plt.plot(fk0, 2e-23*(fk0/1e-3)**(-1/2), color='black',
             ls='dashdot', lw=.7)
    plt.text(6e-4, 3e-24, '$\sim\!f^{-1/2}$')

    ax = plt.gca()
    ytics2 = 10**np.array(np.linspace(-28, -20, 9))
    ytics2 = 10**np.array(np.linspace(-28, -20, 9))
    yticss = ['', '$10^{-27}$', '',  '$10^{-25}$', '',
              '$10^{-23}$', '', '$10^{-21}$', '']
    ax.set_yticks(ytics2)
    ax.set_yticklabels(yticss)
    plot_sets.axes_lines()

    if save: plt.savefig('plots/hc_vs_f_driven.pdf', bbox_inches='tight')
    if not show: plt.close()
    
def plot_OmMK_OmGW_vs_t(runs, save=True, show=True):

    """
    Function that generates the plots of the total magnetic/kinetic energy
    density as a function of time ('OmM_vs_t.pdf') and the GW energy density
    as a function of time ('OmGW_vs_t.pdf').

    It corresponds to figure 5 of A. Roper Pol, S. Mandal, A. Brandenburg,
    T. Kahniashvili, and A. Kosowsky, "Numerical simulations of gravitational
    waves from early-universe turbulence," Phys. Rev. D 102, 083512 (2020),
    https://arxiv.org/abs/1903.08585.

    Arguments:
        runs -- dictionary that includes the run variables
        save -- option to save the resulting figure as
                plots/OmGW_vs_t.pdf (default True)
        show -- option to show the resulting figure (default True)
    """

    # chose the runs to be shown
    rrs = ['ini1', 'ini2', 'ini3', 'hel1', 'hel2', 'ac1']
    # chose the colors of each run
    col = ['black', 'red', 'blue', 'red', 'blue', 'black']
    # chose the line style for the plots
    ls = ['solid']*6
    ls[3] = 'dashed'
    ls[4] = 'dashed'
    ls[5] = 'dashed'

    plt.figure(1, figsize=(10,6))
    plt.figure(2, figsize=(10,6))

    j = 0
    for i in rrs:
        run = runs.get(i)
        k = run.spectra.get('k')[1:]
        EGW_stat = run.spectra.get('EGW_stat_sp')[1:]
        t = run.ts.get('t')[1:]
        indst = np.argsort(t)
        t = t[indst]
        EEGW = run.ts.get('EEGW')[1:][indst]
        if run.type == 'EEM': EEM = run.ts.get('EEM')[1:][indst]
        if run.type == 'EEK': EEM = run.ts.get('EEK')[1:][indst]

        plt.figure(1)
        plt.plot(t, EEGW, color=col[j], lw=.8, ls=ls[j])
        # text with run name
        if i=='ini1': plt.text(1.02, 5e-8, i, color=col[j])
        if i=='ini2': plt.text(1.07, 5e-11, i, color=col[j])
        if i=='ini3': plt.text(1.2, 6e-9, i, color=col[j])
        if i=='hel1': plt.text(1.15, 2e-9, i, color=col[j])
        if i=='hel2': plt.text(1.12, 7e-10, i, color=col[j])
        if i=='ac1': plt.text(1.2, 1e-7, i, color=col[j])

        plt.figure(2)
        plt.plot(t, EEM, color=col[j], lw=.8, ls=ls[j])
        # text with run name
        if i=='ini1': plt.text(1.01, 8e-2, i, color=col[j])
        if i=='ini2': plt.text(1.12, 3e-3, i, color=col[j])
        if i=='ini3': plt.text(1.01, 9e-3, i, color=col[j])
        if i=='hel1': plt.text(1.15, 1.3e-2, i, color=col[j])
        if i=='hel2': plt.text(1.02, 1e-3, i, color=col[j])
        if i=='ac1': plt.text(1.17, 1.5e-3, i, color=col[j])

        j += 1

    plt.figure(1)
    plt.yscale('log')
    plt.xlabel('$t$')
    plt.xlim(1, 1.25)
    plt.ylim(2e-11, 2e-7)
    plt.ylabel(r'$\Omega_{\rm GW}$')
    plot_sets.axes_lines()

    if save: plt.savefig('plots/OmGW_vs_t.pdf', bbox_inches='tight')
    if not show: plt.close()

    plt.figure(2)
    plt.yscale('log')
    plt.xlim(1, 1.25)
    plt.ylim(5e-4, 2e-1)
    plt.xlabel('$t$')
    plt.ylabel(r'$\Omega_{\rm M, K}$')
    plot_sets.axes_lines()

    if save: plt.savefig('plots/OmM_vs_t.pdf', bbox_inches='tight')
    if not show: plt.close()
    
def plot_OmGW_vs_OmMK(runs, save=True, show=True):

    """
    Function that generates the plot of the total (saturated) GW energy
    density integrated over wave numbers as a function of the
    magnetic/kinetic energy density.

    It corresponds to figure 7 of A. Roper Pol, S. Mandal, A. Brandenburg,
    T. Kahniashvili, and A. Kosowsky, "Numerical simulations of gravitational
    waves from early-universe turbulence," Phys. Rev. D 102, 083512 (2020),
    https://arxiv.org/abs/1903.08585.

    Arguments:
        runs -- dictionary that includes the run variables
        save -- option to save the resulting figure as
                plots/OmGW_vs_OmMK.pdf (default True)
        show -- option to show the resulting figure (default True)
    """

    # chose the runs to be shown
    rrs = [s for s in runs]
    # chose the colors of each run
    col = ['darkorange', 'darkorange', 'lime', 'red', 'red', 'red',
           'lime', 'red', 'red', 'blue', 'blue', 'blue']

    plt.figure(figsize=(8,5))
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel(r'$\Omega_{\rm M,K}^{\rm max}$')
    plt.ylabel(r'$\Omega_{\rm GW}^{\rm sat}$')
    plt.xlim(1e-3, 2e-1)
    plt.ylim(7e-12, 1e-7)

    j = 0
    for i in runs:
        run = runs.get(i)
        if 'noh' in i:
            plt.scatter(run.Ommax, run.OmGWsat,
                        edgecolors='red', facecolors='none')
        else:
            plt.plot(run.Ommax, run.OmGWsat, 'o', color=col[j])
        j += 1

    OmMs = np.logspace(-2.2, -.8)
    plt.plot(OmMs, 1.7e-6*OmMs**2, color='darkorange')
    OmMs = np.logspace(-2.5, -1.5)
    plt.plot(OmMs, 9e-6*OmMs**2, color='red')
    OmMs = np.logspace(-2.5, -1.7)
    plt.plot(OmMs, 1.5e-5*OmMs**2, color='red', ls='dashed', lw=.7)
    OmMs = np.logspace(-2.8, -1.8)
    plt.plot(OmMs, 3.5e-4*OmMs**2, color='blue')

    plot_sets.axes_lines()
    plt.yticks([1e-11, 1e-10, 1e-9, 1e-8, 1e-7])
    plt.text(1.5e-2, 2e-10, 'ini', color='darkorange')
    plt.text(1.5e-2, 2e-10/2, '$i$=M', color='darkorange')
    plt.text(2e-2, 2e-8, 'hel', color='red')
    plt.text(2e-2, 2e-8/2, '$i$=M', color='red')
    plt.text(2e-3, 1e-8, 'ac', color='blue')
    plt.text(2e-3, 1e-8/2, '$i$=K', color='blue')
    plt.text(7e-3, 4e-9, '(ini3)', color='lime')
    plt.text(5e-3, 1.7e-11, '(hel4)', color='lime')
    plt.text(3.5e-3, 7e-10, '(noh)', color='red')

    if save: plt.savefig('plots/OmGW_vs_OmMK.pdf', bbox_inches='tight')
    if not show: plt.close()
    
def plot_EGW_vs_k_initial_ts(runs, rr='ini2', save=True, show=True):

    """
    Function that generates the plot of the amplitudes of the GW energy
    density at the first time step, to show the evolution from k^2 to k^0
    slopes.

    It corresponds to a plot generated for a presentation, related to the
    work A. Roper Pol, S. Mandal, A. Brandenburg, T. Kahniashvili,
    and A. Kosowsky, "Numerical simulations of gravitational waves from
    early-universe turbulence," Phys. Rev. D 102, 083512 (2020),
    https://arxiv.org/abs/1903.08585.

    Arguments:
        runs -- dictionary that includes the run variables
        rr -- string that selects which run to plot (default 'ini2')
        save -- option to save the resulting figure as
                plots/EGW_vs_k_initial_ts.pdf (default True)
        show -- option to show the resulting figure (default True)
    """
    
    import spectra as sp

    run = runs.get(rr)
    EGW = run.spectra.get('EGW')[:, 1:]
    k = run.spectra.get('k')[1:]
    t = run.spectra.get('t_GWs')
    EGW_av = run.spectra.get('EGW_integ')
    ind_max = np.argmax(EGW_av)
    tmax = t[ind_max]
    max_sp_EGW = EGW[np.argmin(abs(t - tmax)), :]

    plt.figure(figsize=(10,6))
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('$k$')
    plt.ylabel(r'$\Omega_{\rm GW} (k)/k$')
    # initial time step GW spectrum
    plt.plot(k, EGW[0, :], color='black', lw=.7, ls='dashed')
    # local maximum amplitudes of oscillations at next time step
    plt.plot(k, max_sp_EGW, color='grey', lw=.7, ls='dashed', alpha=.7)
    kmax, EGWmax = sp.local_max(k, EGW[1, :])
    plt.plot(kmax, EGWmax, color='black', lw=.7, ls='dashed')
    plt.plot(k[:10], EGW[1, :10], color='black', lw=.7, ls='dashed')
    # maximum over all times of the GW energy density (amplitude of the
    # oscillations)
    kmax, max_sp_EGW2 = sp.local_max(k, max_sp_EGW)
    kks = np.append(k[np.where(k < 8e2)], kmax[np.where(kmax > 8e2)])
    EGW_kks = np.append(max_sp_EGW[np.where(k < 8e2)], max_sp_EGW2[np.where(kmax > 8e2)])
    plt.plot(kks, EGW_kks, color='black', lw=1.5)
    # 6th time step
    kmax, EGWmax = sp.local_max(k, EGW[5, :])
    plt.plot(kmax, EGWmax, color='black', lw=.7, ls='dashed')
    plt.plot(k[:4], EGW[5, :4], color='black', lw=.7, ls='dashed')

    # line k^(2.5)
    k0s = np.linspace(250, 1000)
    plt.plot(k0s, 3e-19*(k0s/100)**2.5, color='black',
             ls='dashdot', lw=.7)
    plt.text(500, 5e-18, '$k^{2.5}$', fontsize=20)

    # line k^(1/2)
    k0s = np.linspace(200, 800)
    plt.plot(k0s, 7e-14*(k0s/100)**0.5, color='black',
             ls='dashdot', lw=.7)
    plt.text(500, 3e-13, '$k^{0.5}$', fontsize=20)

    # line k^(-11/3)
    k0s = np.linspace(3000, 10000)
    plt.text(5000, 1e-14, '$k^{-11/3}$', fontsize=20)
    plt.plot(k0s, 1e-14*(k0s/4e3)**(-11/3), color='black',
             ls='dashdot', lw=.7)

    plt.text(3000, 1e-17, r'$t-1=4 \times 10^{-5}$', fontsize=14)
    plt.text(250, 2e-15, '$t-1=10^{-3}$', fontsize=14)
    plt.text(130, 8e-15, r'$5 \times 10^{-3}$', fontsize=14)

    plt.ylim(1e-19, 2e-12)
    yticks = np.logspace(-19, -12, 8)
    ax = plt.gca()
    ax.set_yticks(yticks)
    plot_sets.axes_lines()
    ax.tick_params(axis='x', pad=12)

    if save: plt.savefig('plots/EGW_vs_k_initial_ts.pdf',
                         bbox_inches='tight')
    if not show: plt.close()

def plot_efficiency(runs, save=True, show=True, sqrt=False):

    """
    Function that generates the plot of the total GW energy density
    compensated by EM^2/kf^2 (efficiency).

    It corresponds to the runs presented in A. Roper Pol, S. Mandal,
    A. Brandenburg, T. Kahniashvili, and A. Kosowsky, "Numerical simulations of gravitational waves from
    early-universe turbulence," Phys. Rev. D 102, 083512 (2020),
    https://arxiv.org/abs/1903.08585.

    This plot is analogous to figure 10 of A. Roper Pol, S. Mandal,
    A. Brandenburg, and T. Kahniashvili, "Polarization of gravitational waves
    from helical MHD turbulent sources", https://arxiv.org/abs/2107.05356.

    Arguments:
        runs -- dictionary that includes the run variables
        save -- option to save the resulting figure as plots/efficiency.pdf
                or plots/efficiency_sqrt.pdf (if sqrt = True) (default True)
        show -- option to show the resulting figure (default True)
        sqrt -- option to plot the efficiency defined as the square root of
                the compensated GW energy density (default False)
    """

    exp = 1
    if sqrt: exp = 1/2
    plt.figure(1, figsize=(12,8))
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(6e-3, 1)
    plt.ylim(1e-1**exp, 2e2**exp)
    plt.xlabel('$\delta t = t - 1$')
    plt.ylabel(r'$q^2 (t) = k_*^2\, \Omega_{\rm GW}$' + \
               r'$(t)/({\cal E}_{\rm M, K}^{\rm max})^2$')

    if sqrt:
        plt.ylabel(r'$q (t) = k_* \left[\Omega_{\rm GW} (t)\right]^{1/2}$' + \
                   r'$/{\cal E}_{\rm M, K}^{\rm max}$')
    for i in runs:
        # chose colors
        if 'ini' in i: col = 'green'
        elif 'hel4' in i: col = 'red'
        elif 'noh' in i: col = 'black'
        elif 'ac' in i: col = 'blue'
        else: col = 'orange'
        run = runs.get(i)
        t = run.ts.get('t')
        EGW = run.ts.get('EEGW')
        ind_max = np.argmax(run.Om_max)
        kf = run.kf_max[ind_max]
        plt.plot(t-1, (EGW/run.Ommax**2*kf**2)**exp,
                 color=col)

    plot_sets.axes_lines()
    plt.text(3e-1, 70**exp, 'acoustic', color='blue')
    plt.text(3e-1, 1.5**exp, 'initial', color='green')
    plt.text(1.3e-1, 2.8**exp, r'helical 1--3', color='orange')
    plt.text(3e-1, 15**exp, 'helical 4', color='red')
    plt.text(4e-1, 8.6**exp, 'non-helical', color='black')

    if save:
        if sqrt: plt.savefig('plots/efficiency_sqrt.pdf',
                         bbox_inches='tight')
        else: plt.savefig('plots/efficiency.pdf',
                         bbox_inches='tight')
    if not show: plt.close()
