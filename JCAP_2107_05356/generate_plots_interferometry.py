"""
generate_plots_interferometry.py is a Python routine
to generate the plots of Appendix B (LISA and Taiji interferometry)
of A. Roper Pol, S. Mandal, A. Brandenburg, and T. Kahniashvili,
"Polarization of gravitational waves from helical MHD turbulent sources",
https://arxiv.org/abs/2107.05356.
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u

plt.rcParams.update({'xtick.labelsize': 'x-large',
                     'ytick.labelsize': 'x-large',
                     'axes.labelsize': 'x-large'})

# get working directory, where the runs and routines should be stored
dir0 = os.getcwd() + '/'
HOME = dir0 + '/..'
os.chdir(HOME)

import plot_sets
import interferometry as inte

# response functions of LISA and Taiji
fs, MAs, MAEs, MTs, DAEs, MAs_Tai, MAEs_Tai, MTs_Tai, DAEs_Tai = \
        inte.read_response_LISA_Taiji()

# Noise power spectral density of LISA channel A and Taiji channel C
PnA, PnT = inte.Pn_TDI(fs)
PnX, PnXY = inte.Pn_f(fs)
L = 3e6*u.km
PnC, PnS = inte.Pn_TDI(fs, P=8, A=3, L=L)
PnX_Tai, PnXY_Tai = inte.Pn_f(fs, P=8, A=3, L=L)
# strain sensitivities
SnA = PnA/MAs
SnT = PnT/MTs
SnX = PnX/MAs
SnC = PnC/MAs_Tai
SnX_Tai = PnX_Tai/MAs_Tai
SnS = PnS/MTs_Tai

# Dipole response sensitivities
v = 1.23e-3
SnAE_Xi = PnA/v/abs(DAEs)
SnCD_Xi = PnC/v/abs(DAEs_Tai)

f_AC, M_AC = inte.read_MAC()
f_AD, M_AD = inte.read_MAC(M='MAD')
f_EC, M_EC = inte.read_MAC(M='MEC')
f_ED, M_ED = inte.read_MAC(M='MED')
f_ED_I, M_ED_I = inte.read_MAC(M='MED', V='I')

# refine and interpolate in fs
f_ED_I, M_ED_I = inte.refine_M(f_ED_I, M_ED_I, A=0.028277196782809974)
M_ED_I = np.interp(fs, f_ED_I, M_ED_I)
M_ED_I[np.where(fs > f_ED_I[-1])] = 1e-50
f_AC, M_AC = inte.refine_M(f_AC, M_AC, A=15, exp=1)
M_AC = np.interp(fs, f_AC, M_AC)
M_AC[np.where(fs > f_AC[-1])] = 1e-50
f_AD, M_AD = inte.refine_M(f_AD, M_AD, A=10, exp=1)
M_AD = np.interp(fs, f_AD, M_AD)
M_AD[np.where(fs > f_AD[-1])] = 1e-50
f_EC, M_EC = inte.refine_M(f_EC[2:], M_EC[2:], A=-9, exp=1)
M_EC = np.interp(fs, f_EC, M_EC)
M_EC[np.where(fs > f_EC[-1])] = 1e-50
f_ED, M_ED = inte.refine_M(f_ED[2:], M_ED[2:], A=-13, exp=1)
M_ED = np.interp(fs, f_ED, M_ED)
M_ED[np.where(fs > f_ED[-1])] = 1e-50

# multiply R functions by 2 to get monopole response function
M_AC*=2.
M_AD*=2.
M_EC*=2.
M_ED*=2.
M_ED_I*=2.

# Sensitivity of the cross-correlated channels of the LISA-Taiji network
# Stokes parameter I
Sn_ED_I = np.sqrt(PnA*PnC)/abs(M_ED_I)
Sn_AC_V = np.sqrt(PnA*PnC)/abs(M_AC)
Sn_AD_V = np.sqrt(PnA*PnC)/abs(M_AD)
Sn_EC_V = np.sqrt(PnA*PnC)/abs(M_EC)
Sn_ED_V = np.sqrt(PnA*PnC)/abs(M_ED)

# noise functions
Pacc = inte.Pacc_f(fs)
Pacc_Tai = inte.Pacc_f(fs, A=3, L=L)
Poms = inte.Poms_f(fs)
Poms_Tai = inte.Poms_f(fs, P=8, L=L)

# GW energy density sensitivity
OmSA = inte.Oms(fs, SnA)
OmST = inte.Oms(fs, SnT)
OmSC = inte.Oms(fs, SnC)
OmSS = inte.Oms(fs, SnS)
OmS_ED_I = .5*inte.Oms(fs, Sn_ED_I)
OmS_comb = 1/np.sqrt(1/OmSA**2 + 1/OmSC**2)

###### GW polarization sensitivity
# From LISA and Taiji dipole response functions
XiSAE = .5*inte.Oms(fs, SnAE_Xi)
XiSCD = .5*inte.Oms(fs, SnCD_Xi)
# From monopole response functions of cross-correlated channels of the
# LISA-Taiji network
XiSAC = inte.Oms(fs, Sn_AC_V)
XiSAD = inte.Oms(fs, Sn_AD_V)
XiSEC = inte.Oms(fs, Sn_EC_V)
XiSED = inte.Oms(fs, Sn_ED_V)
XiS_comb = 1/np.sqrt(1/XiSAD**2 + 1/XiSAC**2 + \
                     1/XiSEC**2 + 1/XiSED**2)
T = 1*u.yr
T = T.to(u.s)
Xiflat = .5/np.sqrt(np.trapz(1/XiSAE**2, fs.value)*T.value)
Xiflat_Tai = .5/np.sqrt(np.trapz(1/XiSCD**2, fs.value)*T.value)
Xiflat_comb = .5/np.sqrt(np.trapz(1/XiS_comb**2, fs.value)*T.value)

os.chdir(dir0)

def get_response():
    return fs, MAs, MAEs, MTs, DAEs, MAs_Tai, MAEs_Tai, MTs_Tai, DAEs_Tai
def get_response_LT():
    return fs, M_AC, M_AD, M_EC, M_ED
def get_sensitivity():
    return fs, SnX, SnA, SnT, SnX_Tai, SnC, SnS
def get_ED_I():
    return fs, M_ED_I, Sn_ED_I, OmS_ED_I
def get_sensitivity_V():
    return fs, Sn_AC_V, Sn_AD_V, Sn_EC_V, Sn_ED_V
def get_noise():
    return fs, Pacc, Poms, PnX, PnA, PnT, Pacc_Tai, Poms_Tai, PnX_Tai, PnC, PnS
def get_Omega_sensitivity():
    return fs, OmSA, OmST, OmSC, OmSS, OmS_comb
def get_Xi_sensitivity():
    return fs, XiSAE, XiSCD, XiSAD, XiSAC, XiSEC, XiSED, XiS_comb
def get_Xiflat():
    return Xiflat, Xiflat_Tai, Xiflat_comb

def get_Omega_PLS(beta):

    """
    Function that returns the PLS of the GW energy density
    for SNR = 1 and 1 year.

    Arguments:
        beta -- array of values of beta to compute the PLS
    """

    OmPLS = inte.OmPLS(fs, OmSA, beta)
    OmPLS_Tai = inte.OmPLS(fs, OmSC, beta)
    OmPLS_comb = inte.OmPLS(fs, OmS_comb, beta)
    return fs, OmPLS, OmPLS_Tai, OmPLS_comb

def get_Xi_PLS_dip(beta, Xi=.25):

    """
    Function that returns the PLS helical GW energy density using the dipole
    response function for SNR = 1 and 1 year.

    Arguments:
        beta -- array of values of beta to compute the PLS
        Xi -- factor to be included in the computation of the PLS
              (default is 0.25, which corresponds to the dipole response,
               Xi=0 can be used for PLS obtained using monopole response)
    Returns:
        fs -- frequency array
        XiPLS -- helical PLS from LISA dipole response function
        XiPLS_Tai -- helical PLS from Taiji dipole response function
    """

    XiPLS = inte.OmPLS(fs, XiSAE, beta, Xi=Xi)
    XiPLS_Tai = inte.OmPLS(fs, XiSCD, beta, Xi=Xi)
    return fs, XiPLS, XiPLS_Tai

def get_Xi_PLS(beta):

    """
    Function that returns the PLS helical GW energy density using the monopole
    response function of cross-correlated channels of the LISA-Taiji network
    for SNR = 1 and 1 year.

    Arguments:
        beta -- array of values of beta to compute the PLS
    Returns:
        fs -- frequency array
        XiPLS_comb -- helical PLS from the combination of the cross-correlated
                      channels of the LISA-Taiji network
    """

    XiPLS_comb = inte.OmPLS(fs, XiS_comb, beta)
    return fs, XiPLS_comb

def plot_Mf(save=True, ED_I=False):

    """
    Function that generates the plot of the LISA and Taiji monopole
    response functions.

    It corresponds to left panel of figure 15 of A. Roper Pol, S. Mandal,
    A. Brandenburg, and T. Kahniashvili, "Polarization of gravitational waves
    from helical MHD turbulent sources," submitted to JCAP,
    https://arxiv.org/abs/2107.05356.

    Arguments:
        save -- option to save the figure in "plots/LISA_Taiji_response.pdf"
                (default True)
        ED_I -- option to plot the response function of cross-correlating
                channels E and D of the LISA-Taiji network (default False)
    """

    plt.figure(figsize=(12,10))
    plt.rc('font', size=30)
    plt.plot(fs, MAs, color='black', label=r'${\cal M}_{\rm AA} (f)$')
    plt.plot(fs, MAs_Tai, color='black', ls='-.',
             label=r'${\cal M}_{\rm CC} (f)$')
    plt.plot(fs, MTs, color='blue', label=r'${\cal M}_{\rm TT} (f)$')
    plt.plot(fs, MTs_Tai, color='blue', ls='-.',
             label=r'${\cal M}_{\rm SS} (f)$')
    aaux = ''
    if ED_I == True:
        gg = np.where(abs(M_ED_I) > 1e-48)
        plt.plot(fs[gg], abs(M_ED_I[gg]), color='purple', alpha=.5)
        plt.text(3e-3, 6e-3, r'$|{\cal M}^I_{\rm ED}|$', color='purple')
        aaux = 'ED'
    Rf = inte.R_f(fs)
    L = 3e6*u.km
    Rf_Tai = inte.R_f(fs, L=L)
    plt.plot(fs, Rf, color='black', ls='dashed', lw=.7)
    plt.plot(fs, Rf_Tai, color='black', ls='dashed', lw=.7)
    plt.text(5e-2, 7e-2, r'$\tilde {\cal R}^{\rm A, C}(f)$', fontsize=34)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(1e-3, 1)
    plt.ylim(1e-7, 1e0)
    plt.legend(fontsize=28, frameon=False, loc='center left')
    plt.xlabel(r'$f$ [Hz]')
    plt.ylabel(r'${\cal M} (f)$')
    plot_sets.axes_lines()
    plt.yticks(np.logspace(-7, 0, 8))

    ax = plt.gca()
    ax.tick_params(axis='x', pad=20)
    ax.tick_params(axis='y', pad=10)

    if save: plt.savefig('plots/LISA_Taiji_response' + aaux + '.pdf',
                         bbox_inches='tight')

def DAE_posneg(f, D, lw=1, ls='solid'):

    """
    Function to plot negative values in red and positive values in blue
    of a function D.

    Arguments:
        f -- frequency array
        D -- dipole response function
        lw -- line width of lines (default 1)
        ls -- line style of lines (default 'solid')
    """

    sgn = np.sign(D)
    converge = False
    sgn0 = sgn[0]
    i = 0
    while not converge:
        sign = False
        i0 = i
        while not sign and not converge:
            if sgn0 == 1: col='blue'
            else: col='red'
            if i==len(sgn) - 2: converge=True
            if sgn[i] != sgn0:
                sign = True
                sgn0 = sgn[i]
                i += 1
            else: i += 1
        plt.plot(f[i0:i + 1], abs(D[i0:i + 1]),
                 color=col, ls=ls, lw=lw)
        if i == len(sgn) - 2: converge=True

def plot_Df(save=True):

    """
    Function that generates the plot of the LISA and Taiji dipole
    response functions.

    It corresponds to right panel of figure 15 of A. Roper Pol, S. Mandal,
    A. Brandenburg, and T. Kahniashvili, "Polarization of gravitational waves
    from helical MHD turbulent sources," submitted to JCAP,
    https://arxiv.org/abs/2107.05356.

    Arguments:
        save -- option to save the figure in "plots/DEA_LISA_Taiji.pdf"
                (default True)
    """

    fig, ax0 = plt.subplots(figsize=(12,10))
    # plt.rc('font', size=30)
    DAE_posneg(fs, DAEs)
    DAE_posneg(fs, DAEs_Tai, ls='-.')
    plt.xlim(1e-3, 1)
    plt.ylim(1e-7, 1)
    plt.xlabel(r'$f$ [Hz]')
    plt.ylabel(r'${\cal D}(f)$')
    plt.xscale('log')
    plt.yscale('log')
    plot_sets.axes_lines()
    plt.text(4.7e-2, 2e-2, r'${\cal D}_{AE}(f)$', fontsize=30)
    plt.text(1.3e-2, 7e-3, r'${\cal D}_{CD}(f)$', fontsize=30)
    #plt.text(2.1e-2, 3e-3, r'${\cal D}_{CD}(f)$', fontsize=20)

    line_pos, = ax0.plot([], [], color='blue', lw=.7,
                           label=r'positive values')
    line_neg, = ax0.plot([], [], color='red', lw=.7,
                           label=r'negative values')
    handles = [line_pos, line_neg]
    lgd = ax0.legend(handles=handles, loc='lower left',
                     fontsize=34, frameon=False)

    ax = plt.gca()
    ax.tick_params(axis='x', pad=20)
    ax.tick_params(axis='y', pad=10)

    plt.yticks(np.logspace(-7, 0, 8))

    if save: plt.savefig('plots/DEA_LISA_Taiji.pdf', bbox_inches='tight')

def plot_sensitivity(save=True):

    """
    Function that generates the plot of LISA and Taiji strain sensitivities.

    Arguments:
        save -- option to save the figure in "plots/sensitivity_LISA_Taiji.pdf"
                (default True)
    """

    plt.figure(figsize=(12,8))
    plt.rc('font', size=20)
    plt.plot(fs, np.sqrt(SnX), color='red', label='channel X')
    plt.plot(fs, np.sqrt(SnA), color='blue', label='TDI channel A')
    plt.plot(fs, np.sqrt(SnT), color='orange', label='TDI channel T',
             alpha=.6)

    plt.plot(fs, np.sqrt(SnX_Tai), color='red', ls='-.')
    plt.plot(fs, np.sqrt(SnC), color='blue', ls='-.')
    plt.plot(fs, np.sqrt(SnS), color='orange', ls='-.', alpha=.6)
    gg = np.where(abs(M_ED_I) > 1e-48)
    plt.plot(fs[gg], np.sqrt(Sn_ED_I[gg]), color='purple', alpha=.5,
             label='TDI channels ED')

    plt.legend(fontsize=20)
    plt.xlabel('$f$ [Hz]')
    plt.ylabel(r'Strain sensitivity $\left[{\rm Hz}^{-1/2}\right]$')
    plt.xscale('log')
    plt.yscale('log')
    plt.text(1e-1, 3e-19, 'LISA', fontsize=20)
    plt.text(1e-1, 8e-21, 'Taiji', fontsize=20)
    plt.ylim(1e-21, 1e-12)
    plt.xlim(1e-5, 1e0)
    plot_sets.axes_lines()
    ax = plt.gca()
    ytics = 10**np.array(np.linspace(-21, -12, 10))
    ax.set_yticks(ytics)

    if save: plt.savefig('plots/sensitivity_LISA_Taiji.pdf',
                         bbox_inches='tight')

def plot_Omega_sensitivity(OmPLS, OmPLS_Tai, OmPLS_comb,
                           SNR=10, T=4, save=True):

    """
    Function that generates the plot of LISA and Taiji GW energy density
    sensitivities and PLS.

    It corresponds to left panel of figure 16 of A. Roper Pol, S. Mandal,
    A. Brandenburg, and T. Kahniashvili, "Polarization of gravitational waves
    from helical MHD turbulent sources," submitted to JCAP,
    https://arxiv.org/abs/2107.05356.

    Arguments:
        OmPLS -- power law sensitivity (PLS) of LISA for T = 1yr and SNR = 1
        OmPLS_Tai -- power law sensitivity (PLS) of Taiji for T = 1yr
                     and SNR = 1
        OmPLS_comb -- power law sensitivity (PLS) of the LISA-Taiji network
                      for T = 1yr and SNR = 1
        SNR -- signal-to-noise ratio (SNR) for the plotted PLS (default 10)
        T -- duration of the observations for the plotted PLS in
             years (default 4)
        save -- option to save the figure in "plots/Omega_LISA_Taiji.pdf"
                (default True)
    """

    import pandas as pd

    fact = SNR/np.sqrt(T)
    fig, ax0 = plt.subplots(figsize=(12,10))
    line_OmPLS, = ax0.plot([], [], color='green',
                           label=r'$\Omega_{\rm PLS}^{\rm A}$')
    line_OmPLS_Tai, = ax0.plot([], [], color='green', ls='-.',
                               label=r'$\Omega_{\rm PLS}^{\rm C}$')
    line_OmPLS_comb, = ax0.plot([], [], color='green', lw=.8,
                               label=r'$\Omega_{\rm PLS}^{\rm comb}$')

    line_OmSA, = ax0.plot([], [], color='blue',
                          label=r'$\Omega_{\rm s}^{\rm A}$')
    line_OmSC, = ax0.plot([], [], color='blue', ls='-.',
                          label=r'$\Omega_{\rm s}^{\rm C}$')
    line_OmS_comb, = ax0.plot([], [], color='blue', lw=.8,
                              label=r'$\Omega_{\rm s}^{\rm comb}$')

    line_OmST, = ax0.plot([], [], color='orange', alpha=.5,
                          label=r'$\Omega_{\rm s}^{\rm T}$')
    line_OmSS, = ax0.plot([], [], color='orange', ls='-.', alpha=.5,
                          label=r'$\Omega_{\rm s}^{\rm S}$')
    line_OmS_ED_I, = ax0.plot([], [], color='purple', alpha=.5,
                              label=r'$\Omega_{\rm s}^{\rm ED}$')

    handles = [line_OmSA, line_OmSC, line_OmS_comb]
    lgd1 = ax0.legend(handles=handles, loc='lower right',
                      fontsize=28, frameon=False)
    handles2 = [line_OmST, line_OmSS, line_OmS_ED_I]
    lgd2 = ax0.legend(handles=handles2, loc='upper left',
                      fontsize=28, frameon=False)
    handles3 = [line_OmPLS, line_OmPLS_Tai, line_OmPLS_comb]
    lgd3 = ax0.legend(handles=handles3, loc='lower left',
                      fontsize=28, frameon=False)
    ax0.add_artist(lgd1)
    ax0.add_artist(lgd2)

    # plt.rc('font', size=30)
    plt.plot(fs, OmSA, color='blue')
    plt.plot(fs, OmST, color='orange', alpha=.5)
    plt.plot(fs, OmPLS*fact, color='green')
    plt.plot(fs, OmSC, '-.', color='blue')
    plt.plot(fs, OmSS, '-.', color='orange', alpha=.5)
    plt.plot(fs, OmPLS_Tai*fact, '-.', color='green')
    gg = np.where(abs(M_ED_I) > 1e-48)
    plt.plot(fs[gg], OmS_ED_I[gg], color='purple', alpha=.5)
    plt.plot(fs, OmS_comb, color='blue', lw=.8)
    plt.plot(fs, OmPLS_comb*fact, color='green', lw=.8)

    # Add reference PLS from Caprini et al 2019
    dir = '../detector_sensitivity/'
    df = pd.read_csv(dir + 'OmegaPLS_Caprinietal19.csv')
    f = np.array(df['f'])
    ff = np.logspace(np.log10(f[0]), np.log10(f[-1]), 70)
    OmGW_PLS_LISA = np.array(df['Omega'])
    OmGW_PLS_LISA = np.interp(ff, f, OmGW_PLS_LISA)
    plt.plot(ff, OmGW_PLS_LISA, '.', color='green', alpha=.4)

    plot_sets.axes_lines()
    plt.ylim(1e-14, 1e-2)
    plt.yticks(np.logspace(-14, -2, 8))
    plt.xlim(1e-5, 1e0)
    plt.xlabel('$f$ [Hz]')
    plt.ylabel(r'$h_0^2\,\Omega_{\rm s} (f)$')
    plt.xscale('log')
    plt.yscale('log')
    ax = plt.gca()
    ax.tick_params(axis='x', pad=20)
    ax.tick_params(axis='y', pad=10)
    plt.yticks(np.logspace(-14, -2, 7))
    plt.xticks(np.logspace(-5, 0, 6))

    if save: plt.savefig('plots/Omega_LISA_Taiji.pdf',
                         bbox_inches='tight')

def plot_Xi_sensitivity(XiPLSa, XiPLSb, XiPLSa_Tai, XiPLSb_Tai, XiPLS_comb,
                        SNR=10, T=4, save=True):

    """
    Function that generates the plot of LISA and Taiji GW energy density
    polarization sensitivities and PLS.

    It corresponds to right panel of figure 16 of A. Roper Pol, S. Mandal,
    A. Brandenburg, and T. Kahniashvili, "Polarization of gravitational waves
    from helical MHD turbulent sources," submitted to JCAP,
    https://arxiv.org/abs/2107.05356.

    Arguments:
        XiPLSa -- polarization PLS of the dipole response function of LISA
                  with beta_max = 2 for T = 1yr and SNR = 1
        XiPLSb -- polarization PLS of the dipole response function of LISA
                  with beta_max = 3 for T = 1yr and SNR = 1
        XiPLSa_Tai -- polarization PLS of the dipole response function of Taiji
                      with beta_max = 2 for T = 1yr and SNR = 1
        XiPLSb_Tai -- polarization PLS of the dipole response function of Taiji
                      with beta_max = 3 for T = 1yr and SNR = 1
        XiPLS_comb -- polarization PLS of the LISA-Taiji network
                      for T = 1yr and SNR = 1
        SNR -- signal-to-noise ratio (SNR) for the plotted PLS (default 10)
        T -- duration of the observations for the plotted PLS in
             years (default 4)
        save -- option to save the figure in "plots/Xi_LISA_Taiji.pdf"
                (default True)
    """

    fact = SNR/np.sqrt(T)
    fig, ax0 = plt.subplots(figsize=(12,10))
    #plt.rc('font', size=18)
    plt.plot(fs, XiSAE, color='blue')
    plt.plot(fs, XiSCD, color='blue', ls='-.')

    gg = np.where(XiS_comb < 1e5)
    plt.plot(fs[gg], XiS_comb[gg], color='blue', lw=.6)

    plt.plot(fs, XiPLSa*fact, color='green')
    plt.plot(fs, XiPLSb*fact, color='green')
    plt.plot(fs, XiPLSa_Tai*fact, color='green', ls='-.')
    plt.plot(fs, XiPLSb_Tai*fact, color='green', ls='-.')

    plt.plot(fs, XiPLS_comb*fact, color='green', lw=.6)

    line_XiPLS, = ax0.plot([], [], color='green',
                           label=r'$\Xi_{\rm PLS}^{\rm AE}$')
    line_XiPLS_Tai, = ax0.plot([], [], color='green', ls='-.',
                               label=r'$\Xi_{\rm PLS}^{\rm CD}$')
    line_XiPLS_comb, = ax0.plot([], [], color='green', lw=.8,
                               label=r'$\Xi_{\rm PLS}^{\rm comb}$')

    line_XiSA, = ax0.plot([], [], color='blue',
                          label=r'$\Xi_{\rm s}^{\rm AE}$')
    line_XiSC, = ax0.plot([], [], color='blue', ls='-.',
                          label=r'$\Omega_{\rm s}^{\rm CD}$')
    line_XiS_comb, = ax0.plot([], [], color='blue', lw=.8,
                              label=r'$\Xi_{\rm s}^{\rm comb}$')

    handles = [line_XiSA, line_XiSC, line_XiS_comb]
    lgd1 = ax0.legend(handles=handles, loc='lower right',
                      fontsize=28, frameon=False)
    handles2 = [line_XiPLS, line_XiPLS_Tai, line_XiPLS_comb]
    lgd2 = ax0.legend(handles=handles2, loc='lower left',
                      fontsize=28, frameon=False)
    ax0.add_artist(lgd1)

    #plt.legend(fontsize=30, loc='lower right', frameon=False)
    plt.xlabel(r'$f$ [Hz]')
    plt.ylabel(r'$h_0^2\, \Xi_{\rm s} (f)$')
    plt.xscale('log')
    plt.yscale('log')
    plt.ylim(1e-14, 1e-2)
    plt.xlim(1e-5, 1e0)
    plot_sets.axes_lines()

    plt.text(1e-1, 8e-9, r'$\beta_{\rm max}=2$', color='green', fontsize=30)
    plt.text(9e-2, 6e-4, r'$\beta_{\rm max}=3$', color='green', fontsize=30)

    ax = plt.gca()
    ax.tick_params(axis='x', pad=20)
    ax.tick_params(axis='y', pad=10)
    plt.yticks(np.logspace(-14, -2, 7))
    plt.xticks(np.logspace(-5, 0, 6))

    if save: plt.savefig('plots/Xi_LISA_Taiji.pdf',
                         bbox_inches='tight')

def plot_Xi_sensitivity_dipole(XiPLSa, XiPLSb, XiPLSc, XiPLSd,
                               XiPLSe, XiPLSf, XiPLSg, XiPLS0, XiPLS_comb,
                               interf='LISA', SNR=10, T=4, save=True):

    """
    Function that generates the plot of the GW energy density
    polarization PLS for different values of beta_max.

    It corresponds to figure 17 of A. Roper Pol, S. Mandal,
    A. Brandenburg, and T. Kahniashvili, "Polarization of gravitational waves
    from helical MHD turbulent sources," submitted to JCAP,
    https://arxiv.org/abs/2107.05356.

    Arguments:
        XiPLSa -- polarization PLS of the dipole response function of LISA
                  with beta_max = 2 for T = 1yr and SNR = 1
        XiPLSb -- polarization PLS of the dipole response function of LISA
                  with beta_max = 3 for T = 1yr and SNR = 1
        XiPLSc -- polarization PLS of the dipole response function of LISA
                  with beta_max = 3.7 for T = 1yr and SNR = 1
        XiPLSd -- polarization PLS of the dipole response function of LISA
                  with beta_max = 3.95 for T = 1yr and SNR = 1
        XiPLSe -- polarization PLS of the dipole response function of LISA
                  with beta_max = 3.999 for T = 1yr and SNR = 1
        XiPLSf -- polarization PLS of the dipole response function of LISA
                  with beta_max = 3.999999 for T = 1yr and SNR = 1
        XiPLSg -- polarization PLS of the dipole response function of LISA
                  with beta values in (-20, 0)U(5, 20) for T = 1yr and SNR = 1
        XiPLS0 -- polarization PLS of the dipole response function of LISA
                  omitting 1/(1 - beta/4) term for T = 1yr and SNR = 1
        XiPLS_comb -- polarization PLS of the LISA-Taiji network
                      for T = 1yr and SNR = 1
        interf -- selects interferometer (default 'LISA',
                  also available 'Taiji')
        SNR -- signal-to-noise ratio (SNR) for the plotted PLS (default 10)
        T -- duration of the observations for the plotted PLS in
             years (default 4)
        save -- option to save the figure in "plots/Xi_PLS_interf.pdf"
                (default True)
    """

    import pandas as pd
    fact = SNR/np.sqrt(T)
    plt.figure(figsize=(12,10))
    if interf=='LISA':
        AE = 'AE'
        ybeta_a = 1.3e-8
        ybeta_b = 1.5e-7
        ybeta_c = 3e-8
        ybeta_d = 1e-7
        ybeta_e = 1e-8
        ybeta_f = 3e-8
    if interf=='Taiji':
        AE = 'CD'
        ybeta_a = 3e-9
        ybeta_b = 3.5e-8
        ybeta_c = 8e-9
        ybeta_d = 3e-8
        ybeta_e = 2.5e-9
        ybeta_f = 6e-9
    # plot PLS using dipole for different beta_max
    plt.plot(fs, XiPLSa*fact, color='blue', lw=1,
             label=r'$\Xi_{\rm PLS}^{\rm %s}$'%AE)
    plt.plot(fs, XiPLSb*fact, color='blue', lw=.8, ls='-.')
    plt.plot(fs, XiPLSc*fact, color='blue', lw=.8, ls='-.')
    plt.plot(fs, XiPLSd*fact, color='blue', lw=.8, ls='-.')
    plt.plot(fs, XiPLSe*fact, color='blue', lw=.8, ls='-.')
    plt.plot(fs, XiPLSf*fact, color='blue', lw=.8, ls='-.')
    # plot PLS using dipole ignoring 1/abs(1 - beta/4) term
    plt.plot(fs, XiPLS0*fact, color='red', lw=.8,
             label=r'$\Xi_{\rm PLS}^{0, {\rm %s}}$'%AE)
    # plot PLS using LISA-Taiji combined network
    plt.plot(fs, XiPLS_comb*fact, color='green', lw=.8,
             label=r'$\Xi_{\rm PLS}^{\rm comb}$')

    # text with values of beta_max
    plt.text(2.5e-2, ybeta_a/4, r'$\beta_{\rm max}=2$', color='blue',
             fontsize=22,
             bbox=dict(facecolor='white', edgecolor='none',
                      boxstyle='round,pad=.2'))
    plt.text(3.2e-2, ybeta_b/4, r'$\beta_{\rm max}=3$', color='blue',
             fontsize=22,
             bbox=dict(facecolor='white', edgecolor='none',
                      boxstyle='round,pad=.2'))
    plt.text(5.3e-3, ybeta_c, r'$\beta_{\rm max}=3.7$', color='blue',
             fontsize=22,
             bbox=dict(facecolor='white', edgecolor='none',
                      boxstyle='round,pad=.2'))
    plt.text(5.5e-3, ybeta_d*3, r'$\beta_{\rm max}=3.95$', color='blue',
             fontsize=22,
             bbox=dict(facecolor='white', edgecolor='none',
                      boxstyle='round,pad=.2'))
    plt.text(8e-4, ybeta_e, r'$\beta_{\rm max}=3.999$', color='blue',
            fontsize=22)
    plt.text(2e-4, ybeta_f*5, r'$\beta_{\rm max}=3.999999$',
             fontsize=22, color='blue')

    plt.hlines(Xiflat*fact, 1e-3, 7e-2, color='blue', lw=.7)
    plt.hlines(Xiflat_Tai*fact, 1e-3, 7e-2, color='blue',
               ls='-.', lw=.7)
    plt.hlines(Xiflat_comb*fact, 1.e-3, 7e-2, color='green', lw=.8)
    plt.text(3e-2, 1.5e-10, r'$\Xi_{\rm flat}^{\rm AE}$',
             fontsize=28, color='blue')
    plt.text(3e-2, 1.7e-11, r'$\Xi_{\rm flat}^{\rm CD}$',
             fontsize=28, color='blue')
    plt.text(3e-2, 8e-13, r'$\Xi_{\rm flat}^{\rm comb}$',
             fontsize=28, color='green')

    # read reference Xi PLS from Ellis et al 2019
    if interf=='Taiji':
        dir = '../detector_sensitivity/'
        df = pd.read_csv(dir + 'XiPLS_Ellisetal19.csv')
        f = np.array(df['f'])
        ff = np.logspace(np.log10(f[0]), np.log10(f[-1]), 50)
        XiGW_PLS_LISA = np.array(df['Xi'])
        XiGW_PLS_LISA = np.interp(ff, f, XiGW_PLS_LISA)
        plt.plot(ff, XiGW_PLS_LISA, '.', color='blue', alpha=.4)
        plt.text(1.15e-2, 1.5e-10, r'$\sim\!f^5$',
                 fontsize=22, color='blue')
        gg = np.where(XiPLSg > Xiflat_Tai*1.01)
        hh = np.where(fs.value[gg] > 4e-3)
        plt.plot(fs[gg][hh], XiPLSg[gg][hh]*fact, color='blue', lw=.5)

    plt.legend(fontsize=28, loc='lower left', frameon=False)
    plt.xlabel(r'$f$ [Hz]')
    plt.ylabel(r'$h_0^2\, \Xi_{\rm PLS} (f)$')
    plt.xscale('log')
    plt.yscale('log')
    plt.ylim(1e-13, 1e-6)
    plt.xlim(1e-4, 1e-1)
    plot_sets.axes_lines()

    ax = plt.gca()
    ax.tick_params(axis='x', pad=20)
    ax.tick_params(axis='y', pad=10)
    plt.yticks(np.logspace(-13, -6, 8))

    if save: plt.savefig('plots/Xi_PLS_' + interf + '.pdf',
                         bbox_inches='tight')

def plot_MAC(save=True, log=False):

    """
    Function that generates the plot of the helical (V Stokes parameter)
    monopole responses of the cross-correlated channels of the LISA-Taiji
    network.

    It corresponds to figure 18 of A. Roper Pol, S. Mandal,
    A. Brandenburg, and T. Kahniashvili, "Polarization of gravitational waves
    from helical MHD turbulent sources," submitted to JCAP,
    https://arxiv.org/abs/2107.05356.

    Arguments:
        save -- option to save the figure in
                "plots/Mcross_LISA_Taiji.pdf" (default True)
        log -- option to plot loglog with absolute values of the response
               functions (default False)
    """

    plt.figure(figsize=(12,8))
    #plt.rc('font', size=20)
    lg = ''
    if log:
        MAC = abs(M_AC)
        MAD = abs(M_AD)
        MEC = abs(M_ED)
        MED = abs(M_ED)
        MED_I = abs(M_ED_I)
        plt.xscale('log')
        plt.yscale('log')
        lg = '_log'
    else:
        MAC = M_AC
        MAD = M_AD
        MEC = M_EC
        MED = M_ED
        MED_I = M_ED_I
    plt.plot(fs, MAC, color='black',
             label=r'${\cal M}^V_{\rm AC}$')
    plt.plot(fs, MAD, color='blue', ls='dotted',
             label=r'${\cal M}^V_{\rm AD}$')
    plt.plot(fs, MEC, color='red', ls='-.',
             label=r'${\cal M}^V_{\rm EC}$')
    plt.plot(fs, MED, color='green', ls='--',
             label=r'${\cal M}^V_{\rm ED}$')
    gg = np.where(abs(MED_I) > 1e-48)
    plt.plot(fs[gg], MED_I[gg], color='purple', alpha=.8,
             label=r'${\cal M}^I_{\rm ED}$')
    plt.legend(fontsize=18, frameon=False)
    plt.xlabel('$f$ [Hz]')
    plt.ylabel(r'${\cal M} (f)$')
    plt.xscale('log')
    if log:
        plt.xlim(1e-4, 4e-2)
        plt.ylim(1e-4, 2e-1)
    else:
        plt.xlim(3e-5, 4e-2)
        plt.ylim(-0.08, 0.08)
        plt.yticks(np.linspace(-.075, .075, 7))
    plot_sets.axes_lines()

    if save: plt.savefig('plots/Mcross_LISA_Taiji' + lg + '.pdf',
                         bbox_inches='tight')

def plot_Xi_sensitivity_comb(save=True):

    """
    Function that generates the plot of the GW energy density
    polarization sensitivity obtained by combining the cross-correlated
    channels of the LISA-Taiji network.

    Arguments:
        save -- option to save the figure in "plots/Xi_LISA_Taiji_comb.pdf"
                (default True)
    """

    plt.figure(figsize=(12,8))
    #plt.rc('font', size=18)
    plt.plot(fs, XiSAC, color='black', lw=.8,
             label = r'$h_0^2\, \Xi_{\rm s}^{AC} (f)$')
    plt.plot(fs, XiSAD, color='blue', ls='dotted',
         label = r'$h_0^2\, \Xi_{\rm s}^{AD} (f)$')
    plt.plot(fs, XiSEC, color='red', ls='-.',
             label = r'$h_0^2\, \Xi_{\rm s}^{EC} (f)$')
    plt.plot(fs, XiSED, color='green', ls='--', lw=.6,
             label = r'$h_0^2\, \Xi_{\rm s}^{ED} (f)$')
    plt.plot(fs, XiS_comb, color='blue',
             label = r'$h_0^2\, \Xi_{\rm s}^{\rm comb} (f)$')

    plt.legend(fontsize=20, loc='lower right', frameon=False)
    plt.xlabel(r'$f$ [Hz]')
    plt.ylabel(r'$h_0^2\, \Xi_{\rm s} (f)$')
    plt.xscale('log')
    plt.yscale('log')
    plt.ylim(1e-12, 1e-4)
    plt.xlim(1e-4, 1e-1)
    plot_sets.axes_lines()

    if save: plt.savefig('plots/Xi_LISA_Taiji_comb.pdf',
                         bbox_inches='tight')

def plot_noise_PSD(interf='LISA', save=True):

    """
    Function that generates the plot of LISA (Taiji) noise PSD of the channel
    A (C) and compares with the optical metrology system P_oms and
    mass acceleration P_acc PSD noises.

    Arguments:
        interf -- selects interferometer (default 'LISA',
                  also available 'Taiji')
        save -- option to save the figure in "plots/noise_PSD_interf.pdf"
                (default True)
    """

    plt.figure(figsize=(12,8))
    if interf=='LISA':
        #plt.plot(fs, PnX, color='blue')
        plt.plot(fs, PnA, color='red')
        plt.plot(fs, Pacc, color='blue', ls='-.', lw=.8)
        plt.plot(fs, Poms, color='blue', ls='-.', lw=.8)
        yPn = 1e-40
        yPoms = 1e-40
        yPacc = 6e-43
        A = 'A'

    if interf=='Taiji':
        #plt.plot(fs, PnX_Tai, color='blue')
        plt.plot(fs, PnC, color='red')
        plt.plot(fs, Pacc_Tai, color='blue', ls='-.', lw=.8)
        plt.plot(fs, Poms_Tai, color='blue', ls='-.', lw=.8)
        yPn = 2e-41
        yPoms = 2e-41
        yPacc = 4e-43
        A = 'C'

    plt.xlabel('$f$ [Hz]')
    plt.ylabel('noise PSD $(f)$')
    plt.xscale('log')
    plt.yscale('log')
    plt.ylim(1e-43, 1e-30)
    plt.xlim(1e-5, 1e0)
    plt.text(1e-1, yPoms, r'$P_{\rm oms} (f)$', color='blue')
    plt.text(1e-1, yPacc, r'$P_{\rm acc} (f)$', color='blue')
    plt.text(1e-2, yPn, r'$P_n^{\rm %s} (f)$'%A, color='red')
    plt.text(1e-1, 1e-32, interf, fontsize=24,
             bbox=dict(facecolor='none', edgecolor='black',
                       boxstyle='round,pad=.5'))
    plot_sets.axes_lines()

    if save: plt.savefig('plots/noise_PSD_' + interf + '.pdf',
                         bbox_inches='tight')
