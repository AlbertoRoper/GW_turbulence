"""
pta.py is a Python routine with functions used in the analysis of observations
by pulsar timing array (PTA) collaborations: NANOGrav, PPTA, EPTA, and IPTA;
in the context of GW backgrounds produced by MHD turbulence in the early
universe.

The main reference for the study is:
- A. Roper Pol, C. Caprini, A. Neronov, D. Semikoz,
"The gravitational wave signal from primordial magnetic fields in the
Pulsar Timing Array frequency band," Phys. Rev. D 105, 123502 (2022),
arXiv:2201.05630.

Author: Alberto Roper Pol
Date: 01/12/2021
"""

import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt
import cosmoGW

### duration of PTA observations
T_NG = 12.5
T_P = 15
T_E = 24
T_I = 31

def get_gamma_A(file, beta_b=False, Omega_b=False, fref=0, disc=1000,
                plot=False, color='blue', alpha=.4, return_all=False,
                fill=True):

    """
    Function to read files with results of PTA amplitude vs slope of
    GW backgrounds.

    Arguments:
        file -- name of the file to read
        beta_b -- option to return the slope of the GW background
                  instead of gamma (beta = 5 - gamma) (default False)
        Omega_b -- option to return the amplitude of the GW energy density
                   instead of the amplitude of the characteristic strain
                   (default False)
        fref -- frequency used as reference to define the PL of the GW
                background (default is 1/(1 year))
        disc -- number of discretization points for the reconstructed A vs
                gamma function (default 1000)
        plot -- option to generate the plot (default False)
        color -- option to chose the color of the plot (default 'blue')
        alpha -- option to chose the transparency of the shaded region within
                 the range of allowed amplitudes for every slope
                 (default 0.4)
        return_all -- option to return all functions computed

    Returns:
        gamma -- slope of the power spectral density
        A -- amplitude
        if return_all is selected:
            gammas -- refined equidistant gamma values
            A1, A2 -- maximum and minimum allowed values of the amplitude
                      by the observations
    """

    import pandas as pd

    dir = '../detector_sensitivity/PTA/'
    df = pd.read_csv(dir + file)
    gamma = np.array(df['gamma'])
    A = np.array(df['A'])
    gammas, A1, A2 = gammas_As(gamma, A, disc=disc)
    if beta_b: gammas = 5 - gammas
    if Omega_b:
        if beta_b:
            A = cosmoGW.Omega_A(A=A, fref=fref, beta=gamma)
            A1 = cosmoGW.Omega_A(A=A1, fref=fref, beta=gammas)
            A2 = cosmoGW.Omega_A(A=A2, fref=fref, beta=gammas)
        else:
            A = cosmoGW.Omega_A(A=A, fref=fref, beta=5 - gamma)
            A1 = cosmoGW.Omega_A(A=A1, fref=fref, beta=5 - gammas)
            A2 = cosmoGW.Omega_A(A=A2, fref=fref, beta=5 - gammas)
    else:
        if fref != 0:
            if beta_b:
                alphas = gammas/2 - 1
                alpha0 = gamma/2 - 1
            else:
                alphas = .5*(3 - gammas)
                alpha0 = .5*(3 - gamma)
            fyr = 1/u.yr
            fyr = fyr.to(u.Hz)
            A *= (fref.value/fyr.value)**alpha0
            A1 *= (fref.value/fyr.value)**alphas
            A2 *= (fref.value/fyr.value)**alphas

    if plot:
        plt.plot(gammas, A1, lw=2, color=color)
        plt.plot(gammas, A2, lw=2, color=color)
        plt.vlines(gammas[0], A1[0], A2[0], color=color, lw=2)
        if fill: plt.fill_between(gammas, A1, A2, color=color, alpha=alpha)

    if return_all:
        return gamma, A, gammas, A1, A2
    else: return gamma, A

def gammas_As(gamma, A, disc=1000):

    """
    Function that uses the amplitude vs slope gamma data to divide into upper
    bound and lower bounds.
    It assumes that the data starts with the smallest gamma and it is given
    in counterclockwise order in the amplitude vs slope plot.

    Arguments:
        gamma -- array of slopes of the power spectral density
        A -- array of amplitudes of the characteristic strain of the background
        disc -- number of discretization points for the reconstructed A vs
                gamma function (default 1000)
    """

    # min and max values of gamma
    gam_0 = gamma[0]
    ind = np.argmax(gamma)
    gam_f = gamma[ind]
    # subdivide array of gammas and
    gamma1 = gamma[:ind]
    gamma2 = gamma[ind:]
    A1 = A[:ind]
    A2 = A[ind:]
    # reorder upper limit
    inds = np.argsort(gamma2)
    gamma2 = gamma2[inds]
    A2 = A2[inds]
    # discretize and interpolate for equidistant arrays in gamma
    gammas = np.linspace(gam_0, gam_f, disc)
    A1 = 10**np.interp(np.log10(gammas), np.log10(gamma1), np.log10(A1))
    A2 = 10**np.interp(np.log10(gammas), np.log10(gamma2), np.log10(A2))

    return gammas, A1, A2

def read_PTA_data(beta_b=True, Omega_b=True, fref=0, return_all=False,
                  plot=False, fill_1s=True, fill_2s=True):

    """
    Function that reads the data from the PTA observations (in terms of region
    of the allowed amplitudes and slopes of the GW background).

    h_c(f) = A (f/fyr)^((3 - gamma)/2)

    Files are stored in detector_sensitivity/PTA, see README for the references
    used and description of the files.

    Arguments:
        beta_b -- option to return the slope of the GW background
                  instead of gamma (beta = 5 - gamma) (default True)
        Omega_b -- option to return the amplitude of the GW energy density
                   instead of the amplitude of the characteristic strain
                   (default True)
        fref -- frequency used as reference to define the PL of the GW
                background (default is 1/(1 year))
        return_all -- option to return all the computed amplitudes and slopes
        plot -- option to generate the plots of the GW background allowed by
                the regions of amplitudes A and slopes gamma

    Returns:
        A -- amplitudes of the GW background
        gamma -- slopes of the GW background
        if return_all is selected:
            gamma -- refined equidistant gamma values
            A1, A2 -- maximum and minimum allowed values of the amplitude
                      by the observations

            NG, P, E, and I indicate the NANOGrav, the PPTA, the EPTA, and the
            IPTA collaborations;
            sPL and bPL indicate single and broken power law fits;
            1s and 2s indicate the 1sigma and 2sigma confidence intervals
    """

    # NANOGrav results for a single power law with 1-sigma confidence
    file = 'NANOGrav_singlePL_1s.csv'
    _ = get_gamma_A(file, beta_b=beta_b, Omega_b=Omega_b,
                    fref=fref, disc=300, return_all=return_all,
                    plot=plot, color='green', alpha=.8, fill=fill_1s)
    if return_all:
        gamma_NG_sPL_1s, A1_NG_sPL_1s, A2_NG_sPL_1s = [_[2], _[3], _[4]]
    else: gamma_NG_sPL_1s, A_NG_sPL_1s = [_[0], _[1]]

    # NANOGrav results for a single power law with 1-sigma confidence
    file = 'NANOGrav_singlePL_2s.csv'
    _ = get_gamma_A(file, beta_b=beta_b, Omega_b=Omega_b,
                    fref=fref,  disc=5000, return_all=return_all,
                    plot=plot, color='green', fill=fill_2s)
    if return_all:
        gamma_NG_sPL_2s, A1_NG_sPL_2s, A2_NG_sPL_2s = [_[2], _[3], _[4]]
    else: gamma_NG_sPL_2s, A_NG_sPL_2s = [_[0], _[1]]

    # NANOGrav results for a broken power law with 1-sigma confidence
    file = 'NANOGrav_brokenPL_1s.csv'
    _ = get_gamma_A(file, beta_b=beta_b, Omega_b=Omega_b,
                    fref=fref, disc=300, return_all=return_all,
                    plot=plot, color='blue', alpha=.8, fill=fill_1s)
    if return_all:
        gamma_NG_bPL_1s, A1_NG_bPL_1s, A2_NG_bPL_1s = [_[2], _[3], _[4]]
    else: gamma_NG_bPL_1s, A_NG_bPL_1s = [_[0], _[1]]

    # NANOGrav results for a broken power law with 1-sigma confidence
    file = 'NANOGrav_brokenPL_2s.csv'
    _ = get_gamma_A(file, beta_b=beta_b, Omega_b=Omega_b,
                    return_all=return_all, fill=fill_2s,
                    fref=fref, plot=plot, color='blue')
    if return_all:
        gamma_NG_bPL_2s, A1_NG_bPL_2s, A2_NG_bPL_2s = [_[2], _[3], _[4]]
    else: gamma_NG_bPL_2s, A_NG_bPL_2s = [_[0], _[1]]

    # PPTA results with 1-sigma confidence
    file = 'PPTA_blue_1s.csv'
    _ = get_gamma_A(file, beta_b=beta_b, Omega_b=Omega_b,
                    return_all=return_all, fill=fill_1s,
                    fref=fref, plot=plot, color='red', alpha=.8)
    if return_all:
        gamma_P_b_1s, A1_P_b_1s, A2_P_b_1s = [_[2], _[3], _[4]]
    else: gamma_P_b_1s, A_P_b_1s = [_[0], _[1]]

    # PPTA results with 2-sigma confidence
    file = 'PPTA_blue_2s.csv'
    _ = get_gamma_A(file, beta_b=beta_b, Omega_b=Omega_b,
                    return_all=return_all, fill=fill_2s,
                    fref=fref, plot=plot, color='red')
    if return_all:
        gamma_P_b_2s, A1_P_b_2s, A2_P_b_2s = [_[2], _[3], _[4]]
    else: gamma_P_b_2s, A_P_b_2s = [_[0], _[1]]

    # EPTA results with 1-sigma confidence from old reference
    file = 'EPTA_singlePL_1s_old.csv'
    _ = get_gamma_A(file, beta_b=beta_b, Omega_b=Omega_b,
                    return_all=return_all, fill=fill_1s,
                    fref=fref, plot=False, color='purple', alpha=.8)
    if return_all:
        gamma_E_1s_old, A1_E_1s_old, A2_E_1s_old = [_[2], _[3], _[4]]
    else: gamma_E_1s_old, A_E_1s_old = [_[0], _[1]]

    # EPTA results with 2-sigma confidence from old reference
    file = 'EPTA_singlePL_2s_old.csv'
    _ = get_gamma_A(file, beta_b=beta_b, Omega_b=Omega_b,
                    return_all=return_all, fill=fill_2s,
                    fref=fref, plot=False, color='purple')
    if return_all:
        gamma_E_2s_old, A1_E_2s_old, A2_E_2s_old = [_[2], _[3], _[4]]
    else: gamma_E_2s_old, A_E_2s_old = [_[0], _[1]]

    # EPTA results with 1-sigma confidence
    file = 'EPTA_singlePL_1s_EP.csv'
    _ = get_gamma_A(file, beta_b=beta_b, Omega_b=Omega_b,
                    return_all=return_all, fill=fill_1s,
                    fref=fref, plot=plot, color='purple', alpha=.8)
    if return_all:
        gamma_E_1s, A1_E_1s, A2_E_1s = [_[2], _[3], _[4]]
    else: gamma_E_1s, A_E_1s = [_[0], _[1]]

    # EPTA results with 2-sigma confidence
    file = 'EPTA_singlePL_2s_EP.csv'
    _ = get_gamma_A(file, beta_b=beta_b, Omega_b=Omega_b,
                    return_all=return_all, fill=fill_2s,
                    fref=fref, plot=plot, color='purple')
    if return_all:
        gamma_E_2s, A1_E_2s, A2_E_2s = [_[2], _[3], _[4]]
    else: gamma_E_2s, A_E_2s = [_[0], _[1]]

    # IPTA results (DR2) with 1-sigma confidence
    file = 'IPTA_singlePL_1s.csv'
    _ = get_gamma_A(file, beta_b=beta_b, Omega_b=Omega_b,
                    return_all=return_all, fill=fill_1s,
                    fref=fref, plot=plot, color='black', alpha=.8)
    if return_all:
        gamma_I_1s, A1_I_1s, A2_I_1s = [_[2], _[3], _[4]]
    else: gamma_I_1s, A_I_1s = [_[0], _[1]]

    # IPTA results (DR2) with 2-sigma confidence
    file = 'IPTA_singlePL_2s.csv'
    _ = get_gamma_A(file, beta_b=beta_b, Omega_b=Omega_b,
                    return_all=return_all, fill=fill_2s,
                    fref=fref, plot=plot, color='black')
    if return_all:
        gamma_I_2s, A1_I_2s, A2_I_2s = [_[2], _[3], _[4]]
    else: gamma_I_2s, A_I_2s = [_[0], _[1]]

    if return_all:
        return (gamma_NG_sPL_2s, A1_NG_sPL_2s, A2_NG_sPL_2s, gamma_NG_sPL_1s,
                A1_NG_sPL_1s, A2_NG_sPL_1s, gamma_NG_bPL_2s, A1_NG_bPL_2s,
                A2_NG_bPL_2s, gamma_NG_bPL_1s, A1_NG_bPL_1s, A2_NG_bPL_1s,
                gamma_P_b_2s, A1_P_b_2s, A2_P_b_2s, gamma_P_b_1s, A1_P_b_1s,
                A2_P_b_1s, gamma_E_2s, A1_E_2s_old, A2_E_2s_old, gamma_E_1s_old,
                A1_E_1s_old, A2_E_1s_old, gamma_E_2s, A1_E_2s, A2_E_2s,
                gamma_E_1s, A1_E_1s, A2_E_1s, gamma_I_2s, A1_I_2s, A2_I_2s,
                gamma_I_1s, A1_I_1s, A2_I_1s)
    else:
        return (gamma_NG_sPL_2s, A_NG_sPL_2s, gamma_NG_sPL_1s, A_NG_sPL_1s,
            gamma_NG_bPL_2s, A_NG_bPL_2s, gamma_NG_bPL_1s, A_NG_bPL_1s,
            gamma_P_b_2s, A_P_b_2s, gamma_P_b_1s, A_P_b_1s, gamma_E_2s_old,
            A_E_2s_old, gamma_E_1s_old, A_E_1s_old, gamma_E_2s, A_E_2s,
            gamma_E_1s, A_E_1s, gamma_I_2s, A_I_2s, gamma_I_1s, A_I_1s)
    
def single_CP(f, gamma, A1, A2, beta=0, broken=True, plot=False, alpha=.1,
              T=12.5, color='blue'):

    """
    Function to compute the GW spectral maxima and minima at each frequency
    from the allowed amplitudes and slopes from the PTA results.

    Arguments:
        f -- array of frequencies to compute the background
        gamma -- slopes of the power spectral density allowed by PTA data
        A1 -- minimum amplitude allowed by PTA data for the slope gamma
        A2 -- maximum amplitude allowed by PTA data for the slope gamma
        beta -- slope used to compute the GW energy density spectrum
                (default is 0)
        broken -- option to plot the broken power law fit (default is True)
        plot -- option to plot the resulting spectra (default False)
        alpha -- transparency of the plot (default 1)
        color -- color of the plot (default 'blue')
        T -- duration of the mission (default is 12.5 years for NANOGrav)

    Returns:
        Sf -- power spectral density as a function of frequency
        hc -- characteristic strain spectrum
        OmGW -- GW energy density spectrum
            a, b indicate minima and maxima at each frequency
    """

    T = T*u.yr
    T = T.to(u.s)
    f1yr = 1/u.yr
    f1yr = f1yr.to(u.Hz)
    gam = 5 - beta
    inside = True
    if gam < min(gamma) or gam > max(gamma): inside = False
    As = np.linspace(np.interp(gam, gamma, A1), np.interp(gam, gamma, A2), 10)
    if not inside: As = np.zeros(10)
    ACP = np.sqrt(As[0]*As[-1])
    Sf_a, hc_a, OmGW_a = Sf_PL_PTA(As[0], f, gam, broken=broken)
    Sf_b, hc_b, OmGW_b = Sf_PL_PTA(As[-1], f, gam, broken=broken)
    Sf_c, hc_c, OmGW_c = Sf_PL_PTA(ACP, f, gam, broken=broken)

    if plot:
        if alpha < 1:
            plt.plot(f, np.sqrt(Sf_c/T), color=color)
            plt.plot(f, np.sqrt(Sf_a/T), color=color, ls='dashed', lw=1)
            plt.plot(f, np.sqrt(Sf_b/T), color=color, ls='dashed', lw=1)
        plt.fill_between(f, np.sqrt(Sf_a/T), np.sqrt(Sf_b/T),
                         alpha=alpha, color=color)

    return Sf_a, hc_a, OmGW_a, Sf_b, hc_b, OmGW_b

def CP_delay(betas, colors=[], obs='NANOGrav_brokenPL_1s',
             plot=False, alpha=.1):

    # Compute the data of NANOGrav and PPTA
    _ = read_PTA_data(beta_b=False, Omega_b=False, plot=False,
                      return_all=True)
    gamma_NG_sPL_2s, A1_NG_sPL_2s, A2_NG_sPL_2s = [_[0], _[1], _[2]]
    gamma_NG_sPL_1s, A1_NG_sPL_1s, A2_NG_sPL_1s = [_[3], _[4], _[5]]
    gamma_NG_bPL_2s, A1_NG_bPL_2s, A2_NG_bPL_2s = [_[6], _[7], _[8]]
    gamma_NG_bPL_1s, A1_NG_bPL_1s, A2_NG_bPL_1s = [_[9], _[10], _[11]]
    gamma_P_b_2s, A1_P_b_2s, A2_P_b_2s = [_[12], _[13], _[14]]
    gamma_P_b_1s, A1_P_b_1s, A2_P_b_1s = [_[15], _[16], _[17]]
    gamma_E_2s_old, A1_E_2s_old, A2_E_2s_old = [_[18], _[19], _[20]]
    gamma_E_1s_old, A1_E_1s_old, A2_E_1s_old = [_[21], _[22], _[23]]
    gamma_E_2s, A1_E_2s, A2_E_2s = [_[24], _[25], _[26]]
    gamma_E_1s, A1_E_1s, A2_E_1s = [_[27], _[28], _[29]]
    gamma_I_2s, A1_I_2s, A2_I_2s = [_[30], _[31], _[32]]
    gamma_I_1s, A1_I_1s, A2_I_1s = [_[33], _[34], _[35]]

    # time of data considered for each of the PTA collaborations (in years)
    if 'NANOGrav' in obs: Tdata = T_NG
    if 'PPTA' in obs: Tdata = T_P
    if 'EPTA' in obs: Tdata = T_E
    if 'IPTA' in obs: Tdata = T_I

    # duration of NANOGrav observations
    broken = True
    obsss = ['NANOGrav_brokenPL_1s', 'NANOGrav_brokenPL_2s',
             'NANOGrav_singlePL_1s', 'NANOGrav_singlePL_2s', 'PPTA_blue_1s',
             'PPTA_blue_2s', 'EPTA_1s_old', 'EPTA_2s_old', 'EPTA_1s',
             'EPTA_2s', 'IPTA_1s', 'IPTA_2s']
    if obs == 'NANOGrav_brokenPL_1s':
        gamma = gamma_NG_bPL_1s
        A1 = A1_NG_bPL_1s
        A2 = A2_NG_bPL_1s
    elif obs == 'NANOGrav_brokenPL_2s':
        gamma = gamma_NG_bPL_2s
        A1 = A1_NG_bPL_2s
        A2 = A2_NG_bPL_2s
    elif obs == 'NANOGrav_singlePL_1s':
        broken = False
        gamma = gamma_NG_sPL_1s
        A1 = A1_NG_sPL_1s
        A2 = A2_NG_sPL_1s
    elif obs == 'NANOGrav_singlePL_2s':
        broken = False
        gamma = gamma_NG_sPL_2s
        A1 = A1_NG_sPL_2s
        A2 = A2_NG_sPL_2s
    elif obs == 'PPTA_blue_1s':
        broken = False
        gamma = gamma_P_b_1s
        A1 = A1_P_b_1s
        A2 = A2_P_b_1s
    elif obs == 'PPTA_blue_2s':
        broken = False
        gamma = gamma_P_b_2s
        A1 = A1_P_b_2s
        A2 = A2_P_b_2s
    elif obs == 'EPTA_1s_old':
        broken = False
        gamma = gamma_E_1s_old
        A1 = A1_E_1s_old
        A2 = A2_E_1s_old
    elif obs == 'EPTA_2s_old':
        broken = False
        gamma = gamma_E_2s_old
        A1 = A1_E_2s_old
        A2 = A2_E_2s_old
    elif obs == 'EPTA_1s':
        broken = False
        gamma = gamma_E_1s
        A1 = A1_E_1s
        A2 = A2_E_1s
    elif obs == 'EPTA_2s':
        broken = False
        gamma = gamma_E_2s
        A1 = A1_E_2s
        A2 = A2_E_2s
    elif obs == 'IPTA_1s':
        broken = False
        gamma = gamma_I_1s
        A1 = A1_I_1s
        A2 = A2_I_1s
    elif obs == 'IPTA_2s':
        broken = False
        gamma = gamma_I_2s
        A1 = A1_I_2s
        A2 = A2_I_2s
    else:
        gamma = 0
        A1 = 0
        A2 = 0
        print('The selected observation data \' should be one of',
                ' the available ', obsss)

    fmin = 1/Tdata/u.yr
    fmin = fmin.to(u.Hz)
    f = np.logspace(np.log10(fmin.value), np.log10(7e-8), 300)*u.Hz

    if len(colors) == 0: colors = ['blue']*len(betas)
    Sf_a = np.zeros((len(betas), len(f)))
    hc_a = np.zeros((len(betas), len(f)))
    OmGW_a = np.zeros((len(betas), len(f)))
    Sf_b = np.zeros((len(betas), len(f)))
    hc_b = np.zeros((len(betas), len(f)))
    OmGW_b = np.zeros((len(betas), len(f)))
    Sf_c = np.zeros((len(betas), len(f)))
    hc_c = np.zeros((len(betas), len(f)))
    OmGW_c = np.zeros((len(betas), len(f)))
    for i in range(0, len(betas)):
        Sf_a[i, :], hc_a[i, :], OmGW_a[i, :], Sf_c[i, :], \
                    hc_c[i, :], OmGW_c[i, :] = \
                            single_CP(f, gamma, A1, A2, beta=betas[i],
                                      color=colors[i], plot=plot,
                                      broken=broken, alpha=alpha,
                                      T=Tdata)

    return f, Sf_a, hc_a, OmGW_a, Sf_c, hc_c, OmGW_c

def OmGW_PTA(betas, ff='Om'):

    """
    Function that returns the GW energy density minima and maxima for a range
    of slopes using the reported amplitudes by NANOGrav and PPTA.

    Arguments:
        betas -- array with values of the slopes beta
    """

    import spectra as sp

    # NANOGrav broken PL
    _ = CP_delay(betas, obs='NANOGrav_brokenPL_2s')
    f_bPL_NG = _[0]
    if ff == 'Om':
        OmGW_a = _[3]
        OmGW_b = _[6]
    if ff == 'Sf':
        OmGW_a = _[1]
        OmGW_b = _[4]
    if ff == 'hc':
        OmGW_a = _[2]
        OmGW_b = _[5]
    min_OmGW_bPL_NG, max_OmGW_bPL_NG = sp.get_min_max(f_bPL_NG, OmGW_a,
                                                      OmGW_b)

    # NANOGrav single PL
    _ = CP_delay(betas, obs='NANOGrav_singlePL_2s')
    f_sPL_NG = _[0]
    if ff == 'Om':
        OmGW_a = _[3]
        OmGW_b = _[6]
    if ff == 'Sf':
        OmGW_a = _[1]
        OmGW_b = _[4]
    if ff == 'hc':
        OmGW_a = _[2]
        OmGW_b = _[5]
    min_OmGW_sPL_NG, max_OmGW_sPL_NG = sp.get_min_max(f_sPL_NG, OmGW_a,
                                                      OmGW_b)

    # PPTA single PL
    _ = CP_delay(betas, obs='PPTA_blue_2s')
    f_sPL_P = _[0]
    if ff == 'Om':
        OmGW_a = _[3]
        OmGW_b = _[6]
    if ff == 'Sf':
        OmGW_a = _[1]
        OmGW_b = _[4]
    if ff == 'hc':
        OmGW_a = _[2]
        OmGW_b = _[5]
    min_OmGW_sPL_P, max_OmGW_sPL_P = sp.get_min_max(f_sPL_P, OmGW_a,
                                                    OmGW_b)

    # EPTA single PL
    _ = CP_delay(betas, obs='EPTA_2s_old')
    f_sPL_E_old = _[0]
    if ff == 'Om':
        OmGW_a = _[3]
        OmGW_b = _[6]
    if ff == 'Sf':
        OmGW_a = _[1]
        OmGW_b = _[4]
    if ff == 'hc':
        OmGW_a = _[2]
        OmGW_b = _[5]
    min_OmGW_sPL_E_old, max_OmGW_sPL_E_old = sp.get_min_max(f_sPL_E_old,
                                                            OmGW_a, OmGW_b)

    # EPTA single PL
    _ = CP_delay(betas, obs='EPTA_2s')
    f_sPL_E = _[0]
    if ff == 'Om':
        OmGW_a = _[3]
        OmGW_b = _[6]
    if ff == 'Sf':
        OmGW_a = _[1]
        OmGW_b = _[4]
    if ff == 'hc':
        OmGW_a = _[2]
        OmGW_b = _[5]
    min_OmGW_sPL_E, max_OmGW_sPL_E = sp.get_min_max(f_sPL_E, OmGW_a,
                                                    OmGW_b)

    # IPTA single PL
    _ = CP_delay(betas, obs='IPTA_2s')
    f_sPL_I = _[0]
    if ff == 'Om':
        OmGW_a = _[3]
        OmGW_b = _[6]
    if ff == 'Sf':
        OmGW_a = _[1]
        OmGW_b = _[4]
    if ff == 'hc':
        OmGW_a = _[2]
        OmGW_b = _[5]
    min_OmGW_sPL_I, max_OmGW_sPL_I = sp.get_min_max(f_sPL_I, OmGW_a,
                                                    OmGW_b)

    return (f_bPL_NG, min_OmGW_bPL_NG, max_OmGW_bPL_NG, f_sPL_NG,
            min_OmGW_sPL_NG, max_OmGW_sPL_NG, f_sPL_P, min_OmGW_sPL_P,
            max_OmGW_sPL_P, f_sPL_E_old, min_OmGW_sPL_E_old, max_OmGW_sPL_E_old,
            f_sPL_E, min_OmGW_sPL_E, max_OmGW_sPL_E, f_sPL_I, min_OmGW_sPL_I,
            max_OmGW_sPL_I)


def Sf_PL_PTA(A, f, gamma, fbend=0, kappa=0.1, broken=True):

    """
    Function to compute the power spectral density Sf from the amplitude A
    of the background of GW characteristic strain.

    See A. Roper Pol, C. Caprini, A. Neronov, D. Semikoz,
    "The gravitational wave signal from primordial magnetic fields in the Pulsar
    Timing Array frequency band," https://arxiv.org/abs/2201.05630 (2022).

    h_c(f) = A (f/fyr)^((3 - gamma)/2)

    Arguments:
        A -- amplitudes of the GW background
        f -- frequency array
        gamma -- negative slope of the power spectral density
        fbend -- bend frequency of the broken power law model (default is
                 1.035e-8 Hz)
        kappa -- smoothing parameter of the broken power law model (default 0.1)
        broken -- allows to use a broken power law model, instead of a single
                  power law (default True)

    Returns:
        Sf -- power spectral density as a function of frequency
        hc -- characteristic strain spectrum
        OmGW -- GW energy density spectrum
    """

    f1yr = 1/u.yr
    f1yr = f1yr.to(u.Hz)
    if fbend == 0: fbend = 1.035e-8*u.Hz
    if broken:
        Sf = A**2/12/np.pi**2*(f.value/f1yr.value)**(-gamma)* \
            (1 + (f.value/fbend.value)**(1/kappa))**(kappa*gamma)*f1yr**(-3)
    else: Sf = A**2/12/np.pi**2*(f.value/f1yr.value)**(-gamma)*f1yr**(-3)
    hc = cosmoGW.hc_Sf(f, Sf)
    OmGW = cosmoGW.hc_OmGW(f, hc, d=-1)

    return Sf, hc, OmGW

def plot_PTA_all(ff='tdel', betas=[], lines=True, alp_bl=0.3, alp_E=0.2, alp_g=0.1,
                 alp_P=0.1, alp_I=0.1, plot=True, ret=False):

    """
    Function that overplots the GW spectra Omega_GW (f) = Omyr (f/fyr)^beta
    for the values of beta and Omyr reported by the PTA collaborations.

    Arguments:
        ff -- option to chose what function to plot (default 'tdel' for time
              delay in seconds, 'Sf' for power spectral density, 'hc' for
              characteristic amplitude spectrum, and 'Om' for GW energy density
              spectrum, from -2 to 5)
        betas -- range of slopes considered for the plot (default is all
                 possible values)
        lines -- option to explicitly plot the boundaring lines on top of the
                 allowed region.
        alp_bl, _E, _g, _P, _I -- alpha of the contour plots for NG broken PL
                                  (bl for blue), NG single PL (g for green),
                                  EPTA (E), PPTA (P) and IPTA (I)
        plot -- option to plot the results
        ret -- option to return the resulting OmGW spectra
    """

    if len(betas) == 0: betas = np.linspace(-2, 5, 100)

    if ff == 'tdel': ff2 = 'Sf'
    else: ff2 = ff

    (f_bPL_NG, min_OmGW_bPL_NG, max_OmGW_bPL_NG, f_sPL_NG, min_OmGW_sPL_NG,
    max_OmGW_sPL_NG, f_sPL_P, min_OmGW_sPL_P, max_OmGW_sPL_P, f_sPL_E_old,
    min_OmGW_sPL_E_old, max_OmGW_sPL_E_old, f_sPL_E, min_OmGW_sPL_E,
    max_OmGW_sPL_E, f_sPL_I, min_OmGW_sPL_I, max_OmGW_sPL_I) = \
                    OmGW_PTA(betas, ff=ff2)

    # duration of observations for NANOGrav
    TNG = T_NG*u.yr
    TNG = TNG.to(u.s)
    # duration of observations for PPTA
    TP = T_P*u.yr
    TP = TP.to(u.s)
    # duration of observations for EPTA
    TE = T_E*u.yr
    TE = TE.to(u.s)
    # duration of observations for IPTA
    TI = T_I*u.yr
    TI = TI.to(u.s)

    if ff == 'tdel':
        min_OmGW_bPL_NG = np.sqrt(min_OmGW_bPL_NG/T_NG)
        max_OmGW_bPL_NG = np.sqrt(max_OmGW_bPL_NG/T_NG)
        min_OmGW_sPL_NG = np.sqrt(min_OmGW_sPL_NG/T_NG)
        max_OmGW_sPL_NG = np.sqrt(max_OmGW_sPL_NG/T_NG)
        min_OmGW_sPL_P = np.sqrt(min_OmGW_sPL_P/T_P)
        max_OmGW_sPL_P = np.sqrt(max_OmGW_sPL_P/T_P)
        min_OmGW_sPL_E_old = np.sqrt(min_OmGW_sPL_E_old/T_E)
        max_OmGW_sPL_E_old = np.sqrt(max_OmGW_sPL_E_old/T_E)
        min_OmGW_sPL_E = np.sqrt(min_OmGW_sPL_E/T_E)
        max_OmGW_sPL_E = np.sqrt(max_OmGW_sPL_E/T_E)
        min_OmGW_sPL_I = np.sqrt(min_OmGW_sPL_I/T_I)
        max_OmGW_sPL_I = np.sqrt(max_OmGW_sPL_I/T_I)

    flim = 1.25e-8*u.Hz
    good = np.where(f_bPL_NG.value < flim.value)
    
    if plot:
        # NG single PL
        plt.fill_between(f_sPL_NG.value, min_OmGW_sPL_NG, max_OmGW_sPL_NG,
                         color='darkgreen', alpha=alp_g)
        if lines:
            plt.plot(f_sPL_NG.value, min_OmGW_sPL_NG, color='darkgreen', lw=2)
            plt.plot(f_sPL_NG.value, max_OmGW_sPL_NG, color='darkgreen', lw=2)
        
        # NG broken PL
        plt.fill_between(f_bPL_NG[good].value, min_OmGW_bPL_NG[good],
                         max_OmGW_bPL_NG[good], color='blue', alpha=alp_bl,
                         label='NG bPL')
        if lines:
            plt.plot(f_bPL_NG[good].value, min_OmGW_bPL_NG[good], color='blue', lw=2)
            plt.plot(f_bPL_NG[good].value, max_OmGW_bPL_NG[good], color='blue', lw=2)
    
        # PPTA
        plt.fill_between(f_sPL_P.value, min_OmGW_sPL_P, max_OmGW_sPL_P,
                         color='red', alpha=alp_P)
        if lines:
            plt.plot(f_sPL_P.value, min_OmGW_sPL_P, color='red', lw=2)
            plt.plot(f_sPL_P.value, max_OmGW_sPL_P, color='red', lw=2)

        # EPTA
        plt.fill_between(f_sPL_E.value, min_OmGW_sPL_E, max_OmGW_sPL_E,
                         color='purple', alpha=alp_E)
        if lines:
            plt.plot(f_sPL_E.value, min_OmGW_sPL_E, color='purple', lw=2)
            plt.plot(f_sPL_E.value, max_OmGW_sPL_E, color='purple', lw=2)

        # IPTA
        plt.fill_between(f_sPL_I.value, min_OmGW_sPL_I, max_OmGW_sPL_I,
                         color='black', alpha=alp_I)
        if lines:
            plt.plot(f_sPL_I.value, min_OmGW_sPL_I, color='black', lw=2)
            plt.plot(f_sPL_I.value, max_OmGW_sPL_I, color='black', lw=2)

    if ret:
        return (f_bPL_NG, min_OmGW_bPL_NG, max_OmGW_bPL_NG, f_sPL_NG,
                min_OmGW_sPL_NG, max_OmGW_sPL_NG, f_sPL_P, min_OmGW_sPL_P,
                max_OmGW_sPL_P, f_sPL_E_old, min_OmGW_sPL_E_old, max_OmGW_sPL_E_old,
                f_sPL_E, min_OmGW_sPL_E, max_OmGW_sPL_E, f_sPL_I, min_OmGW_sPL_I,
                max_OmGW_sPL_I)