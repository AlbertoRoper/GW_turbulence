"""
GW_analytical.py is a Python routine that contains the functions to make
calculations related to the analytical description of GW backgrounds produced
by MHD turbulence in the early universe.
It also includes some general mathematical functions that are useful.
"""

import numpy as np
import matplotlib.pyplot as plt
import plot_sets
import os

dir0 = os.getcwd()

def compute_Pi_from_numerical(k, mag, Np=3000, Nk=60, plot=False,
                              extend=False, largek=3, smallk=-3):

    """
    Function that computes the projected stress spectrum Pi(k) from the
    magnetic spectrum under the assumption of Gaussianity.

    Arguments:
        k -- array of wave numbers
        mag -- array of magnetic spectrum values
        Np -- number of discretizations in the wave number p to be numerically
              integrated
        Nk -- number of discretizations of k to be used for the computation of
              the final spectrum
        plot -- option to plot the interpolated magnetic spectrum for debugging
                purposes (default False)
        extend -- option to extend the array of wave numbers of the resulting
                  Pi spectrum compared to that of the given magnetic spectrum
                  (default False)

    Returns:
        PiM -- spectrum of the magnetic stress
        kp -- final array of wave numbers
    """

    p = np.logspace(np.log10(k[0]), np.log10(k[-1]), Np)
    kp = np.logspace(np.log10(k[0]), np.log10(k[-1]), Nk)
    if extend:
        Nk = int(Nk/6)
        kp = np.logspace(-smallk, np.log10(k[0]), Nk)
        kp = np.append(kp, np.logspace(np.log10(k[0]), np.log10(k[-1]),
                       4*Nk))
        kp = np.append(kp, np.logspace(np.log10(k[-1]), largek, Nk))

    Nphi = 300
    phi = np.linspace(0, np.pi, Nphi)
    kij, pij, phij = np.meshgrid(kp, p, phi, indexing='ij')
    PB_p = np.interp(p, k, mag)
    if plot:
        plt.plot(p, PB_p)
        plt.xscale('log')
        plt.yscale('log')

    kmp_mod_2 = pij**2 + kij**2 - 2*pij*kij*np.cos(phij)
    kmp_mod = np.sqrt(kmp_mod_2)
    kmp_mod_2[np.where(kmp_mod_2==0)] = 1e-20
    PB_kmp = np.interp(kmp_mod, k, mag)/kmp_mod_2
    Pi0 = 1 + np.cos(phij)**2
    Pi1 = 1 + (kij - pij*np.cos(phij))**2/kmp_mod_2
    Piph_B = Pi0*Pi1*PB_kmp*np.sin(phij)
    Pi_imB = np.trapz(Piph_B, phi, axis=2)
    Pi_p_BB = Pi_imB*PB_p
    PiM = .5*np.trapz(Pi_p_BB, p, axis=1)

    return kp, PiM

def compute_Pihel_from_numerical(k, mag, Np=3000, Nk=60, plot=False,
                                 extend=False, largek=3, smallk=-3):

    """
    Function that computes the projected stress spectrum Pi(k) from the
    magnetic helicity spectrum under the assumption of Gaussianity.

    Arguments:
        k -- array of wave numbers
        mag -- array of magnetic spectrum values, it should be given in the
               form (1/2) k H_M (k)
        Np -- number of discretizations in the wave number p to be numerically
              integrated
        Nk -- number of discretizations of k to be used for the computation of
              the final spectrum
        plot -- option to plot the interpolated magnetic spectrum for debugging
                purposes (default False)
        extend -- option to extend the array of wave numbers of the resulting
                  Pi spectrum compared to that of the given magnetic spectrum
                  (default False)

    Returns:
        PiMhel -- spectrum of the magnetic stress
        kp -- final array of wave numbers
    """

    p = np.logspace(np.log10(k[0]), np.log10(k[-1]), Np)
    kp = np.logspace(np.log10(k[0]), np.log10(k[-1]), Nk)
    if extend:
        Nk = int(Nk/6)
        kp = np.logspace(-smallk, np.log10(k[0]), Nk)
        kp = np.append(kp, np.logspace(np.log10(k[0]), np.log10(k[-1]),
                       4*Nk))
        kp = np.append(kp, np.logspace(np.log10(k[-1]), largek, Nk))

    Nphi = 300
    phi = np.linspace(0, np.pi, Nphi)
    kij, pij, phij = np.meshgrid(kp, p, phi, indexing='ij')
    PB_p = np.interp(p, k, mag)
    if plot:
        plt.plot(p, PB_p)
        plt.xscale('log')
        plt.yscale('log')

    kmp_mod_2 = pij**2 + kij**2 - 2*pij*kij*np.cos(phij)
    kmp_mod = np.sqrt(kmp_mod_2)
    kmp_mod_2[np.where(kmp_mod_2==0)] = 1e-20
    kmp_mod[np.where(kmp_mod==0)] = 1e-20
    PB_kmp = np.interp(kmp_mod, k, mag)/kmp_mod_2
    Pi0 = 4*np.cos(phij)
    Pi1 = (pij*np.cos(phij) - kij)/kmp_mod
    Piph_B = Pi0*Pi1*PB_kmp*np.sin(phij)
    Pi_imB = np.trapz(Piph_B, phi, axis=2)
    Pi_p_BB = Pi_imB*PB_p
    PiM = .5*np.trapz(Pi_p_BB, p, axis=1)

    return kp, PiM

def fit_smoothed_bPL(x, A=1, a=4, b=5/3, alp=2, xc=1, Omega=False):

    """
    Function that returns the value of the smoothed broken power law (PL) model
    for the magnetic spectrum, which has the form
        y = A(1 + D)^(1/alp)*(k/kp)^a/(1 + D(k/kp)^(a + b))^(1/alp),
    where D = a/b

    Arguments:
        x -- values of x
        A -- amplitude of the spectrum at the peak
        a -- slope of the smoothed broken PL at low wave numbers, k^a
        b -- slope of the smoothed broken PL at high wave numbers, k^(-b)
        alp -- smoothness of the transition from one power law to the other
        xc -- spectral peak, i.e., position of the break from k^a to k^(-b)
        Omega -- option to use the integrated energy desnity as the input A

    Returns:
        y -- value of the spectrum
    """

    if Omega:
        A = A/xc/A_alpha(alp=alp, a=a, b=b)

    if b == 0:
        D1 = 0
        D2 = 1.
    else:
        D1 = a/b
        D2 = a/b
    if a == 0:
        D2 = 1.

    y = A*(1 + D1)**(1/alp)*(x/xc)**a/(1 + D2* \
            (x/xc)**((a + b)*alp))**(1/alp)

    return y

def value_Pi0_sbPL(alp=2):

    """
    Function that computes the value of the Pi spectrum at k = 0 analytically
    using the smoothed broken power law model defined in fit_smoothed_bPL()
    with the slopes a = 4 and b = 5/3 and the amplitude A = 1 at the spectral
    peak.

    Arguments:
        alpha -- smoothing parameter

    Returns:
        C -- amplitude of the Pi spectrum at k = 0
    """

    import math as m

    C = 17**(2/alp)*5**(-13/17/alp - 1)*3**(-21/17/alp - 1)* \
            2**(-42/17/alp + 2)
    C *= m.gamma(1 + 21/17/alp)
    C *= m.gamma(13/17/alp)
    C *= 1/m.gamma(2/alp)

    return C

def A_alp(alp=2):

    import math as m

    A = (1 + 12/5)**(1/alp)*2**(-30/17/alp - 1)* \
            3**(-15/17/alp + 1)*5**(15/17/alp)
    A *= m.gamma(1 + 2/17/alp)*m.gamma(15/17/alp)/m.gamma(1/alp)

    return A

def Iab_n_alpha(alp=2, a=4, b=5/3, n=0):

    import math as m

    D = a/b
    Iab = D**(-(a + 1 + n)/alp/(a + b))/alp/(a + b)/m.gamma(1/alp)
    Iab *= m.gamma((b - 1 - n)/alp/(a + b))*m.gamma((a + 1 + n)/alp/(a + b))

    return Iab

# def A_alpha(alp=2, a=4, b=5/3):
#
#     import math as m
#
#     D = a/b
#     A = ((1 + D)*D**(-(a + 1)/(a + b)))**(1/alp)
#     A *= m.gamma((b - 1)/alp/(a + b))
#     A *= m.gamma((a + 1)/alp/(a + b))
#     A *= 1/alp/(a + b)/m.gamma(1/alp)
#
#     return A

def A_alpha(alp=2, a=4, b=5/3):

    D = a/b
    A = (1 + D)**(1/alp)*Iab_n_alpha(alp=alp, a=a, b=b, n=0)

    return A

# def B_alpha(alp=2, a=4, b=5/3):
#
#     import math as m
#
#     D = a/b
#     B = D**(1/alp/(a + b))*m.gamma(b/alp/(a + b))
#     B *= m.gamma(a/alp/(a + b))
#     B *= 1/m.gamma((b - 1)/alp/(a + b))
#     B *= 1/m.gamma((a + 1)/alp/(a + b))
#
#     return B

def B_alpha(alp=2, a=4, b=5/3):

    D = a/b
    B = Iab_n_alpha(alp=alp, a=a, b=b, n=-1)
    B = B/Iab_n_alpha(alp=alp, a=a, b=b, n=0)

    return B

def AB2_alpha(alp=2, a=4, b=5/3):

    import math as m

    D = a/b
    alp2 = 1/alp/(a + b)
    A2 = m.gamma(alp2*(a + b))/alp2/m.gamma(alp2*(b - 1))/m.gamma(alp2*(a + 1))
    B2 = m.gamma(alp2*(b - 1))*m.gamma(alp2*(a + 1))
    B2 *= 1/m.gamma(alp2*b)/m.gamma(alp2*a)

    return A2, B2

def C_alpha(alp=2, a=4, b=5/3):

    D = a/b
    C = 28/15*(1 + D)**(2/alp)*Iab_n_alpha(a=2*a, b=2*b, alp=alp/2, n=-2)

    return C

def Cp_alpha(alp=2, a=4, b=5/3):

    D = a/b
    A = A_alpha(alp=alp, a=a, b=b)

    return 8/3*((1 + D)/D)**(1/alp)*A

def Cp_alpha_2(alp=2, a=4, b=5/3):

    D = a/b
    Cp = 10/7/Iab_n_alpha(a=2*a, b=2*b, alp=alp/2, n=-2)
    Cp *= Iab_n_alpha(a=a, b=b, alp=alp, n=0)*D**(-1/alp)

    return Cp

def K_Cp(alp=2, a=4, b=5/3):

    Cp = Cp_alpha_2(alp=alp, a=a, b=b)

    return Cp**(1/(b + 2))

def alpha2(alp=2, a=4, b=5/3):

    return 1./alp/(a + b)

def compute_Piref(ks=[], alp=2, save=False, str_alp='0'):

    """
    Function that computes the spectrum Pi using a magnetic spectrum
    with amplitude and peak values of 1 and stores it as Pi_ref_alpha''.csv

    Arguments:
        ks -- array of wave numbers to define the magnetic spectrum
              (default is 1000 points logarithmically distributed
              from 10^{-3} to 10^3)
        alp -- alpha (smoothing parameter of the magnetic field spectrum)
               value used for the computation
        save -- option to save the resulting Pi in a file

    Returns:
        kf -- array of wave numbers for the spectrum Pi
        Pi -- values of Pi compensated by the factor C = Pi(k = 0)
    """

    import pandas as pd

    if len(ks) == 0: ks = np.logspace(-3, 3, 1000)
    EM_test = fit_smoothed_bPL(ks, alp=alp)
    kfPi, Pi = compute_Pi_from_numerical(ks, EM_test, Np=1500,
                                         Nk=300, plot=True)
    C = value_Pi0_sbPL(alp=alp)
    Pi /= C

    if save:
        df = pd.DataFrame({'k': kfPi, 'Pi':Pi})
        fl = dir0 + 'analytical/Pi_ref_alpha%.1f.csv'%alp
        if str_alp != '0': fl = dir0 + \
                '/analytical/Pi_ref_alpha%s.csv'%str_alp
        df.to_csv(fl)
        print('Pi_ref saved in %s'%fl)

    return kfPi, Pi

def read_Pi(dir='', str_alp='2'):

    """
    Function that reads the spectrum Pi from a file of the type
    Pi_ref_alpha'str_alpha'.csv previously generated with compute_Piref.

    Arguments:
        str_alp -- string that defines the name of the specific file to be
                   read
        dir -- directory where 'analytical' directory is located with the
               files in it (default uses dir0, which corresponds to the same
               directory as where the routine GW_analytical.py is stored)

    Returns:
        kf -- array of wave numbers for the spectrum Pi
        Pi -- values of Pi compensated by the factor C = Pi(k = 0)
    """

    import pandas as pd

    if dir == '': dir = dir0
    file = 'Pi_ref_alpha' + str_alp
    df = pd.read_csv(dir + '/analytical/' + file + '.csv')
    kf_Pi = np.array(df['k'])
    Pi_ref = np.array(df['Pi'])
    kf1 = np.logspace(3.001, 9, 1000)
    kf0 = np.logspace(-9, -2.999, 1000)
    Pi_ref1 = Pi_ref[-1]*(kf1/kf_Pi[-1])**(-11/3)
    Pi_ref0 = Pi_ref[0]*kf0**0
    kf_Pi = np.append(kf0, kf_Pi)
    kf_Pi = np.append(kf_Pi, kf1)
    Pi_ref = np.append(Pi_ref0, Pi_ref)
    Pi_ref = np.append(Pi_ref, Pi_ref1)

    return kf_Pi, Pi_ref

def shift_Pi(k, Pi, kpeak):

    """
    Function that returns Pi by shifting the position of its peak, i.e.,
    shifts a function in the x-axis.

    Arguments:
        k -- original wave number array
        Pi -- original Pi values
        kpeak -- position of the peak k (or value of x used for shifting)

    Returns:
        Pinew -- new Pi values after the shift
    """

    ks = k/kpeak
    Pinew = 10**np.interp(np.log10(ks), np.log10(k), np.log10(Pi))
    return Pinew

def plot_Pi_max(dir0='', ymin=1e-3, ymax=1e1, xmin=1e-1, xmax=20,
                str_alp='2', plot=True, txt=True):

    """
    Function that plots the reference Pi (i.e., for amplitude and position of
    the peak being 1) and shows the maximum values and positions of k*Pi,
    k^2*Pi, and k^3*Pi.

    Arguments:
        dir0 -- directory where the Pi_ref files are stored
                (dir0/analytical/Pi_ref_''.csv)
        ymin, ymax -- minimum and maximum y limits of the plot
        xmin, xmax -- minimum and maximum x limits of the plot
        str_alp -- string that indicates the name of the Pi_ref file to be
                   read (default is str_alp = '2', indicating a smoothing
                   parameter alpha = 2)
        plot -- option to plot the resulting Pi function (default True)
        txt -- option to print out the text with the results (default True)

    Returns:
        max_ks, Pi_max_ks -- k and Pi that correspond to maximum k*Pi
        max_ks2, Pi_max_ks2 -- k and Pi that correspond to maximum k^2*Pi
        max_ks3, Pi_max_ks3 -- k and Pi that correspond to maximum k^3*Pi
    """

    import spectra as sp

    kf_Pi, Pi_ref = read_Pi(str_alp=str_alp)
    max_ks, Pi_max_ks = sp.max_E_kf(kf_Pi, Pi_ref, exp=1)
    max_ks2, Pi_max_ks2 = sp.max_E_kf(kf_Pi, Pi_ref, exp=2)
    max_ks3, Pi_max_ks3 = sp.max_E_kf(kf_Pi, Pi_ref, exp=3)

    if plot:
        plt.figure(figsize=(8, 5))
        plt.plot(kf_Pi, Pi_ref, color='black')
        plt.vlines(1, ymin, ymax, color='black', ls='dashed', lw=.6)
        plt.hlines(1, xmin, xmax, color='black', ls='dashed', lw=.6)

        plt.vlines(max_ks, ymin, ymax, color='black', ls='dashed', lw=.6)
        plt.vlines(max_ks2, ymin, ymax, color='black', ls='dashed', lw=.6)
        plt.vlines(max_ks3, ymin, ymax, color='black', ls='dashed', lw=.6)

        plt.hlines(Pi_max_ks, xmin, xmax, color='black', ls='dashed', lw=.6)
        plt.hlines(Pi_max_ks2, xmin, xmax, color='black', ls='dashed', lw=.6)
        plt.hlines(Pi_max_ks3, xmin, xmax, color='black', ls='dashed', lw=.6)

        plt.xlim(xmin, xmax)
        plt.ylim(ymin, ymax)
        plt.xscale('log')
        plt.yscale('log')
        plot_sets.axes_lines()
        plt.yticks(np.logspace(-3, 1, 5))

    if txt:
        print('Maximum of k*Pi corresponds to Pi = ', Pi_max_ks,
              ' at k = ', max_ks)
        print('Maximum of k^2*Pi corresponds to Pi = ', Pi_max_ks2,
              ' at k = ', max_ks2)
        print('Maximum of k^3*Pi corresponds to Pi = ', Pi_max_ks3,
              ' at k = ', max_ks3)

    return max_ks, Pi_max_ks, max_ks2, Pi_max_ks2, max_ks3, Pi_max_ks3

def function_D(k, t, tfin=1e4, tini=1):

    """
    Function that computes the value of the function D(k, t) used in the
    analytical calculations of the GW energy density spectrum when assuming
    a constant sourcing stress spectrum, i.e., Pi(k, t1, t2) = Pi(k)

    Arguments:
        k -- array of wave numbers
        t -- array of times
        tini -- initial time of the turbulence sourcing (default 1)
        tfin -- final time of the turbulence sourcing

    Returns:
        D -- function D(k, t)
    """

    import scipy.special as spe

    tij, kij = np.meshgrid(t, k, indexing='ij')
    cost = np.cos(kij*tij)
    sint = np.sin(kij*tij)
    tij[np.where(tij>tfin)] = tfin
    si_t, ci_t = spe.sici(kij*tij)
    si_tini, ci_tini = spe.sici(kij*tini)
    aux1 = cost*(ci_t - ci_tini)
    aux2 = sint*(si_t - si_tini)
    D = aux1 + aux2

    return D

def function_D2_av(k, tfin=1e4, tini=1):

    """
    Function that computes the value of the function D(k, t) used in the
    analytical calculations of the GW energy density spectrum when assuming
    a constant sourcing stress spectrum, i.e., Pi(k, t1, t2) = Pi(k).
    It takes the average over oscillations considering very large times, i.e.,
    by shifting to present time.

    Arguments:
        k -- array of wave numbers
        t -- array of times
        tini -- initial time of the turbulence sourcing (default 1)
        tfin -- final time of the turbulence sourcing

    Returns:
        D2_av -- average value of the square of the function, D^2_av (k)
    """

    import scipy.special as spe

    si_t, ci_t = spe.sici(k*tfin)
    si_tini, ci_tini = spe.sici(k*tini)
    aux1 = (ci_t - ci_tini)**2
    aux2 = (si_t - si_tini)**2
    D2_av = aux1 + aux2

    return D2_av

def factor_FF(tini=1, tfin=1e4):

    FF = (1 + tini/tfin)**2
    return FF

def function_D2_env(k, tini=1, tfin=1e4, FF=0, A=1):

    # envelope at low wave numbers
    f1 = np.log(tfin)**2
    # envelope at high wave numbers
    f2 = A*np.log(tini + 1/k)**2
    # overshooting factor that appears for short tfin
    # (if FF is given as an input parameter, then not computed)
    if FF == 0: FF = factor_FF(tini=tini, tfin=tfin)

    # envelope
    f = np.minimum(f1, FF*f2)

    # position of the break
    diff = f1 - FF*f2
    ind = np.where(diff > 0)[0][0]
    k_br = k[ind]

    return f, k_br

def model_EM(k, E, ks=[], indt=0):

    """
    Function that fits the numerical results of the magnetic spectrum EM(k, t)
    over times to a smoothed broken power law (alpha = 2) using
    'fit_smoothed_bPL' and returns the fitting parameters.

    Arguments:
        k -- array of wave numbers of the original spectra
        E -- 2d array of spectral values (first index is time and second is
             wave number)
        ks -- new array of wave numbers in which to define the modelled spectra
              (default is to use same as k)
        indt -- allows to limit the times in which the modelling is performed,
                from time 0 to indt (default 0, i.e, only at initial time)

    Returns:
        mod_EMs -- modelled magnetic spectra using the numerical wave numbers k
        err_M -- error of the fitting model
        mod_EMs_ks -- modelled magnetic spectra using the new wave numbers ks
        As_EM -- amplitudes of the fitting model
        alps_EM -- smoothing parameter of the fitting model
        xcs_EM -- position of the peak of the fitting model
    """

    import scipy.optimize as opt

    if len(ks) == 0: ks = k

    mod_EMs_ks = np.zeros((indt + 1, len(ks)))
    mod_EMs = np.zeros((indt + 1, len(k)))
    As_EM = np.zeros(indt + 1)
    as_EM = np.zeros(indt + 1)
    bs_EM = np.zeros(indt + 1)
    alps_EM = np.zeros(indt + 1)
    xcs_EM = np.zeros(indt + 1)

    def fit_test(x, A, b, alp, xc):
        y = fit_smoothed_bPL(x, A=A, b=b, alp=alp, xc=xc)
        return y

    for i in range(0, indt + 1):

        bound = ((-np.inf, 0, -np.inf, 0),
                 (np.inf, np.inf, np.inf, 30))

        popt, pcov = opt.curve_fit(fit_test, k, E[i, :],
                                   p0=(1e-5, 5/3, 2, 10), bounds=bound)
        As_EM[i] = popt[0]
        as_EM[i] = 4
        bs_EM[i] = popt[1]
        alps_EM[i] = popt[2]
        xcs_EM[i] = popt[3]
        mod_EMs_ks[i, :] = fit_test(ks, popt[0], popt[1], popt[2], popt[3])
        mod_EMs[i, :] = fit_test(k, popt[0], popt[1], popt[2], popt[3])

    err_M = np.sqrt(np.array(np.trapz((mod_EMs - E[:indt + 1, :])**2, k,
                                       axis=1)/\
                                       np.array(np.trapz(E[:indt + 1, :]**2,
                                                         k, axis=1)),
                                                dtype='float'))

    return mod_EMs, err_M, mod_EMs_ks, As_EM, alps_EM, xcs_EM

def factM(k1, k2, E1, E2):

    """
    Function that computes the ratio between the amplitude of two magnetic
    spectra at initial time.

    Arguments:
        k1, k2 -- wave number arrays of runs 1 and 2
        E1, E2 -- spectral values arrays of runs 1 and 2

    Returns:
        factM -- ratio of amplitudes A2/A1
    """

    _ = model_EM(k1, E1)
    A1 = _[3]
    _ = model_EM(k2, E2)
    A2 = _[3]
    return A2/A1

def OmGW_from_Pi0(Pi0, k, t, tfin=1e4, tini=1):

    D = function_D(k, t, tfin=tfin, tini=tini)
    tij, kij = np.meshgrid(t, k, indexing='ij')
    tij, Pi0 = np.meshgrid(t, Pi0, indexing='ij')
    OmGW = 3*kij**3*D**2*Pi0
    EGW = 3*kij**2*D**2*Pi0

    return EGW, OmGW

def OmGW_from_Pi0_av(Pi0, k, tfin=1e4, tini=1):

    D2 = function_D2_av(k, tfin=tfin, tini=tini)
    OmGW = 3*k**3*D2*Pi0
    EGW = 3*k**2*D2*Pi0

    return EGW, OmGW

def OmGW_from_Pi0_env(Pi0, k, tini=1, tfin=1e4, FF=0, A=1):

    D2, k_br = function_D2_env(k, tini=tini, tfin=tfin, FF=FF, A=A)
    OmGW = 3*k**3*D2*Pi0
    EGW = 3*k**2*D2*Pi0

    return EGW, OmGW, k_br

def OmGW_from_OmM_kM_env(OmM, kM, k, tini=1, tfin=1e4, FF=0, A=1, multi=False):

    EM = OmM/A_alp()/kM
    if multi:
        NN = np.shape(OmM)
        OmGW = np.zeros((len(k), NN[0], NN[1]))
        EGW = np.zeros((len(k), NN[0], NN[1]))
        k_br = np.zeros((NN[0], NN[1]))
        kf_Pi, Pi_ref = read_Pi(str_alp='2')
        Pi0 = value_Pi0_sbPL()
        for i in range(0, NN[0]):
            Pi_ref_i = shift_Pi(kf_Pi, Pi_ref, kM[i, 0])
            for j in range(0, NN[1]):
                PiM = Pi0*Pi_ref_i*EM[i, j]**2/kM[i, j]
                PiM_k = 10**np.interp(np.log10(k), np.log10(kf_Pi),
                                      np.log10(PiM))
                EGW[:, i, j], OmGW[:, i, j], k_br[i, j] = \
                        OmGW_from_Pi0_env(PiM_k, k, tini=tini,
                                          tfin=tfin[i, j], FF=FF, A=A)
    else:
        kf_Pi, Pi = compute_PiM(EM, kM)
        PiM = 10**np.interp(np.log10(k), np.log10(kf_Pi),
                            np.log10(Pi))
        EGW, OmGW, k_br = OmGW_from_Pi0_env(PiM, k, tini=tini, tfin=tfin, FF=FF,
                                            A=A)

    return EGW, OmGW, k_br

def OmGW_from_OmM_kM_env_tfin_fit(OmM, kM, k, tini=1, FF=0, A=1, multi=False):

    if multi: kM_ij, OmM_ij = np.meshgrid(kM, OmM, indexing='ij')
    else:
        kM_ij = kM
        OmM_ij = OmM
    vA = np.sqrt(1.5*OmM_ij)
    dte = 1/kM_ij/vA
    dtfin = dtfin_dte(dte)
    EGW, OmGW, k_br = OmGW_from_OmM_kM_env(OmM_ij, kM_ij, k, tini=tini,
                                           tfin=tini + dtfin, FF=FF, A=A,
                                           multi=multi)

    return EGW, OmGW, k_br

def compute_PiM(EM, kp, dir='', alpha=2, alp_str='2'):

    import spectra as sp

    # read the reference Pi
    kf_Pi, Pi_ref = read_Pi(str_alp='2', dir=dir)
    # shift values of Pi by the factor kp
    Pi = shift_Pi(kf_Pi, Pi_ref, kp)
    # amplitude of Pi at k = 0
    Pi0 = value_Pi0_sbPL(alp=alpha)
    Pi = Pi0*Pi*EM**2/kp

    return kf_Pi, Pi

def OmGW_from_Pi_coh(t, k, Pi, NNt=100):

    Nt = len(t)
    Nk = len(k)
    OmGW = np.zeros((Nt, Nk))
    for j in range(0, Nt):
        tint = np.logspace(0, np.log10(t[j]), NNt)
        tij_int, kij = np.meshgrid(tint, k, indexing='ij')
        Pi_interp = np.zeros((NNt, Nk))
        for i in range(0, Nk):
            Pi_interp[:, i] = np.interp(tint, t, Pi[:, i])
        Pisqr = np.sqrt(Pi_interp)
        func_aux = Pisqr/tij_int*np.cos(kij*((t[j] - tij_int)))
        func = np.trapz(func_aux, tint, axis=0)
        OmGW[j, :] = func**2
    tij, kij = np.meshgrid(t, k, indexing='ij')
    EGW = 3*kij**2*OmGW
    OmGW = 3*kij**3*OmGW
    return EGW, OmGW

def dte_Om_k(OmM=.1, ks=2*np.pi, multi=False):

    """
    Function that returns the eddy turnover time as a function of the
    spectral peak and energy density (assumes radiation-dominated era to
    compute the Alfven speed).
    """

    if multi == True:
        ks, OmM = np.meshgrid(ks, OmM, indexing='ij')

    vA = np.sqrt(3/2*OmM)

    return 1/ks/vA

def dtfin_dte(dte):

    """
    Function that uses the empirical fit from numerical simulations to give
    a value of \delta \tfin as a function of the eddy turnover time.

    From A. Roper Pol et al., "The gravitational wave signal from primordial
    magnetic fields in the Pulsar Timing Array frequency band,"
    https://arxiv.org/abs/2201.05630 (2022).
    """
    dtfin = 0.184 + 1.937*dte
    return dtfin

def Af_dte(dte):

    """
    Function that uses the empirical fit from numerical simulations to give
    a value of the numerical ratio at the peak as a function of the eddy
    turnover time.

    From A. Roper Pol et al., "The gravitational wave signal from primordial
    magnetic fields in the Pulsar Timing Array frequency band,"
    https://arxiv.org/abs/2201.05630 (2022).
    """

    Af = 1.317 - .097*dte
    return Af

def slope_log_OmGW_model(k):

    """
    Function that returns the slope of the k^3 ln^2 (1 + 1/k) that appears
    in the analytical description of the GW spectra, see function_D2_env()
    and OmGW_from_Pi0().

    From A. Roper Pol et al., "The gravitational wave signal from primordial
    magnetic fields in the Pulsar Timing Array frequency band,"
    https://arxiv.org/abs/2201.05630 (2022).
    """

    beta = 3 - 2/(1 + k)/abs(np.log(1 + 1/k))

    return beta
