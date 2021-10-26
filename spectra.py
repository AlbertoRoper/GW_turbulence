"""
spectra.py is a Python routine that contains postprocessing routines for
spectral functions.
"""

import numpy as np

def compute_kpeak(k, E, tol=.01, quiet=False):

    """
    Function that computes the maximum of the spectrum E and its spectral
    peak.

    Arguments:
        k -- array of wave numbers
        E -- array of the spectral values
        tol -- factor to avoid faulty maxima due to nearly flat spectrum
               (default 1%)
        quiet -- option to print out the result if quiet is False
                 (default False)

    Return:
        kpeak -- position of the spectral peak
        Emax -- maximum value of the spectrum
    """

    max1 = np.argmax(E)
    indmax = max1

    if E[max1] == 0:
        Emax = 0
        kpeak = 0
    else:
        max2 = np.argmax(k*E)
        # if the maximum of the spectrum is within tol of the maximum of k*E,
        # then we take as the maximum value where k*E is maximum, to take into
        # account flat and nearly flat spectra
        if abs(E[max1] - E[max2])/E[max1] < tol: indmax = max2
        Emax = E[indmax]
        kpeak = k[indmax]

    if not quiet:
        print('The maximum value of the spectrum is ', Emax,
              ' and the spectral peak is ', kpeak)

    return kpeak, Emax

def max_E_kf(k, E, exp=0):

    """
    Function that computes the maximum of a spectrum compensated by the
    wave number, i.e., max(E*k^exp)

    Arguments:
        k -- array of wave numbers
        E -- array of spectral values
        exp -- exponent of k (default 0)
    """

    indmax = np.argmax(k**exp*E)
    max_k = k[indmax]
    max_E = E[indmax]

    return max_k, max_E

def characteristic_k(k, E, exp=1):

    """
    Function that computes the characteristic wave number.

    Arguments:
        k -- array of wave numbers
        E -- array of spectral values
        exp -- exponent used to define the characteristic wave number
               k_ch ~ (\int k^exp E dk/\int E dk)^(1/exp)
               (default 1)

    Returns:
        kch -- characteristic wave number defined with the power 'exp'
    """

    k = k[np.where(k != 0)]
    E = abs(E[np.where(k != 0)])
    spec_mean = np.trapz(E, k)
    int = np.trapz(E*k**exp, k)
    # avoid zero division
    if exp >= 0 and spec_mean == 0: spec_mean = 1e-30
    if exp < 0 and int == 0: int = 1e-30
    kch = (int/spec_mean)**(1/exp)

    return kch

def min_max_stat(t, k, E, indt=0, plot=False, hel=False):

    """
    Function that computes the minimum, the maximum, and the averaged
    functions over time of a spectral function.

    Arguments:
        t -- time array
        k -- wave number array
        E -- spectrum 2d array (first index t, second index k)
        indt -- index of time array to perform the average
                from t[indt] to t[-1]
        plot -- option to overplot all spectral functions
                from t[indt] to t[-1]
        hel -- option for helical spectral functions where positive and
               negative values can appear (default False)
               It then returns min_E_pos, min_E_neg, max_E_pos, max_E_neg
               referring to the maximum/minimum absolute values of the positive
               and negative values of the helical funtion.

    Returns:
        min_E -- maximum values of the spectral function over time
        max_E -- maximum values of the spectral function over time
        stat_E -- averaged values of the spectral function over time
                   from t[indt] to t[-1]
    """

    import matplotlib.pyplot as plt

    if hel:
        min_E_neg = np.zeros(len(k)) + 1e30
        max_E_neg = np.zeros(len(k))
        min_E_pos = np.zeros(len(k)) + 1e30
        max_E_pos = np.zeros(len(k))
    else:
        min_E = np.zeros(len(k)) + 1e30
        max_E = np.zeros(len(k))
    for i in range(indt, len(t)):
        if hel:
            if plot: plt.plot(k, abs(E[i,:]))
            # split between positive and negative values
            x_pos, x_neg, f_pos, f_neg, color = red_blue_func(k, E[i, :])
            for j in range(0, len(k)):
                if k[j] in x_pos:
                    indx = np.where(x_pos == k[j])[0][0]
                    min_E_pos[j] = min(min_E_pos[j], f_pos[indx])
                    max_E_pos[j] = max(max_E_pos[j], f_pos[indx])
                else:
                    indx = np.where(x_neg == k[j])[0][0]
                    min_E_neg[j] = min(min_E_neg[j], abs(f_neg[indx]))
                    max_E_neg[j] = max(max_E_neg[j], abs(f_neg[indx]))
        else:
            E = abs(E)
            if plot: plt.plot(k, E[i,:])
            min_E = np.minimum(E[i,:], min_E)
            max_E = np.maximum(E[i,:], max_E)
    # averaged spectrum
    stat_E = np.trapz(E[indt:,:], t[indt:], axis=0)/(t[-1] - t[indt])
    if hel:
        min_E_pos[np.where(min_E_pos == 1e30)] = \
                abs(stat_E[np.where(min_E_pos == 1e30)])
        min_E_neg[np.where(min_E_neg == 1e30)] = \
                abs(stat_E[np.where(min_E_neg == 1e30)])
        max_E_pos[np.where(max_E_pos == 0)] = \
                abs(stat_E[np.where(max_E_pos == 0)])
        max_E_neg[np.where(max_E_neg == 0)] = \
                abs(stat_E[np.where(max_E_neg == 0)])
    else:
        min_E[np.where(min_E == 1e30)] = \
                abs(stat_E[np.where(min_E == 1e30)])
        max_E[np.where(min_E == 0)] = \
                abs(stat_E[np.where(min_E == 0)])
    if plot:
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel('$k$')

    if hel: return min_E_pos, min_E_neg, max_E_pos, max_E_neg, stat_E
    else: return min_E, max_E, stat_E

def local_max(k, E, order=1):

    """
    Function that computes the local maxima of a spectrum.

    Arguments:
        k -- array of wave numbers
        E -- spectrum E
        order -- order of the local maximum solver, which uses
                 scipy.signal.argrelextrema

    Returns:
        kmax -- position of the local maxima
        Emax -- values of the local maxima
    """

    from scipy.signal import argrelextrema

    inds_model_max = argrelextrema(E,
                 np.greater, order=order)
    kmax = k[inds_model_max]
    Emax = E[inds_model_max]

    return kmax, Emax

def compute_yks(k, E, N):

    """
    Function that computes N power law fittings, logarithmically
    equidistributed in k, of the spectrum E.

    Arguments:
        k -- array of wave numbers
        E -- array of spectrum values
        N -- number of power law fittings to discretize E
    """

    # compute N number of single power law fits around the model
    kps = np.logspace(np.log10(k[0]),
                      np.log10(k[-1]), N + 2)
    Ess = np.interp(kps, k, E)
    akss = np.zeros(N)
    akss[0], c = slope(Ess[1], Ess[2],
                       kps[1], kps[2])
    kss = np.logspace(np.log10(kps[0]),
                      np.log10(kps[2]), 5)
    Ekss = kss**akss[0]*10**c
    for i in range(2, N+1):
        akss[i - 1], c = slope(Ess[i], Ess[i + 1],
                               kps[i], kps[i + 1])
        ksss = np.logspace(np.log10(kps[i]),
                           np.log10(kps[i + 1]), 5)
        Eksss = ksss**akss[i - 1]*10**c
        kss = np.append(kss, ksss)
        Ekss = np.append(Ekss, Eksss)
    km, Em = mean_pos_loglog(np.append(kps[0], kps[2:]),
                             np.append(Ess[0], Ess[2:]))

    return kss, Ekss, akss, km, Em

def mean_pos_loglog(k, E):

    """
    Function that computes the loglog middle values km, EM of the intervals
    of the arrays k, E

    Arguments:
        k -- array of wave numbers
        E -- array of spectrum values

    Returns:
        km -- array of middle log values of the k intervals
        Em -- array of middle log values of the E intervals
    """

    N = len(k)
    km = np.zeros(N - 1)
    Em = np.zeros(N - 1)
    for i in range(1, N):
        km[i-1] = np.sqrt(k[i-1]*k[i])
        Em[i-1] = np.sqrt(E[i-1]*E[i])

    return km, Em

def slope(y1, y2, x1, x2):

    """
    Function that computes the slope between points 1 and 2

    Arguments:
        x1 -- x coordinate of point 1
        x2 -- x coordinate of point 2
        y1 -- y coordinate of point 1
        y2 -- y coordinate of point 2

    Returns:
        a -- slope between points 1 and 2
        c -- y-intercept of the straight line joining points 1 and 2
    """

    a = np.log10(y1/y2)/np.log10(x1/x2)
    c = np.log10(y1) - a*np.log10(x1)

    return a, c

def red_blue_func(x, f, col=0):

    """
    Function that splits an array into positive and negative values, and
    assigns colors (red to positive and blue to negative).

    Arguments:
        x -- array of x
        f -- array of the function values
        col -- option to choose blue and red (default 0 is red for positive
               and blue for negative, 1 is swapped)

    Returns:
        x_pos -- array of x values where f is positive
        x_neg -- array of x values where f is negative
        f_pos -- array of f values where f is positive
        f_neg -- array of f values where f is negative
        color -- array of colors assigned (blue and red)
    """

    N = len(f)
    color = []
    f_pos=[]; x_pos=[]
    f_neg=[]; x_neg=[]
    for i in range(0, N):
        sgn = np.sign(f[i])
        if sgn > 0:
            if col == 0: color.append('red')
            if col == 1: color.append('blue')
            f_pos.append(f[i])
            x_pos.append(x[i])
        else:
            if col == 0: color.append('blue')
            if col == 1: color.append('red')
            f_neg.append(f[i])
            x_neg.append(x[i])
    f_pos = np.array(f_pos)
    f_neg = np.array(f_neg)
    x_pos = np.array(x_pos)
    x_neg = np.array(x_neg)

    return x_pos, x_neg, f_pos, f_neg, color

def plot_neg_pos(x, f, ls1='solid', lw1=1, ls2=':', lw2=2, col='black'):

    """
    Function that splits an array into positive and negative values, and
    plots them with different line styles.

    Arguments:
        x -- array of x
        f -- array of the function values
        col -- option to choose blue and red (default 0 is red for positive
               and blue for negative, 1 is swapped)
    """

    import matplotlib.pyplot as plt

    # plot positive and negative values with different line styles
    sgn = np.sign(f)
    converge = False
    sgn0 = sgn[0]
    i = 0
    lw = 1
    while not converge:
        sign = False
        i0 = i
        while not sign and not converge:
            if sgn0 == 1:
                ls = ls1
                lw = lw1
            else:
                ls = ls2
                lw = lw2
            if i==len(sgn) - 2: converge=True
            if sgn[i] != sgn0:
                sign = True
                sgn0 = sgn[i]
            i += 1
        plt.plot(x[i0-1:i], abs(f[i0-1:i]),
                 color=col, ls=ls, lw=lw)

def str_exp(exp, ak, den, diff=0.05):

    """
    Function that returns a string k^(a/den) if the absolute difference between
    the value a/den and the exponent ak is below diff.

    Arguments:
        exp -- initial string (given by the previous best estimation)
        ak -- slope
        den -- denominator of the fractions to be tested
        diff -- difference used to accept the corrected fraction approximating
                the slope

    Returns:
        exp -- updated string of the fractional slope
        diff -- difference updated
    """

    test = np.array(range(1, int(20*den)))
    test_i = test/den
    difft = abs(test_i - abs(ak))
    ind_min = np.argmin(difft)
    if difft[ind_min] < diff:
        m = test[ind_min]
        if ak > 0:
            if den == 1: exp = '$k^{%i}$'%m
            else: exp = '$k^{%i/%i}$'%(m, den)
        else:
            if den == 1: exp = '$k^{-%i}$'%m
            else: exp = '$k^{-%i/%i}$'%(m, den)
        diff = difft[ind_min]

    return exp, diff
