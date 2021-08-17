"""
spectra.py is a Python routine that contains postprocessing routines for
spectral functions.
"""

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

    import numpy as np

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

    import numpy as np

    k=k[np.where(k != 0)]
    E = abs(E[np.where(k != 0)])
    spec_mean = np.trapz(E, k)
    int = np.trapz(E*k**exp, k)
    # avoid zero division
    if exp >= 0 and spec_mean == 0: spec_mean = 1e-30
    if exp < 0 and int == 0: int = 1e-30
    kch = (int/spec_mean)**(1/exp)

    return kch

def min_max_stat(t, k, E, indt=0, plot=False):

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

    Returns:
        min_E -- maximum values of the spectral function over time
        max_E -- maximum values of the spectral function over time
        stat_E -- averaged values of the spectral function over time
                   from t[indt] to t[-1]
    """

    import numpy as np
    import matplotlib.pyplot as plt

    min_E = np.zeros(len(k)) + 1
    max_E = np.zeros(len(k))
    for i in range(indt, len(t)):
        if plot: plt.plot(k, E[i,:])
        min_E = np.minimum(E[i,:], min_E)
        max_E = np.maximum(E[i,:], max_E)
    # averaged spectrum
    stat_E = np.trapz(E[indt:,:], t[indt:], axis=0)/(t[-1] - t[indt])
    if plot:
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel('$k$')
        
    return min_E, max_E, stat_E
