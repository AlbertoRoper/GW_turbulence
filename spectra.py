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
