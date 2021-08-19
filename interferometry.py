"""
interferometry.py is a Python routine that includes some functions
used to generate and read LISA and Taiji sensitivities to the detection
of stochastic gravitational wave backgrounds.
"""

def read_sens():
    """
    Function that reads the LISA and Taiji power law sensitivities
    (PLS), both Omega and Xi, previously generated and saved as pickle
    variables.

    Returns:
        f -- array of frequencies
        f_Xi_comb -- array of frequencies used for the PLS LISA_Taiji_XiPLS
        LISA_sensitivity -- LISA sensitivity expressed as a
                            GW energy density spectrum
        LISA_OmPLS -- LISA PLS to the GW energy density spectrum of a
                      cosmological GW background
        LISA_XiPLS -- LISA PLS to the helical GW energy density spectrum of a
                      cosmological GW background
        Taiji_OmPLS -- Taiji PLS to the GW energy density spectrum of a
                      cosmological GW background
        Taiji_XiPLS -- Taiji PLS to the helical GW energy density spectrum of a
                      cosmological GW background.
        LISA_Taiji_XiPLS -- PLS to the helical GW energy density spectrum of a
                            cosmological GW background obtained by combining
                            LISA and Taiji
    """

    import numpy as np

    dir = 'detector_sensitivity'
    f = read(dir, 'fs_PLS_Xi_LISA_SNR10.pckl')
    f_Xi_comb = read(dir, 'fs_PLS_Xi_LISA_Taiji_SNR10.pckl')
    LISA_sensitivity = read(dir, 'Om_LISA_sensitivity.pckl')
    LISA_OmPLS = read(dir, 'PLS_Om_LISA_SNR10.pckl')
    LISA_XiPLS = read(dir, 'PLS_Xi_LISA_SNR10.pckl')
    Taiji_OmPLS = read(dir, 'PLS_Om_Taiji_SNR10.pckl')
    Taiji_XiPLS = read(dir, 'PLS_Xi_Taiji_SNR10.pckl')
    LISA_Taiji_XiPLS = read(dir, 'PLS_Xi_LISA_Taiji_SNR10.pckl')

    return (f, f_Xi_comb, np.real(LISA_sensitivity), LISA_OmPLS, LISA_XiPLS,
            Taiji_OmPLS, Taiji_XiPLS, LISA_Taiji_XiPLS)

def read(dir, file):

    """
    Function that reads a pickle variable file and returns it.

    Arguments:
        dir -- directory that contains the file
        file -- name of the pickle file

    Returns:
        x -- pickle variable contained in the file
    """

    import pickle

    f = open(dir + '/' + file, 'rb')
    x = pickle.load(f)
    f.close()

    return x
