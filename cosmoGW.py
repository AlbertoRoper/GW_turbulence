"""
cosmoGW.py is a Python routine that contains functions relevant for
cosmological gravitational waves (GW).
"""

import astropy.constants as const
import astropy.units as u
import numpy as np

def shift_OmGW_today(k, OmGW, T, g):

    """
    Function that shifts the GW energy density spectrum from the time of
    generation to the present time.
    It assumes that the time of generation is within the radiation dominated
    era.

    Arguments:
        k -- array of wave numbers (normalized by the Hubble scale,
             computed from the Pencil Code)
        OmGW -- GW energy density spectrum per logarithmic interval
                (normalized by the radiation energy density, computed
                from the Pencil Code)
        T -- temperature scale at the time of generation
             (should have units of MeV or equivalent)
        g -- number of relativistic degrees of freedom at the time of
             generation

    Returns:
        f -- shifted wave number to frequency as a present time observable
        OmGW0 -- shifted OmGW spectrum to a present time observable
    """

    Hs_aux = np.sqrt(4*np.pi**3*const.G/45/const.c**5/const.hbar**3)
    Hs_aux = Hs_aux.to(u.Hz/u.MeV**2)*u.MeV**2
    T0 = 2.72548*u.K*const.k_B
    T0 = T0.to(u.MeV)
    g0 = 3.91
    as_aux = T0*g0**(1/3)/u.MeV
    H0_aux = 100*u.km/u.s/u.Mpc
    H0_aux = H0_aux.to(u.Hz)
    OmGW_aux = Hs_aux**2/H0_aux**2*as_aux**4
    f_aux = Hs_aux*as_aux/2/np.pi
    T = T.to(u.MeV)
    OmGW0 = OmGW*OmGW_aux*g**(-1/3)
    f = k*f_aux*g**(1/6)*T.value

    return f, OmGW0

def shift_hc_today(k, hc, T, g):

    """
    Function that shifts the characteristic amplitude spectrum from the time
    of generation to the present time.

    It assumes that the time of generation is within the radiation dominated
    era.

    Arguments:
        k -- array of wave numbers (normalized by the Hubble scale,
             computed from the Pencil Code)
        OmGW -- GW energy density spectrum per logarithmic interval
                (normalized by the radiation energy density, computed
                from the Pencil Code)
        T -- temperature scale at the time of generation
             (should have units of MeV or equivalent)
        g -- number of relativistic degrees of freedom at the time of
             generation

    Returns:
        f -- shifted wave number to frequency as a present time observable
        OmGW0 -- shifted OmGW spectrum to a present time observable
    """

    Hs_aux = np.sqrt(4*np.pi**3*const.G/45/const.c**5/const.hbar**3)
    Hs_aux = Hs_aux.to(u.Hz/u.MeV**2)*u.MeV**2
    T0 = 2.72548*u.K*const.k_B
    T0 = T0.to(u.MeV)
    g0 = 3.91
    as_aux = T0*g0**(1/3)/u.MeV
    H0_aux = 100*u.km/u.s/u.Mpc
    H0_aux = H0_aux.to(u.Hz)
    OmGW_aux = Hs_aux**2/H0_aux**2*as_aux**4
    f_aux = Hs_aux*as_aux/2/np.pi
    T = T.to(u.MeV)
    hc0 = hc*as_aux*g**(-1/3)/T.value
    f = k*f_aux*g**(1/6)*T.value

    return f, hc0

def fac_hc_OmGW(d=1):

    """
    Function that returns the factor to transform the strain function
    hc(f) to the GW energy density Omega_GW (f).

    Arguments:
        d -- option to give the factor to convert from energy density to
             strain if set to -1 (default 1)

    Returns:
        fac -- factor to convert from the strain function
               hc(f) to the GW energy density Omega_GW (f)
    """

    H0 = 100*u.km/u.s/u.Mpc
    H0 = H0.to(u.Hz)
    fac = H0.value*np.sqrt(3/2)/np.pi
    if d == -1: fac = 1/fac**2

    return fac

def hc_OmGW(f, OmGW, d=1):

    """
    Function that transforms the strain function hc(f) to the GW energy
    density Omega_GW (f).

    Arguments:
        f -- frequency array (given in Hz units)
        OmGW -- GW energy density spectrum Omega_GW (f)
        d -- option to convert from energy density to strain if set
             to -1 (default 1)

    Returns:
        hc -- strain spectrum
    """

    hc = fac_hc_OmGW()/f.value*np.sqrt(OmGW)
    if d==-1: hc = fac_hc_OmGW(d=-1)*f.value**2*OmGW**2

    return hc
