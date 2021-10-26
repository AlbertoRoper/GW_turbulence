"""
cosmoGW.py is a Python routine that contains functions relevant for
cosmological gravitational waves (GW).
"""

import astropy.constants as const
import astropy.units as u
import numpy as np

def Hs_fact():

    """
    Function that computes the factor used in the calculation of the Hubble
    parameter at the time of generation (within the radiation-dominated era).

    The factor is given in units of Hz/MeV^2, such that the actual Hubble rate
    is given after multiplying fact*sqrt(g)*T^2
    """

    fact = np.sqrt(4*np.pi**3*const.G/45/const.c**5/const.hbar**3)
    fact = fact.to(u.Hz/u.MeV**2)

    return fact

def Hs_val(g, T):

    """
    Function that computes the Hubble parameter at the time of generation
    (within the radiation-dominated era).

    Arguments:
        g -- number of relativistic degrees of freedom
        T -- temperature scale in eV units
    Returns:
        Hs -- Hubble rate in Hz
    """

    Hs_f = Hs_fact()
    T = T.to(u.MeV)
    Hs = Hs_f*T**2*np.sqrt(g)

    return Hs

def H0_val(h0=1):

    """
    Function that computes the Hubble parameter rate at the present time
    as H0 = 100 h0 km/s/Mpc.

    Arguments:
        h0 -- parameterizes the uncertainties in the value
              of the Hubble rate (default 1), actual value is h0 ~ 0.7
    Returns:
        H0 -- Hubble rate in Hz
    """

    H0 = 100*u.km/u.s/u.Mpc
    H0 = H0.to(u.Hz)

    return H0

def as_fact():

    """
    Function that computes the factor used in the calculation of the ratio
    between the scale factor at the time of generation and the present
    time assuming adiabatic expansion of the universe.

    The factor is in units of MeV and the ratio is obtained by multiplying
    fact*g^(-1/3)/T.
    """

    T0 = 2.72548*u.K*const.k_B
    T0 = T0.to(u.MeV)
    g0 = 3.91
    fact = T0*g0**(1/3)

    return fact

def as_a0_rat(g, T):

    """
    Function that computes the ratio between the scale factor at the time
    of generation and the present time assuming adiabatic expansion of the
    universe.
    """

    as_f = as_fact()
    T = T.to(u.MeV)
    as_a0 = as_f*g**(-1/3)/T

    return as_a0

def shift_OmGW_today(k, OmGW, T, g, d=1):

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
        d -- option to reverse the transformation if set to -1, i.e.,
             to return the normalized k and OmGW(k) from the shifted to present
             time f, OmGW(f) (default 1)

    Returns:
        f -- shifted wave number to frequency as a present time observable
        OmGW0 -- shifted OmGW spectrum to a present time observable
    """

    Hs_f = Hs_fact()*u.MeV**2
    H0 = H0_val(h0=1.)
    as_f = as_fact()/u.MeV
    OmGW_f = Hs_f**2/H0**2*as_f**4
    f_f = Hs_f*as_f/2/np.pi
    T = T.to(u.MeV)
    OmGW0 = OmGW*OmGW_f*g**(-1/3)
    f = k*f_f*g**(1/6)*T.value
    if d==-1:
        f = k/f_f/g**(1/6)/T.value
        OmGW0 = OmGW/OmGW_f/g**(-1/3)

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

    Hs_f = Hs_fact()*u.MeV**2
    H0 = H0_val(h0=1.)
    as_f = as_fact()/u.MeV
    OmGW_f = Hs_f**2/H0**2*as_f**4
    f_f = Hs_f*as_f/2/np.pi
    T = T.to(u.MeV)
    hc0 = hc*as_f*g**(-1/3)/T.value
    f = k*f_f*g**(1/6)*T.value

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

    H0 = H0_val(h0=1.)
    fac = H0.value*np.sqrt(3/2)/np.pi
    if d == -1: fac = 1/fac**2

    return fac

def hc_OmGW(f, OmGW, d=1):

    """
    Function that transforms the  GW energy density Omega_GW (f) to the
    characteristc strain spectrum function hc(f).

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

def ks_infla(beta, gamma, eta=1):

    ks = 2*beta/(1 + eta)*(gamma + np.sqrt(1 + gamma**2 + 1/2/beta))

    return ks
