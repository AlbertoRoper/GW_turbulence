"""
cosmoGW.py is a Python routine that contains functions relevant for
cosmological sources of the stochastic gravitational wave background (SGWB).

Author: Alberto Roper Pol
"""

import astropy.constants as const
import astropy.units as u
import numpy as np

######################### Values at present time #########################

def values_0(h0=1.):

    """
    Function that returns the values of some parameters at the present time.

    Arguments:
        h0 -- parameterizes the uncertainties (Hubble tension) in the value
              of the Hubble rate (default 1)
    Returns:
        g0 -- 2 relativistic degrees of freedom (massive neutrinos)
        g0s -- 3.91 entropic/adiabatic degrees of freedom (including neutrinos)
        T0 -- temperature 2.72548 K, returned in energy units (MeV)
        H0 -- Hubble rate at the present time H0 = 100 h0 km/s/Mpc in frequency
              units (Hz)
    """

    g0 = 2
    g0s = 3.91
    T0 = 2.72548*u.K*const.k_B
    T0 = T0.to(u.MeV)
    H0 = 100*u.km/u.s/u.Mpc*h0
    H0 = H0.to(u.Hz)

    return g0, g0s, T0, H0

########################## SGWB at present time ##########################

def fac_hc_OmGW(d=1, h0=1.):

    """
    Function that returns the factor to transform the strain function
    hc(f) to the GW energy density OmGW (f) away from the source.

    Arguments:
        d -- option to give the factor to convert from energy density to
             strain if set to -1 (default 1)
        h0 -- parameterizes the uncertainties (Hubble tension) in the value
              of the Hubble rate (default 1); see values_0 function

    Returns:
        fac -- factor to convert from the strain function hc(f) to the GW
               energy density OmGW (f) in frequency units (Hz)
               
    Reference: M. Maggiore, "Gravitational wave experiments and early universe cosmology,"
    https://arxiv.org/pdf/gr-qc/9909001.pdf (2000); eq. 17.
    """

    g0, g0s, T0, H0 = values_0(h0=h0)
    fac = H0*np.sqrt(3/2)/np.pi
    if d == -1: fac = 1/fac**2

    return fac

def hc_OmGW(f, OmGW, d=1, h0=1.):

    """
    Function that transforms the  GW energy density OmGW (f) to the
    characteristc strain spectrum function hc(f) away from the source.

    Arguments:
        f -- frequency array (in units of frequency, e.g. Hz)
        OmGW -- GW energy density spectrum OmGW (f)
        d -- option to convert from energy density to strain if set
             to -1 (default 1)
        h0 -- parameterizes the uncertainties (Hubble tension) in the value
              of the Hubble rate (default 1); see values_0 function

    Returns:
        hc -- strain spectrum
        
    Reference: M. Maggiore, "Gravitational wave experiments and early universe cosmology,"
    https://arxiv.org/pdf/gr-qc/9909001.pdf (2000); eq. 17.
    """

    f = f.to(u.Hz)
    fac = fac_hc_OmGW(d=d, h0=h0)
    hc = fac/f*np.sqrt(OmGW)
    if d==-1: hc = fac*f**2*OmGW**2

    return hc

def hc_Sf(f, Sf, d=1):

    """
    Function that transforms the power spectral density Sf (f) to the
    characteristic strain spectrum function hc(f).

    Arguments:
        f -- frequency array (in units of frequency, e.g. Hz)
        Sf -- power spectral density Sf (f) (in units of 1/Hz^3)
        d -- option to convert from strain to power spectral density if set
             to -1 (default 1)

    Returns:
        hc -- strain spectrum
        
    Reference: A. Roper Pol, C. Caprini, A. Neronov, D. Semikoz, "The gravitational wave
    signal from primordial magnetic fields in the Pulsar Timing Array frequency band,"
    https://arxiv.org/pdf/2201.05630.pdf (2022); eq. 42.
    """

    f = f.to(u.Hz)
    Sf = Sf.to(1/u.Hz**3)
    hc = np.sqrt(12*np.pi**2)*np.sqrt(Sf*f**3)
    if d==-1: hc = Sf**2/12/np.pi**2/f**(3/2)

    return hc

########################  Radiation-dominated era ########################

def Hs_fact():

    """
    Function that computes the factor used in the calculation of the Hubble
    parameter at the time of generation (within the radiation-dominated era).

    The factor is given in units of Hz/MeV^2, such that the actual Hubble rate
    is given after multiplying fact*sqrt(g)*T^2, being g the number of dof and
    T the temperature scale (in MeV).
    
    Reference: A. Roper Pol, C. Caprini, A. Neronov, D. Semikoz, "The gravitational wave
    signal from primordial magnetic fields in the Pulsar Timing Array frequency band,"
    https://arxiv.org/pdf/2201.05630.pdf (2022); eq. 29.
    """

    fact = np.sqrt(4*np.pi**3*const.G/45/const.c**5/const.hbar**3)
    fact = fact.to(u.Hz/u.MeV**2)

    return fact

def Hs_val(g=10, T=100*u.MeV):

    """
    Function that computes the Hubble parameter at the time of generation
    (within the radiation-dominated era).

    Arguments:
        g -- number of relativistic degrees of freedom (dof) at the time of generation
             (default 10, i.e., ~QCD scale)
        T -- temperature scale at the time of generation in energy units
             (convertible to MeV) (default 100 MeV, i.e., ~QCD scale)
    Returns:
        Hs -- Hubble rate in frequency units (Hz)
        
    Reference: A. Roper Pol, C. Caprini, A. Neronov, D. Semikoz, "The gravitational wave
    signal from primordial magnetic fields in the Pulsar Timing Array frequency band,"
    https://arxiv.org/pdf/2201.05630.pdf (2022); eq. 29.
    """

    Hs_f = Hs_fact()
    T = T.to(u.MeV)
    Hs = Hs_f*T**2*np.sqrt(g)

    return Hs

def Om_rad(h0=1.):

    """
    Function that computes the ratio of radiation energy density to critical
    energy density at present time, Omrad^0.
    
    Arguments:
        h0 -- parameterizes the uncertainties (Hubble tension) in the value
              of the Hubble rate (default 1); see values_0 function
              
    Returns:
        Om_rad -- ratio of energy density given in natural units (MeV^4)
        
    #### TO REVIEW!! ####
    """

    Hs = Hs_fact()
    g0, g0s, T0, H0 = values_0(h0=h0)
    Om_rad = Hs**2/H0**2*T0**4*g0

    return Om_rad

def as_fact():

    """
    Function that computes the factor used in the calculation of the ratio
    between the scale factor at the time of generation and the present
    time assuming adiabatic expansion of the universe.

    The factor is in units of MeV and the ratio is obtained by multiplying
    fact*g^(-1/3)/T, being g the number of dof and T the temperature scale
    (in MeV).
    
    Reference: A. Roper Pol, C. Caprini, A. Neronov, D. Semikoz, "The gravitational wave
    signal from primordial magnetic fields in the Pulsar Timing Array frequency band,"
    https://arxiv.org/pdf/2201.05630.pdf (2022); eq. 28.
    """

    g0, g0s, T0, _ = values_0()
    fact = T0*g0s**(1/3)

    return fact

def as_a0_rat(g=10, T=100*u.MeV):

    """
    Function that computes the ratio between the scale factor at the time
    of generation (a*) and the present time (a_0) assuming adiabatic expansion of the
    universe.
    
    Arguments:
        g -- number of relativistic degrees of freedom (dof) at the time of generation
             (default 10, i.e., ~QCD scale)
        T -- temperature scale at the time of generation in energy units
             (convertible to MeV) (default 100 MeV, i.e., ~QCD scale)
    Returns:
        as_a0 -- ratio of scale factors (a*/a0)
        
    Reference: A. Roper Pol, C. Caprini, A. Neronov, D. Semikoz, "The gravitational wave
    signal from primordial magnetic fields in the Pulsar Timing Array frequency band,"
    https://arxiv.org/pdf/2201.05630.pdf (2022); eq. 28.
    """

    as_f = as_fact()
    T = T.to(u.MeV)
    as_a0 = as_f*g**(-1/3)/T

    return as_a0

def shift_onlyOmGW_today(OmGW, g=10, d=1, h0=1.):

    """
    Function that shifts the GW energy density spectrum from the time of
    generation to the present time.

    Arguments:
        OmGW -- GW energy density spectrum per logarithmic interval
                (normalized by the radiation energy density)
        g -- number of relativistic degrees of freedom (dof) at the time of generation
             (default 10, i.e., ~QCD scale)
        d -- option to reverse the transformation if set to -1, i.e.,
             to return OmGW(k) from the shifted to present time OmGW(f)
             (default 1)
        h0 -- parameterizes the uncertainties (Hubble tension) in the value
              of the Hubble rate (default 1); see values_0 function

    Returns:
        OmGW0 -- shifted spectrum OmGW to present time
        
    Reference: A. Roper Pol, C. Caprini, A. Neronov, D. Semikoz, "The gravitational wave
    signal from primordial magnetic fields in the Pulsar Timing Array frequency band,"
    https://arxiv.org/pdf/2201.05630.pdf (2022); eq. 27.
    """

    Hs_f = Hs_fact()*u.MeV**2
    as_f = as_fact()/u.MeV
    g0, g0s, T0, H0 = values_0(h0=h0)
    OmGW_f = Hs_f**2/H0**2*as_f**4
    OmGW0 = OmGW*OmGW_f*g**(-1/3)
    if d==-1: OmGW0 = OmGW/OmGW_f/g**(-1/3)

    return OmGW0

def shift_frequency_today(k, g=10, T=100*u.MeV, d=1):

    """
    Function that transforms the normalized wave number at the time of
    generation by the Hubble rate H_* to the present time frequency.

    Arguments:
        k -- array of wave numbers (normalized by the Hubble scale)
        g -- number of relativistic degrees of freedom (dof) at the time of generation
             (default 10, i.e., ~QCD scale)
        T -- temperature scale at the time of generation in energy units
             (convertible to MeV) (default 100 MeV, i.e., ~QCD scale)
        d -- option to reverse the transformation if set to -1, i.e.,
             to return the normalized k from the frequency shifted to present
             time f (default 1)

    Returns:
        f -- shifted wave number to frequency as a present time observable
             (in Hz)
             
    Reference: A. Roper Pol, A. Brandenburg, T. Kahniashvili, A. Kosowsky, S. Mandal,
    "The timestep constraint in solving the gravitational wave equations sourced by
    hydromagnetic turbulence," https://arxiv.org/pdf/1807.05479.pdf (2020); eq. B.13.
    """

    Hs_f = Hs_fact()
    as_f = as_fact()
    f_f = Hs_f*as_f/2/np.pi
    T = T.to(u.MeV)
    f = k*f_f*g**(1/6)*T
    if d==-1:
        k = k.to(u.Hz)
        f = k/f_f/g**(1/6)/T
        
    return f

def shift_OmGW_today(k, OmGW, g=10, T=100*u.MeV, d=1, h0=1.):

    """
    Function that shifts the GW energy density spectrum from the time of
    generation to the present time.
    It assumes that the time of generation is within the radiation dominated
    era.

    Arguments:
        k -- array of wave numbers (normalized by the Hubble scale)
        OmGW -- GW energy density spectrum per logarithmic interval
                (normalized by the radiation energy density)
        g -- number of relativistic degrees of freedom (dof) at the time of generation
             (default 10, i.e., ~QCD scale)
        T -- temperature scale at the time of generation in energy units
             (convertible to MeV) (default 100 MeV, i.e., ~QCD scale)
        d -- option to reverse the transformation if set to -1, i.e.,
             to return the normalized k from the frequency shifted to present
             time f (default 1)
        h0 -- parameterizes the uncertainties (Hubble tension) in the value
              of the Hubble rate (default 1); see values_0 function

    Returns:
        f -- shifted wave number to frequency as a present time observable
             (in Hz)
        OmGW0 -- shifted spectrum OmGW to present time
        
    Reference: see functions shift_onlyOmGW_today and shift_frequency_today
    """

    # shift Omega_GW
    OmGW0 = shift_onlyOmGW_today(OmGW, g=g, d=d, h0=h0)
    # shift frequency
    f = shift_frequency_today(k, g=g, T=T, d=d)
    
    return f, OmGW0

def shift_hc_today(k, hc, g=10, T=100*u.MeV, d=1):

    """
    Function that shifts the characteristic amplitude spectrum from the time
    of generation to the present time.

    It assumes that the time of generation is within the radiation dominated
    era.

    Arguments:
        k -- array of wave numbers (normalized by the Hubble scale)
        hc -- spectrum of GW characteristic amplitude per logarithmic interval
        g -- number of relativistic degrees of freedom (dof) at the time of generation
             (default 10, i.e., ~QCD scale)
        T -- temperature scale at the time of generation in energy units
             (convertible to MeV) (default 100 MeV, i.e., ~QCD scale)
        d -- option to reverse the transformation if set to -1, i.e.,
             to return the normalized k from the frequency shifted to present
             time f, and hc(k) from the shifted to present time hc0(f) (default 1)

    Returns:
        f -- shifted wave number to frequency as a present time observable
             (in Hz)
        hc0 -- shifted hc spectrum to present time
        
    Reference: A. Roper Pol, A. Brandenburg, T. Kahniashvili, A. Kosowsky, S. Mandal,
    "The timestep constraint in solving the gravitational wave equations sourced by
    hydromagnetic turbulence," https://arxiv.org/pdf/1807.05479.pdf (2020); eq. B.12.
    """

    as_f = as_fact()
    T = T.to(u.MeV)
    hc0 = hc*as_f*g**(-1/3)/T
    if d == -1: hc0 = hc/as_f/g**(-1/3)*T
    f = shift_frequency_today(k, g=g, T=T, d=d)

    return f, hc0
