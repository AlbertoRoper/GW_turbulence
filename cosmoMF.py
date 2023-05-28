"""
cosmoMF.py is a Python routine that contains functions relevant for the
cosmological magnetic fields: bounds from different experiments, observations
or projected sensitivities, and expectations from theory, among others.

Author: Alberto Roper Pol
Date: 27/11/2021
"""

import astropy.constants as const
import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt
import plot_sets
import cosmology as co

# reference values and constants (EWPT)
Tref = 100*u.GeV # EWPT
gref = 100       # EWPT
cs2 = 1/3        # eos during RD era

# value used for recombination
Trec = 0.32*u.eV

##################### COMOVING COSMOLOGICAL MAGNETIC FIELDS #####################

def OmM_to_B_G(OmM, g=gref):

    """
    Function that transforms from magnetic energy density (as a ratio to the
    total energy density at the time of generation) to a comoving magnetic
    field strength in cgs units (Gauss).

    Arguments:
        OmM -- magnetic energy density fraction (from 0 to 1)
        g -- number of relativistic degrees of freedom (dof) at the time of generation
             (default is 100)
    Returns:
        BB -- strength of the comoving magnetic field in Gauss
        
    Reference: A. Roper Pol, C. Caprini, A. Neronov, D. Semikoz, "The gravitational wave
    signal from primordial magnetic fields in the Pulsar Timing Array frequency band,"
    Phys. Rev. D 105, 123502 (2022), arXiv:2201.05630, eq. 50
    """

    mu0 = 40*np.pi/u.J*u.m**3               # mu0 to convert from J/m^3 to G^2
    # compute radiation energy density at the time of generation
    # (in MeV, to be cancelled later)
    EErad = co.rho_radiation(T=1*u.MeV, g=g)
    EEM = OmM*EErad                         # magnetic energy density
    BB = np.sqrt(2*mu0*EEM)                 # magnetic field strength
    # factor from expansion a^2
    BB = BB.to(1)*co.as_a0_rat(T=1*u.MeV, g=g)**2*u.G

    return BB

def k_to_ll_Mpc(kf, g=gref, T=Tref, kk=True):

    """
    Function that converts the characteristic wave number k_*/HH_*
    (normalized with the conformal Hubble rate H_*) to the comoving
    length scale in Megaparsec (Mpc).

    Arguments:
        kf -- characteristic scale (minimum is 2pi, which is the Hubble wave
              number)
        g -- number of relativistic degrees of freedom (dof) at the time of generation
             (default is 100)
        T -- temperature scale at the time of generation in energy units
             (convertible to MeV) (default is 100 GeV)
        kk -- if kk is True, then kf corresponds to k_* HH_*, otherwise it refers
              to the length in terms of the Hubble size HH_* l_* = 2 \pi/(k_*/HH_*)
              
    Returns:
        ll -- comoving length scale of the magnetic field in Mpc
        
    Reference: A. Roper Pol, C. Caprini, A. Neronov, D. Semikoz, "The gravitational wave
    signal from primordial magnetic fields in the Pulsar Timing Array frequency band,"
    Phys. Rev. D 105, 123502 (2022), arXiv:2201.05630, eq. 51
    """

    HHs = co.Hs_val(g=g, T=T)*co.as_a0_rat(g=g, T=T)
    
    if kk: ll = 2*np.pi/kf*const.c/HHs
    else: ll = const.c/HHs*kf
    ll = ll.to(u.Mpc)

    return ll

def Alfven_velocity(OmM, cs2=1/3):
    
    """
    Function that computes the Alfvén velocity, following A. Brandenburg, K. Enqvist,
    P. Olesen, "Large scale magnetic fields from hydromagnetic turbulence in the very
    early universe," Phys.Rev.D 54, 1291 (1996).
    
    Reference: A. Roper Pol, C. Caprini, A. Neronov, D. Semikoz, "The gravitational wave
    signal from primordial magnetic fields in the Pulsar Timing Array frequency band,"
    Phys. Rev. D 105, 123502 (2022), arXiv:2201.05630, footnote 4
    
    Arguments:
        OmM -- magnetic energy density fraction (\rho_M/\rho_rad)
        cs2 -- equation of state p = \rho cs2 (default is radiation domination, i.e. 1/3)
    """
    
    vA = np.sqrt(2*OmM/(1 + cs2))
    
    return vA

def B_vs_ll_LPE(ll, T=Tref, cs2=cs2, d=1, A=1, exp=1):

    """
    Function that computes the comoving magnetic field strength (in Gauss) as
    a function of the comoving coherence length (in Mpc) in the case that the
    scale is that of the largest processed eddies of the MHD turbulence
    dominated by magnetic energy.
    
    Reference: A. Roper Pol, C. Caprini, A. Neronov, D. Semikoz, "The gravitational wave
    signal from primordial magnetic fields in the Pulsar Timing Array frequency band,"
    Phys. Rev. D 105, 123502 (2022), arXiv:2201.05630, eq. 52
    
    01/05/23, Alberto Roper Pol:
    
    added option to use reconnection time instead of Alfven time to relate B vs \lambda,
    following D. Hosking, A. Schekochihin, "Cosmic-void observations reconciled with
    primordial magnetogenesis," arXiv:2203.03573 (2022),
    
    used in A. Roper Pol, T. Boyer, C. Caprini, A. Neronov, D. Semikoz, "LISA and CTA
    multi-messenger probe of first-order cosmological phase transitions,"
    arXiv:2206.XXXX (2023)
    
    It takes that the time scale at H_* t_* ~ 1 (valid assuming RD) to be
    a generic:
    
    H_* t_* ~ 1 ~ (A H_* l_*)^exp/vA (not how it is defined below yet, to be corrected)
    
    where A = exp = 1 recovers the Alfvenic speed, vA = \sqrt(2/(1 + cs2) OmM) = \sqrt(3/2 OmM).
    OmM and H_rec l_rec = H_* l_* are constant during RD era.
    
    We use OmM_to_B_G to relate B with OmM and k_to_ll_Mpc to relate l_* H_* with l
    
    Arguments:
        ll -- comoving length scale of the magnetic field in Mpc
        T -- temperature scale at which B is computed
        cs2 -- equation of state p = \rho cs2 (default is radiation domination, i.e. 1/3)
        d -- option to return ll vs B instead of B vs ll if d = -1 (default is d = 1)
        A -- prefactor used to characterize different options from reconnection time
        exp -- exponent used to characterize different options from reconnection time
        
    Returns:
        B -- strength of the comoving magnetic field in Gauss
    """

    # final result does not depend on g*
    BB0 = OmM_to_B_G(1., g=gref)                 
    ll0 = k_to_ll_Mpc(1., T=T, g=gref, kk=False) # equivalent to c/H_*
    vA0 = Alfven_velocity(1., cs2=cs2)           # prefactor from vA vs OmM

    if d == 1:
        ll = ll.to(u.Mpc)
        B = (BB0/ll0/vA0*ll*A)**exp
    if d == -1:
        ll = ll.to(u.G)
        B = (ll0/BB0*vA0*ll/A)**(1/exp)

    return B

def B_vs_ll_rec(ll, T=Tref, cs2=cs2, d=1, A=1, exp=1):

    """
    Function that computes the comoving magnetic field strength (in Gauss) at
    the epoch of recombination as a function of the comoving coherence length
    (in Mpc), assuming that the scale is that of the largest processed eddies
    developed or that determined by the reconnection time, 
    following MHD turbulence decay dominated by magnetic energy.
    
    It uses function B_vs_ll_LPE setting T = Trec
    
    Reference: A. Roper Pol, C. Caprini, A. Neronov, D. Semikoz, "The gravitational wave
    signal from primordial magnetic fields in the Pulsar Timing Array frequency band,"
    Phys. Rev. D 105, 123502 (2022), arXiv:2201.05630, eq. 53
    
     Arguments:
        ll -- comoving length scale of the magnetic field in Mpc
        cs2 -- equation of state p = \rho cs2 (default is radiation domination, i.e. 1/3)
        d -- option to return ll vs B instead of B vs ll if d = -1 (default is d = 1)
        A -- prefactor used to characterize different options from reconnection time
        exp -- exponent used to characterize different options from reconnection time
        
    Returns:
        B -- strength of the comoving magnetic field at recombination in Gauss
    """
    
    B = B_vs_ll_LPE(ll, T=Trec, cs2=cs2, d=1, A=1, exp=1)

    return B

######################### LIMITS FROM IGMF OBSERVATIONS #########################

def gr_limits(x, exp='Fermi_timing'):

    """
    Function that returns the lower limits on the primordial magnetic field
    amplitude in G (y1 and y2) derived from different observations as a function
    of the coherence length in Mpc (x).

    The reference is 

    Arguments:
        x -- coherence length in Mpc
        exp -- reference experiment for the limits:
        
            - 'Fermi': From [Fermi-LAT Collaboration], "Search for Spatial Extension
                       in High-Latitude Sources Detected by the Fermi Large Area Telescope,"
                       Astrophys. J. Suppl. 237, 32 (2018), arXiv:1804.08035
                       
                  - 'Fermi_timing': from timing of the blazar signal
                  - 'Fermi_extended': from the search of extended emission
                  
            - 'CTA': expected sensitivity of CTA, from A. Korochkin, O. Kalashev, A. Neronov, and
                     D. Semikoz, "Sensitivity reach of gamma-ray measurements for strong
                     cosmological magnetic fields," Astrophys. J. 906, 116 (2021),
                     arXiv:2007.14331
                     
            - 'Faraday': upper limit from measures of Faraday rotation, from M. S. Pshirkov,
                         P. G. Tinyakov, and F. R. Urban, "New limits on extragalactic magnetic
                         fields from rotation measures," Phys. Rev. Lett. 116, 191302 (2016),
                         arXiv:1504.06546
            - 'UHECR': observations from ultra-high energy cosmic rays (UHECR), from
                       A. Neronov, D. Semikoz, and O. Kalashev, "Limit on intergalactic magnetic
                       field from ultra-high-energy cosmic ray hotspot in Perseus-Pisces region,"
                       arXiv:2112.08202 (2021)

    Returns:
        y -- amplitudes of the magnetic field in Gauss
    """

    x = x.to(u.Mpc)
    a = -.5                   # slope
    
    if 'Fermi_timing':
        A = 1e-16
        x_br = 3e-2
    if 'Fermi_extended':
        A = 1e-15
        x_br = 3e-2
    if 'CTA':
        A = 1.5e-12
        x_br = 3e-2
    if 'Faraday':
        A = 1e-9
        x_br = 1
    if 'UHECR':
        A = 1e-10
        x_br = 72.68
        
    A = A*u.G
    x_br = x_br*u.Mpc
    y = A*(x/x_br)**a
    
    # position of break to flat sensitivity
    if np.size(x) == 1:
        if x > x_br: y = A
    else: y[(x > x_br)] = A

    return y

def CMB_Galli():

    """
    Function that returns the upper limit on post-recombination magnetic
    field strength obtained from CMB experiments and the lower limit obtained
    assuming a SHOES-based posterior distribution to alleviate the Hubble
    tension.

    The reference is S. Galli, L. Pogosian, K. Jedamzik, and L. Balkenhol,
    "Consistency of Planck, ACT and SPT constraints on magnetically assisted
    recombination and forecasts for future experiments," Phys.Rev.D 105 2 (2022),
    023513, arXiv:2109.03816.
    """

    B_up = 1.04e-10*u.G
    B_bot = .76e-10*u.G

    return B_up, B_bot
