"""
cosmoGW.py is a Python routine that contains functions relevant for
cosmological gravitational waves (GW).
"""

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

    import astropy.constants as const
    import astropy.units as u
    import numpy as np

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

    import astropy.constants as const
    import astropy.units as u
    import numpy as np

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
