"""
cosmology.py is a Python routine that contains functions relevant
for cosmological calculations, including a solver to Friedmann equations.
Author: Alberto Roper Pol
Date: 27/11/2022

The code has been developed and used for the results of Y. He, A. Roper Pol,
A. Brandenburg, "Modified propagation of gravitational waves from the early
radiation era," submitted to JCAP (2022).
"""

import astropy.constants as const
import astropy.units as u
import pandas as pd
import numpy as np
import spectra as sp

######################### Values at present time #########################

def values_0(h0=1., neut=False, Neff=3, ret_rad=False):

    """
    Function that returns the values of some parameters at the present time.
    
    Arguments:
        h0 -- parameterizes the uncertainties (Hubble tension) in the value
              of the Hubble rate (default 1)
        neut -- option to add neutrinos in the calculation of today's relativistic
                dofs
        Neff -- effective number of neutrino species (default is 3)
        ret_rad -- option to return radiation energy density at present time
        
    Returns:
        g0 -- 2 relativistic degrees of freedom (massive neutrinos)
        g0s -- 3.91 entropic/adiabatic degrees of freedom (including neutrinos)
        T0 -- temperature 2.72548 K, returned in energy units (MeV)
        H0 -- Hubble rate at the present time H0 = 100 h0 km/s/Mpc in frequency
              units (Hz)
        rho_rad0, Om_rad0 -- radiation energy density (and ratio) at present time
        
    Reference: Notes on Cosmology by Syksy Räsänen, chapter 4:
    'Thermal history of the early universe'
    (http://www.courses.physics.helsinki.fi/teor/cos1/cosmo2015_05.pdf).
    """

    # relativistic dofs from photons
    g0 = 2
    # contribution from neutrinos if neut is True
    if neut: g0 = 2*(1 + Neff*7/8*(4/11)**(4/3))
    # adiabatic dofs
    g0s = 2*(1 + Neff*7/8*4/11)
    T0 = 2.72548*u.K*const.k_B
    T0 = T0.to(u.MeV)
    H0 = 100*u.km/u.s/u.Mpc*h0
    H0 = H0.to(u.Hz)
    
    if ret_rad:
        rho_rad0 = rho_radiation(T=T0, g=g0)
        rho_c = rho_critical(H0)
        Om_rad0 = rho_rad0/rho_c
        return g0, g0s, T0, H0, rho_rad0, Om_rad0
    else:
        return g0, g0s, T0, H0

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
        g -- number of adiabatic degrees of freedom (dof) at the time of generation
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

def rho_radiation(T=100*u.GeV, g=100):
    
    """
    Function that computes the radiation energy density at different epochs of the
    universe.
    
    Arguments:
        T -- temperature in energy units (default is EWPT ~ 100 GeV)
        g -- relativistic dofs (default is EWPT ~ 100)
        
    Returns:
        rho_rad -- energy density in units GeV/m^3
    """
    
    rho_rad = np.pi**2/30*g*T**4/(const.hbar*const.c)**3
    rho_rad = rho_rad.to(u.GeV/u.m**3)
    
    return rho_rad

def rho_critical(H):
    
    """
    Function that computes the critical energy density at different epochs of the
    universe.
    
    Arguments:
        H -- Hubble rate in units of frequency
        
    Returns:
        rho_c -- energy density in units GeV/m^3
    """
    
    rho_c = 3*H**2*const.c**2/8/np.pi/const.G
    rho_c = rho_c.to(u.GeV/u.m**3)
    
    return rho_c

############################ COSMOLOGY CALCULATIONS ############################

def thermal_g(dir0='', T=100*u.MeV, s=0, file=True):

    """
    Returns the relativistic dof g_* as a function of T_* according to the
    thermal history of the Universe. Note that for T > 0.5 MeV, after neutrino
    decoupling, entropic and relativistic g are approximately equal.

    Arguments:
        T -- temperature given in enery units (convertible to MeV)
             (default 100 MeV, i.e., ~QCD scale)
        s -- option to return adiabatic (s=1) dof instead of relativistic
             (default 0)
        file -- option to read g_* or gS from a file with numerical values,
                based on numerical calculations (default False)

    Returns:
        g -- relativistic degrees of freedom

    Reference: Notes on Cosmology by Syksy Räsänen, chapter 5:
    'Thermal history of the early universe'
    (http://www.courses.physics.helsinki.fi/teor/cos1/cosmo2015_05.pdf);
    see table 3 and figure 1 (stored in file).
    """
    
    g0, g0s, T0, H0 = values_0(neut=True)
    # make sure that the values at present time correspond to those at the
    # end of RD (taking relativistic neutrinos)

    if file:
        import pandas as pd
        df = pd.read_csv(dir0 + 'cosmology/T_gs.csv')
        Ts = np.array(df['T [GeV]'])
        if s == 0: gs = np.array(df['g_*'])
        if s == 1: gs = np.array(df['gS'])
        T = T.to(u.GeV)    # values of T from file are in GeV
        g = np.interp(T.value, np.sort(Ts), np.sort(gs))

    else:
        T = T.to(u.MeV)
        T = T.value
        # Check value of T in MeV and assign g_*
        if T < 0.1:
            if s == 0: g = g0
            if s == 1: g = g0s
        elif T < 0.5: g = 7.25
        elif T <= 100: g = 10.75
        elif T <= 150: g = 17.25
        elif T < 1e3: g = 61.75
        elif T < 4e3: g = 75.75
        elif T < 8e4: g = 86.25
        elif T < 1.7e5: g = 96.25
        else: g = 106.75

    return g

############################### FRIEDMANN EQUATIONS ###############################

def RD_dofs(dir0='', Neff=3.):
    
    """
    Function that computes the degrees of freedom (relativistic and adiabatic)
    during the RD era.
    """
    
    #### compute the relativistic degrees of freedom as a function of T during the RD epoch
    #### from a stored file using the function thermal_g

    T = np.logspace(-4, 8, 1000)*u.MeV
    gS = thermal_g(T=T, s=1, file=True, dir0=dir0)
    gs = thermal_g(T=T, s=0, file=True, dir0=dir0)
    
    g0, g0s, T0, H0 = values_0(neut=True, Neff=Neff)
    
    #### the numerical gs and gS have final (smaller) values of 3.363 and 3.909,
    #### which are not necessary the same as g0 and gS0, especially if we set Neff
    #### different than 3, so we interpolate the last values to correct for this.
    if gs[0] < g0:
        inds = np.where(gs < max(gs[0], g0))[0]
        gs[:inds[-1]+1] = g0
    else:
        inds = np.where((gs - gs[0] < 1e-1))[0]
        indd = 10
        vals = np.linspace(g0, gs[inds[-1]], indd)
        gs[:inds[-1]] = g0
        gs[inds[-1]-indd:inds[-1]] = vals

    if gS[0] < g0s:
        inds = np.where((gS < max(gS[0], g0s)))[0]
        gS[:inds[-1]+1] = g0s
    else:
        inds = np.where(gS - gS[0] < 1e-1)[0]
        indd = 10
        vals = np.linspace(g0s, gS[inds[-1]], indd)
        gS[:inds[-1]] = g0s
        gS[inds[-1]-indd:inds[-1]] = vals

    ### convert temperature and degrees of freedom into an array of scale factors using
    ### adiabatic expansion of the universe
    as_T = (T0.to(u.MeV)/T)*(g0s/gs)**(1/3)
    
    return T, as_T, gs, gS

def dofs_vs_a(a, dir0='', Neff=3.):
    
    """
    Function that computes the degrees of freedom (relativistic and adiabatic)
    for an array of scale factors.
    """
    
    #### Read the dofs during RD era
    T, as_T, gs, gS = RD_dofs(dir0=dir0, Neff=Neff)
    
    #### convert arrays of degrees of freedom and interpolate to original array of scale factors a
    #### (needs to be sorted first for interp to work)
    inds = np.argsort(as_T)
    as2 = as_T[inds]
    gs = gs[inds]
    gS = gS[inds]
    gs = np.interp(a, as2, gs)
    gS = np.interp(a, as2, gS)
    
    return gs, gS

def Omega_rad_dof(a, dir0='', Neff=3.):
    
    """
    Function that computes the factor that takes into account the radiation
    energy density ratio variation due to the inclusion of varying dofs during
    the RD era.
    """
    
    gs, gS = dofs_vs_a(a, dir0=dir0, Neff=Neff)
    g0, g0s, T0, H0 = values_0(neut=True, Neff=Neff)
    Om_rat_dof = (gs/g0)*(gS/g0s)**(-4/3)

    #### the numerical results for the factor (gs/g0)*(gS/g0s)**(-4/3)
    #### present some numerical inacuracy, so it requires some smoothing
    amax, Emax = sp.local_max(a[a>1e-10], Om_rat_dof[a>1e-10])
    Emax[Emax > 1] = 1.
    aa = np.append(a[a<=1e-10], amax)
    Om_rat_dof = np.append(Om_rat_dof[a<=1e-10], Emax)
    Om_rat_dof = np.interp(a, aa, Om_rat_dof)
    
    return Om_rat_dof
    
def Omega_vs_a(a, a0=1, h0=0.6732, OmL0=0.6841, dir0='', dofs=True, Neff=3.):
    
    """
    Function that computes the energy density ratio to present-time critical energy
    density as a function of the scale factor a, for a universe composed by matter,
    radiation, and dark energy:
    
    \Omega (a) = OmR0 x a^(-4) + OmM0 x a^(-3) + OmL0
    
    Arguments:
        a -- array of scale factors
        a0 -- reference value of the scale factor at present time (default is 1)
        h0 -- present-time value of the Hubble rate H0 = h0 x 100 km/s/Mpc
              (default is 67.32 km/s/Mpc based on CMB observations)
        OmL0 -- present-time content of dark energy
        dofs -- option to compensate the rad energy density using dofs during RD era
        dir0 -- directory where file for dofs during RD era is stored
        
    Returns:
        Om_tot -- total energy density (normalized to present-time critical)
        Om_rad -- radiation energy density (normalized)
        Om_matt -- matter energy density (normalized)
    
    Reference: Y. He, A. Roper Pol, A. Brandenburg, "Modified propagation of
    gravitational waves from the early radiation era," in preparation.
    """
    
    ########## compute Om_rad0
    g0, g0s, T0, H0, rho_rad0, OmR0 = values_0(h0=h0, neut=True,
                                               Neff=Neff, ret_rad=True)
    
    ######## compute dofs during RD era if dofs is chosen
    #### to compensate the radiation energy density dependence on dofs
    #### during RD era
    if dofs:
        Om_rat_dof=Omega_rad_dof(a, dir0=dir0, Neff=Neff)
        OmR0 = OmR0*Om_rat_dof
        
    ##### compute other contributions Om_Lam0 and Om_mat0, and total Om_tot
    OmM0 = 1 - OmL0 - OmR0
    Om_matt = (a/a0)**(-3)*OmM0
    Om_rad = (a/a0)**(-4)*OmR0
    Om_tot = OmL0 + Om_rad + Om_matt

    return Om_tot, Om_rad, Om_matt

def friedmann(a, a0=1, h0=0.6732, OmL0=0.6841, dofs=True, dir0='', Neff=3.):
    
    """
    Function that uses Friedmann equations to compute the eos (w) and the time derivatives
    of the scale factor for the radiation energy density described in Omega_vs_a().
    
    Arguments:
        a -- array of scale factors
        a0 -- reference value of the scale factor at present time (default is 1)
        h0 -- present-time value of the Hubble rate H0 = h0 x 100 km/s/Mpc
              (default is 67.32 km/s/Mpc based on CMB observations)
        OmL0 -- present-time content of dark energy
        dofs -- option to compensate the rad energy density using dofs during RD era
        dir0 -- directory where file for dofs during RD era is stored
    
    Reference: Y. He, A. Roper Pol, A. Brandenburg, "Modified propagation of
    gravitational waves from the early radiation era," in preparation.
    """
    
    g0, g0s, T0, H0 = values_0(h0=h0, neut=True, Neff=Neff)
    Om_tot, Om_rad, Om_matt = Omega_vs_a(a, a0=a0, h0=h0, OmL0=OmL0, dofs=dofs,
                                         dir0=dir0, Neff=Neff)
    
    # equation of state
    w = (1/3*Om_rad - OmL0)/(Om_rad + OmL0 + Om_matt)
    # time derivatives (cosmic time)
    add = -.5*Om_tot*H0**2*(a/a0)*(1 + 3*w)
    ad = (a/a0)*np.sqrt(Om_tot)*H0
    # time derivatives (conformal time)
    app = .5*a**3*Om_tot*(1 - 3*w)*H0**2
    ap = ad*a
    
    return w, ad, add, ap, app

def friedmann_solver(a, a0=1., h0=0.6732, OmL0=0.6841, dir0='', dofs=True, Neff=3.,
                     return_all=False, save=True, nm_fl=''):
    
    """
    Function that uses Friedmann equations and solve them numerically to obtain
    the evolution of a(\eta) and a(t) using the Omega distribution obtained from
    the function Omega_vs_a.
    
    Reference: Y. He, A. Roper Pol, A. Brandenburg, "Modified propagation of
    gravitational waves from the early radiation era," in preparation.
    """
    
    Om_tot, Om_rad, Om_matt = Omega_vs_a(a, a0=a0, h0=h0, OmL0=OmL0, dir0=dir0,
                                         dofs=dofs, Neff=Neff)
    
    #### numerically compute the arrays of cosmic and conformal times
    #### from the energy density ratio
    difft = np.zeros(len(a))
    diffeta = np.zeros(len(a))
    if a[0] > 1e-18:
        print('minimum a given is ', a[0])
        print('note that the solver assumes that a[0] corresponds to',
              'times t = \eta = 0, so make sure to use an array of a with',
              'small enough values (a[0] ~ 1e-20)')
    difft[0] = 0
    diffeta[0] = 0
    
    print('Entering Friedmann solver')

    for i in range(1, len(a)):
        aas = np.logspace(np.log10(a[0]), np.log10(a[i]), 10000)
        Om_as = np.interp(aas, a, Om_tot)
        fft = 1/aas/np.sqrt(Om_as)
        ffeta = fft/aas
        difft[i] = np.trapz(fft, aas)
        diffeta[i] = np.trapz(ffeta, aas)
        
    g0, g0s, T0, H0 = values_0(h0=h0)
    t = difft/H0
    eta = diffeta/H0
    
    # Using Friedmann equations, we compute the eos (w), and the time derivatives of the
    # scale factor
    w, ad, add, ap, app = friedmann(a, h0=h0, OmL0=OmL0, dofs=dofs, dir0=dir0, Neff=Neff)
    
    if save:
        df = pd.DataFrame({'a' : a, 't' : t, 'eta': eta,
                           'ap/a': ap/a, 'app/a' : app/a})
        df.to_csv('friedmann/solution' + nm_fl + '.csv')
        with open('friedmann/README' + nm_fl + '.txt', 'w') as f:
            f.write("The file solution" + nm_fl + ".csv contains the")
            f.write("solutions from the Friedmann solver using the parameters:\n")
            f.write("a0 = %.2e, h0 = %.4f, OmL0 = %.4f, Neff = %.4f \n"%(a0, h0, OmL0, Neff))
            f.write("The solution contains 'a', 't', 'eta', 'ap/a' and 'app/a'")
        print('The results are saved in the file friedmann/solution' + nm_fl + '.csv')
        print('The input parameters used are stored in friedmann/README' + nm_fl + '.txt')
            
    print('Leaving Friedmann solver')
    
    if return_all:
        return t, eta, Om_tot, Om_rad, Om_matt, w, ad, add, ap, app
    else: return t, eta
