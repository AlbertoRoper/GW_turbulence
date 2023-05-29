"""
horndeski.py is a Python routine that contains functions relevant
for GW production in the context of general theories of modified gravity.

Author: Alberto Roper Pol
Date: 27/11/2022

Reference: Y. He, A. Roper Pol, A. Brandenburg, "Modified propagation
of gravitational waves from the early radiation era," in press, JCAP (2023),
arXiv:2212.06082
"""

import astropy.constants as const
import astropy.units as u
import pandas as pd
import numpy as np
import spectra as sp
import cosmology as co
import matplotlib.pyplot as plt

# reference values and constants
Tref = 100*u.GeV    #EWPT

def parameterizations_alpM(eta, alpM0=1., a=[], Omega=[], Omega_mat=[], OmM0=0, choice='0', n=1):
    
    """
    Function that returns the different parameterizations of the Horndeski parameter
    alpha_M.
    
    Arguments:
        eta -- array of conformal times
        
    Returns:
        alpM, alpM_prime -- time evolution and time derivative of the \alpM parameter
                            according to parameterization of choice #
    
    Choices are:
        - 0: constant \alpha_M
        - I: proportional to a^n
        - II: proportional to the dark-energy density
        - III: proportional to radiation + dark-energy densities
        
    Reference: Y. He, A. Roper Pol, A. Brandenburg, "Modified propagation
    of gravitational waves from the early radiation era," in press, JCAP (2023),
    arXiv:2212.06082, eq. 4.2
    """
    
    if choice=='0': alpM = eta**0
    elif choice=='I':
        if len(a) == 0:
            alpM = eta**0
            print('for choice I you need to include a in the input of ' + \
                  'parameterizations_alpM')
        else: alpM = a**n
    elif choice=='II':
        if len(Omega) == 0:
            alpM = eta**0
            print('for choice II you need to include Omega in the input of ' + \
                 'parameterizations_alpM')
        else: alpM = 1/Omega
    elif choice=='III':
        if len(Omega_mat) == 0:
            alpM = eta**0
            print('for choice III you need to include Omega_mat in the input of ' + \
                 'parameterizations_alpM')
        elif len(Omega) == 0:
            alpM = eta**0
            print('for choice III you need to include Omega in the input of ' + \
                 'parameterizations_alpM')
        else:
            inds = np.argsort(Omega)
            if OmM0 == 0: OmM0 = np.interp(1, Omega[inds], Omega_mat[inds])
            alpM = (1 - Omega_mat/Omega)/(1 - OmM0)

    alpM = alpM*alpM0
    alpM_prime = np.gradient(alpM, eta)
    
    return alpM, alpM_prime

def damping_fact(eta, alpM, HH):
    
    """
    Computes the damping factor that appears in the modified GR solution of
    the GW equation using the WKB approximation.
    
    Arguments:
        eta -- array of times over cosmological propagation of GWs
        alpM -- time evolution of the Horndeski parameter alpha_M
        HH -- conformal Hubble rate
    
    Returns:
        DDs -- damping (boosting) function

    Reference: Y. He, A. Roper Pol, A. Brandenburg, "Modified propagation
    of gravitational waves from the early radiation era," in press, JCAP (2023),
    arXiv:2212.06082, eq. 2.5
    """
    
    DDs = np.zeros(len(eta))
    if len(alpM) == 1: alpM = alpM[0]*eta**0
    for i in range(0, len(eta)):
        DDs[i] = .5*np.trapz(alpM[:i]*HH[:i], eta[:i])

    return DDs

def DeltaT(eta, alpT):
    
    """
    Computes the GW speed delay that appears in the modified GR solution of
    the GW equation using the WKB approximation.

    Arguments:
        eta -- array of times over cosmological propagation of GWs
        alpT -- time evolution of the Horndeski parameter alpha_T
    Returns:
        DeltaT -- speed delay

    Reference: Y. He, A. Roper Pol, A. Brandenburg, "Modified propagation
    of gravitational waves from the early radiation era," in press, JCAP (2023),
    arXiv:2212.06082, eq. 3.4
    """
    
    DDT = np.zeros(len(eta))
    if len(alpT) == 1: alpT = alpT[0]*eta**0
    for i in range(0, len(eta)):
        cT = np.sqrt(1 + alpT[:i])
        DDT[i] = np.trapz(1 - cT, eta[:i])

    return DDT
    
def sol_WKB(k, eta, HH, DD, eta_ini=1, h0=[0], g0=[0], alpM=[0], alpT=[0], alpM0=0):
    
    """
    Function that uses the WKB approximation to obtain the solution to the modified
    GW equation in the propagation regime (absence of sources).
    
    It is calculated for general modified theories of gravity with \cT = 1 at all times.
    
    Arguments:
        k -- array of wave numbers
        eta -- array of conformal times
        HH -- array of Hubble rate
        DD -- array of damping factor due to MG from \alpM \neq 0
        eta_ini -- initial conformal time
        h0 -- initial condition of the strains
        g0 -- initial condition of the time derivative of the strains
        alpM -- values \alpM
        alpT -- values \alpT
        alpM0 -- initial value of \alpM

    Returns:
        SgtGR -- GW spectrum of the time derivative strains in GR
        SgtmodGR -- GW spectrum of the time derivative strains in MG
    
    Reference: Y. He, A. Roper Pol, A. Brandenburg, "Modified propagation
    of gravitational waves from the early radiation era," in press, JCAP (2023),
    arXiv:2212.06082, eq. 3.8.
    """
    
    # default initial conditions are h = 1, g = k^2
    # \alpM = 1 if not given
    if len(h0) == 1 and h0[0] == 0: h0 = k**0
    if len(g0) == 1 and g0[0] == 0: g0 = k**2
    if len(alpM) == 1 and alpM[0] == 0: alpM = eta**0*1

    eta_ij, kij = np.meshgrid(eta, k, indexing='ij')
    eta_ij, h0_ij = np.meshgrid(eta, h0, indexing='ij')
    eta_ij, g0_ij = np.meshgrid(eta, g0, indexing='ij')
    DD_ij, kij = np.meshgrid(DD, k, indexing='ij')
    HH_ij, kij = np.meshgrid(HH, k, indexing='ij')
    alpM_ij, kij = np.meshgrid(alpM, k, indexing='ij')
    
    htGR = h0_ij*np.cos(kij*(eta_ij - eta_ini))
    htGR += g0_ij/kij*np.sin(kij*(eta_ij - eta_ini))
    gtGR = -h0_ij*kij*np.sin(kij*(eta_ij - eta_ini))
    gtGR += g0_ij*np.cos(kij*(eta_ij - eta_ini))
    
    htmodGR = np.exp(-DD_ij)*(htGR + .5*alpM0/kij*h0_ij*np.sin(kij*(eta_ij - eta_ini)))
    gtmodGR = -.5*alpM_ij*HH_ij*htmodGR
    gtmodGR += np.exp(-DD_ij)*(gtGR + .5*alpM0*h0_ij*np.cos(kij*(eta_ij - eta_ini)))
    
    SgtGR = gtGR**2
    SgtmodGR = gtmodGR**2
    
    return SgtGR, SgtmodGR
    
def WKB_envelope(k, eta, HH, DD, alpM=[0], alpT=[0], alpM0=0):
    
    """
    Function that uses the WKB approximation to obtain the solution to the modified
    GW equation in the propagation regime (absence of sources), it takes the envelope
    after averaging over oscillations in k eta (assumes constant \alpT).
    
    The GW spectrum in MG is: OmGW_MG = OmGW_GR x TT
    
    Arguments:
        k -- array of wave numbers
        eta -- array of conformal times
        HH -- array of Hubble rate
        DD -- array of damping factor due to MG from \alpM \neq 0
        alpM -- values \alpM
        alpT -- values \alpT
        alpM0 -- initial value of \alpM

    Returns:
        TT -- factor TT modifying the GW spectrum in MG
    
    Reference: Y. He, A. Roper Pol, A. Brandenburg, "Modified propagation
    of gravitational waves from the early radiation era," in press, JCAP (2023),
    arXiv:2212.06082, eq. 3.20
    """

    # \alpM = 1, \alpT = 0 if not given
    if len(alpM) == 1 and alpM[0] == 0: alpM = eta**0*1
    if len(alpT) == 1 and alpT[0] == 0: alpT = eta**0*0
        
    eta_ij, kij = np.meshgrid(eta, k, indexing='ij')
    DD_ij, kij = np.meshgrid(DD, k, indexing='ij')
    HH_ij, kij = np.meshgrid(HH, k, indexing='ij')
    alpM_ij, kij = np.meshgrid(alpM, k, indexing='ij')
    #alpT_ij, kij = np.meshgrid(alpT, k, indexing='ij')
    
    cT_ij2 = 1 + alpT_ij
    
    TT = 1 + .5*alpT_ij
    TT = TT + alpM0**2*alpM_ij**2*HH_ij**2/32/kij**4/cT_ij2
    TT = TT + alpM0*alpM_ij**2*HH_ij**2/8/kij**3/cT_ij2
    TT = TT + (alpM_ij**2*HH_ij**2*(1 + 1/cT_ij2) + alpM0**2)/8/kij**2
    TT = TT + alpM0/2/kij
    TT = TT*np.exp(-2*DD_ij)
    
    return TT

def WKB_envelope_late_times(k, DD=0, alpT=0, alpM0=0):
    
    """
    Function that uses the WKB approximation to obtain the solution to the modified
    GW equation in the propagation regime (absence of sources), it takes the envelope
    after averaging over oscillations in k \eta at late times. It gets rid of
    the terms decaying in time 1/t compared to WKB_envelope
    
    The GW spectrum in MG is: OmGW_MG = OmGW_GR x TT
    
    Arguments:
        k -- array of wave numbers
        DD -- array of damping factor due to MG from \alpM \neq 0 (default is 0)
        alpT -- constant value \alpT (default is 0)
        alpM0 -- initial value of \alpM (default is 0)

    Returns:
        TT -- factor TT modifying the GW spectrum in MG
    
    Reference: Y. He, A. Roper Pol, A. Brandenburg, "Modified propagation
    of gravitational waves from the early radiation era," in press, JCAP (2023),
    arXiv:2212.06082, eq. 3.21
    """

    TT = 1 + alpT + (1 + alpM0/2/k)**2
    TT = .5*TT*np.exp(-2*DD)
    
    return TT

def damping_amplification(eta, H, a=[], Omega=[], Omega_mat=[], OmM0=0, n=1, ch='0'):
    
    """
    Function that computes the damping or the amplification of a GW signal
    from its time of generation within the RD era until present time using
    one of the four choices of alpM parameterization with time.
    
    Arguments:
        eta -- array of times
        ch -- choice number for the \alpM parameterization, from parameterizations_alpM
    Returns:
        Q -- amplification (damping) factor
        
    Reference: Y. He, A. Roper Pol, A. Brandenburg, "Modified propagation
    of gravitational waves from the early radiation era," in press, JCAP (2023),
    arXiv:2212.06082, eq. 2.5
    """
    
    alpsM0, _ = parameterizations_alpM(eta, alpM0=1, a=a, Omega=Omega, Omega_mat=Omega_mat,
                                       OmM0=OmM0, choice=ch)

    # compute damping term at present time
    D = .5*np.trapz(alpsM0*H, eta)
    Q = np.exp(-2*D)
    
    return Q

def WKB_slopes(k, alpM):
    
    """
    Function that computes the slopes of the spectral modification due to
    MG under the WKB approximation, computed in WKB_envelope_late_times,
    defined such that
    
    OmGW ~ A k^\beta
    
    is tangent to the spectrum at k
    
    Arguments:
        k -- array of wave numbers
        alpM -- value of \alpM MG parameter
        
    Returns:
        beta -- exponent of the tangent power law of the MG spectrum to k
    
    Reference: Y. He, A. Roper Pol, A. Brandenburg, "Modified propagation
    of gravitational waves from the early radiation era," in press, JCAP (2023),
    arXiv:2212.06082, eq. 5.4
    """
    
    beta = 2*alpM*(alpM + 2*k)/(8*k**2 + alpM**2 + 4*k*alpM)
    
    return beta

def compute_Qs(a, eta, ap_a, app_a, T=Tref, OmM0=co.OmM0_ref, h0=co.h0_ref):
    
    """
    Function that computes the amplification for the 4 different choices
    of alpM parameterization from the time of generation until present time
    given the temperature scale within the RD at which the signal was generated.
    
    It takes alpM0 = 1, so for other values of alpM0, one can compute them
    as Q**alpM0.
    
    Arguments:
        a -- scale factor from Friedmann solver (normalized to a0 = 1)
        eta -- conformal time (normalized to a0 = 1)
        ap_a -- a'/a (normalized to a0 = 1)
        app_a -- a''/a (normalized to a0 = 1)
        T -- temperature scale (default 100 GeV)
        OmM0 -- present day amount of matter
        h0 -- present day Hubble rate H0 = 100 h0 km/s/Mpc
        
    Returns:
        Q0, Q1, Q2, Q3 -- amplification (damping) factor of reference
                          (alpM0 = 1) for the 4 choices of alpM parameterization
    """
    
    g = co.thermal_g(T=T, s=0)
    ast = co.as_a0_rat(T=T, g=g)
    a_n, eta_n, HH_n, app_a_n, Omega, w, eta_n_0, aEQ_n, \
    aL_n, a_acc_n, eta_n_EQ, eta_n_L, eta_n_acc = \
        co.normalized_variables(a, eta, ap_a, app_a, T=T, h0=h0)
    Omega_mat = OmM0*a**(-3)
    
    eta_nn, HH_nn, a_nn, Omega_nn, Omega_mat_nn, app_nn, w_nn = \
            co.norm_variables_cut(eta_n, HH_n, a_n, Omega, Omega_mat, 
                                  eta_n_0, T=T, OmM0=OmM0, h0=h0)
    
    Q0 = damping_amplification(eta_nn, HH_nn, ch='0')
    Q1 = damping_amplification(eta_nn, HH_nn, ch='I', a=a_nn*ast, n=1)
    Q2 = damping_amplification(eta_nn, HH_nn, ch='II', Omega=Omega_nn)
    Q3 = damping_amplification(eta_nn, HH_nn, ch='III', Omega=Omega_nn,
                                  Omega_mat=Omega_mat_nn, OmM0=OmM0)
    
    return Q0, Q1, Q2, Q3

# to be moved to plotting routines
def plot_radicand_slopes():
    
    """
    Function that plots the radicand that appears in WKB slopes.
    """
    
    betas = np.linspace(-3, 3, 1000)
    plt.plot(betas, 2*betas - betas**2 + 1, color='blue')
    plt.ylim(-1, 3)
    plt.xlim(-1, 3)
    plot_sets.axes_lines()
    plt.hlines(0, -2, 3, color='black', ls='dashed', lw=.9)
    plt.hlines(2, -2, 3, color='black', ls='dashed', lw=.9)
    plt.vlines(1 - np.sqrt(2), -2, 3, color='black', ls='dashed', lw=.9)
    plt.vlines(1 + np.sqrt(2), -2, 3, color='black', ls='dashed', lw=.9)
    plt.text(-.35, 2.3, '$-1 - \sqrt{2}$')
    plt.text(1.65, 2.3, '$1 + \sqrt{2}$')
    plt.xticks([-1, 0, 1, 2, 3])