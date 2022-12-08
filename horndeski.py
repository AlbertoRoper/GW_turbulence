"""
horndeski.py is a Python routine that contains functions relevant
for Horndeski theories of modified gravity.
Author: Alberto Roper Pol
Date: 27/11/2022

The reference is Y. He, A. Roper Pol, A. Brandenburg, "Modified propagation
of gravitational waves from the early radiation era," submitted to JCAP (2022).
"""

import astropy.constants as const
import astropy.units as u
import pandas as pd
import numpy as np
import spectra as sp

def parameterizations_alpM(eta, alpM0=1., a=[], Omega=[], Omega_mat=[], OmM0=0, choice='0', n=1):
    
    """
    Function that returns the different parameterizations of the Horndeski parameter
    alpha_M.
    
    Choices are:
        - 0: constant \alpha_M
        - I: proportional to a^n
        - II: proportional to the dark-energy density
        - III: proportional to radiation + dark-energy densities
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
    
    Returns:
        - DDs: damping (boosting) function
    Arguments:
        - eta: array of times over cosmological propagation of GWs
        - alpM: time evolution of the Horndeski parameter alpha_M
        - HH: conformal Hubble rate
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
    
    Returns:
        - DeltaT: speed delay
    Arguments:
        - eta: array of times over cosmological propagation of GWs
        - alpT: time evolution of the Horndeski parameter alpha_T
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
    
    Reference is Y. He, A. Roper Pol, A. Brandenburg, "Modified propagation
    of gravitational waves from the early radiation era,"
    submitted to JCAP (2022); see equation 3.8.
    """
    
    if len(h0) == 1 and h0[0] == 0: h0 = k**0
    if len(g0) == 1 and g0[0] == 0: g0 = k**2
    if len(alpM) == 1 and alpM[0] == 0: alpM = eta**0*0
    if len(alpT) == 1 and alpT[0] == 0: alpT = eta**0*0
        
    eta_ij, kij = np.meshgrid(eta, k, indexing='ij')
    eta_ij, h0_ij = np.meshgrid(eta, h0, indexing='ij')
    eta_ij, g0_ij = np.meshgrid(eta, g0, indexing='ij')
    DD_ij, kij = np.meshgrid(DD, k, indexing='ij')
    HH_ij, kij = np.meshgrid(HH, k, indexing='ij')
    alpM_ij, kij = np.meshgrid(alpT, k, indexing='ij')
    alpT_ij, kij = np.meshgrid(alpM, k, indexing='ij')
    
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
    after averaging over oscillations in k eta.
    
    Reference is Y. He, A. Roper Pol, A. Brandenburg, "Modified propagation
    of gravitational waves from the early radiation era," submitted to
    JCAP (2022); see equation 3.19.
    """

    if len(alpM) == 1 and alpM[0] == 0: alpM = eta**0*0
    if len(alpT) == 1 and alpT[0] == 0: alpT = eta**0*0
        
    eta_ij, kij = np.meshgrid(eta, k, indexing='ij')
    DD_ij, kij = np.meshgrid(DD, k, indexing='ij')
    HH_ij, kij = np.meshgrid(HH, k, indexing='ij')
    alpM_ij, kij = np.meshgrid(alpM, k, indexing='ij')
    alpT_ij, kij = np.meshgrid(alpT, k, indexing='ij')
    
    cT_ij2 = 1 + alpT_ij
    
    TT = 1 + .5*alpT_ij
    TT = TT + alpM0**2*alpM_ij**2*HH_ij**2/32/kij**4/cT_ij2
    TT = TT + alpM0*alpM_ij**2*HH_ij**2/8/kij**3/cT_ij2
    TT = TT + (alpM_ij**2*HH_ij**2*(1 + 1/cT_ij2) + alpM0**2)/8/kij**2
    TT = TT + alpM0/2/kij
    TT = TT*np.exp(-2*DD_ij)
    
    return TT

    
def WKB_envelope_late_times(k, eta, HH, DD, alpT=[0], alpM0=0):
    
    """
    Function that uses the WKB approximation to obtain the solution to the modified
    GW equation in the propagation regime (absence of sources), it takes the envelope
    after averaging over oscillations in k \eta at late times.
    
    Reference is Y. He, A. Roper Pol, A. Brandenburg, "Modified propagation
    of gravitational waves from the early radiation era," submitted to
    JCAP (2022); see equation 3.20.
    """

    if len(alpT) == 1 and alpT[0] == 0: alpT = eta**0*0
        
    eta_ij, kij = np.meshgrid(eta, k, indexing='ij')
    DD_ij, kij = np.meshgrid(DD, k, indexing='ij')
    HH_ij, kij = np.meshgrid(HH, k, indexing='ij')
    alpT_ij, kij = np.meshgrid(alpT, k, indexing='ij')
    
    cT_ij2 = 1 + alpT_ij
    
    TT = 1 + alpT_ij + (1 + alpM0/2/kij)**2
    TT = .5*TT*np.exp(-2*DD_ij)
    
    return TT

def WKB_envelope_late_times_const(k, alpT=0, alpM0=0):
    
    """
    Function that uses the WKB approximation to obtain the solution to the modified
    GW equation in the propagation regime (absence of sources), it takes the envelope
    after averaging over oscillations in k eta at late times.
    
    Reference is Y. He, A. Roper Pol, A. Brandenburg, "Modified propagation
    of gravitational waves from the early radiation era," submitted to
    JCAP (2022); see equation 3.20.
    """
    
    TT = .5*(1 + alpT + (1 + alpM0/2/k)**2)
    
    return TT