"""
GW_analytical.py is a Python routine that contains the functions to make
calculations related to the analytical description of GW backgrounds produced
by MHD turbulence in the early universe.
It also includes some general mathematical functions that are useful.

Author: Alberto Roper Pol
created: 01/12/2021
"""

import numpy as np
import matplotlib.pyplot as plt
import plot_sets
import os

def smoothed_dbpl(k, A=1, kb=1, ks=10, a=2, b=0, c=11/3, alpha1=2, alpha2=2):
    
    """
    Function that computes the smoothed double broken power law used in
    Y. He, A. Roper Pol, A. Brandenburg, "Modified propagation of gravitational
    waves from the early radiation era," submitted to JCAP (2022) for analytical
    calculations; see eq. 5.1.
    """
        
    S = A*(k/ks)**a/(1 + (k/kb)**((a - b)*alpha1))**(1/alpha1)
    S = S/(1 + (k/ks)**((b + c)*alpha2))**(1/alpha2)
    S *= (1 + (ks/kb)**(a - b))
    
    return S

def integrated_peak_Sk_smoothed_dbpl(plot=True, ret=True, save=True):
    
    """
    Function that computes the relation between the integral of the spectrum and
    its spectral amplitude and the relation between the position of the maximum of 
    k x S (k) and the parameter k_*, where S(k) is the smoothed double broken power
    law defined in the function smoothed_dbpl.
    
    It plots the dependence of these relations as a function of k_*.
    """
    
    kss = np.logspace(0, np.log10(6000), 1000)
    ints = np.zeros(len(kss))
    krat = np.zeros(len(kss))
    k = np.logspace(-3, 4, 20000)
    for i in range(0, len(kss)):
        S = smoothed_dbpl(k, A=1, kb=1, ks=kss[i], a=2, b=0, c=11/3, alpha1=2, alpha2=2)
        ints[i] = np.trapz(S, k)
        krat[i] = kss[i]/k[np.argmax(k*S)]
    
    if plot:
        plt.plot(kss, ints/kss, color='blue')
        plt.xscale('log')
        plt.xlabel(r'$k_*$')
        plt.ylabel(r'$(\int S (k) {\rm \, d} k)/k_*$')
        plt.ylim(1.05, 1.35)
        plt.xlim(1, 1e3)
        plt.xticks(np.logspace(0, 3, 4))
        plot_sets.axes_lines()
        if save:
            ffl = 'plots/intS_vs_ks.pdf'
            print('Saving figure %s'%ffl)
            plt.savefig(ffl, bbox_inches='tight')
        plt.figure()
        plt.plot(kss, krat, color='blue')
        plt.xscale('log')
        plt.xlabel(r'$k_*$')
        plt.ylabel(r'$k_{\rm max}/k_*$')
        plt.xlim(1, 1e3)
        plt.xticks(np.logspace(0, 3, 4))
        plt.ylim(.98, 1.2)
        plt.yticks([1, 1.05, 1.1, 1.15, 1.2])
        plt.hlines(1.143, 1e0, 1e3, color='black', lw=.8)
        plot_sets.axes_lines()
        if save:
            ffl = 'plots/kpeak_vs_ks.pdf'
            print('Saving figure %s'%ffl)
            plt.savefig(ffl, bbox_inches='tight')
 
    if ret:
        return kss, ints, krat