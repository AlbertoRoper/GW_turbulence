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
    S *= (1 + (ks/kb)**(a - b))*2**(1/alpha2)
    
    return S
