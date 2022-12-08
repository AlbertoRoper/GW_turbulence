"""
pta.py is a Python routine with functions used in the analysis of observations
by pulsar timing array (PTA) collaborations: NANOGrav, PPTA, EPTA, and IPTA;
in the context of GW backgrounds produced by MHD turbulence in the early
universe.

The main reference for the study is:
- A. Roper Pol, C. Caprini, A. Neronov, D. Semikoz,
"The gravitational wave signal from primordial magnetic fields in the
Pulsar Timing Array frequency band," Phys. Rev. D 105, 123502 (2022),
arXiv:2201.05630.

Author: Alberto Roper Pol
Date: 01/12/2021
"""

import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt
