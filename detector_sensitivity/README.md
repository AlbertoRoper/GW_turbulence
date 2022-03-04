This directory contains numerical values of the sensitivities of
different detectors to gravitational wave (GW) signals.

- The sensitivities (given in ratio of GW energy density to
present time critical energy density, Omega_GW) and power
law sensitivities (PLS) for LISA and Taiji are computed
using the routine
[interferometry.py](https://github.com/AlbertoRoper/GW_turbulence/blob/master/interferometry.py),
following the methodology of the appendix B in:
A. Roper Pol, S. Mandal, A. Brandenburg, T. Kahniashvili, "Polarization
of gravitational waves from helical MHD turbulent sources," JCAP, in press
(2022), [arXiv:2107.05356](https://arxiv.org/pdf/2107.05356.pdf).

  - 'LISA_Omega.csv' 
  - 'LISA_Omega_PLS.csv'
  - 'LISA_Xi.csv'
  - 'LISA_Xi_PLS.csv'
  - 'Taiji_Omega.csv'
  - 'Taiji_Omega_PLS.csv'
  - 'Taiji_Xi.csv'
  - 'Taiji_Xi_PLS.csv'
  - 'LISA_Taiji_Xi.csv'
  - 'LISA_Taiji_Xi_PLS.csv'

- The monopole and dipole response functions of LISA and Taiji are computed using the routine
[interferometry.py](https://github.com/AlbertoRoper/GW_turbulence/blob/master/interferometry.py),
following the methodology of the appendix B in:
A. Roper Pol, S. Mandal, A. Brandenburg, T. Kahniashvili, "Polarization
of gravitational waves from helical MHD turbulent sources," JCAP, in press
(2022), [arXiv:2107.05356](https://arxiv.org/pdf/2107.05356.pdf).

  - 'LISA_response_f.csv'
  - 'Taiji_response_f.csv'

- The monopole response functions (I and V Stokes parameters) obtained by cross-correlating channels of the LISA-Taiji
network are stored in the directory
['LISA_Taiji'](https://github.com/AlbertoRoper/GW_turbulence/tree/master/detector_sensitivity/LISA_Taiji),
from G. Orlando, M. Pieroni, A. Ricciardone, "Measuring Parity Violation
in the Stochastic Gravitational Wave Background with the LISA-Taiji network,"
JCAP 03, 069 (2021), [arXiv:2011.07059](https://arxiv.org/abs/2011.07059);
see Figure 2.

  - 'MAC_V.csv'
  - 'MAD_V.csv'
  - 'MEC_V.csv'
  - 'MED_I.csv'
  - 'MED_V.csv'

- The power law integrated sensitivities (PLS) of BBO and DECIGO are stored
in the directory
['power-law-integrated_sensitivities'](https://github.com/AlbertoRoper/GW_turbulence/tree/master/detector_sensitivity/power-law-integrated_sensitivities),
from K. Schmitz, "New Sensitivity Curves for Gravitational-Wave Signals from Cosmological Phase
Transitions," JHEP 01, 097 (2021), [arXiv:2002.04615](https://arxiv.org/abs/2002.04615);
available in [Zenodo](https://doi.org/10.5281/zenodo.3689582).

  - 'plis_BBO.dat'
  - 'plis_DECIGO.dat'
