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

- Under the directory 'power-law-integrated_sensitivities' one finds
the PLS given in K. Schmitz, "New Sensitivity Curves for
Gravitational-Wave Signals from Cosmological Phase Transitions,"
JHEP 01, 097 (2021), [arXiv:2002.04615](https://arxiv.org/pdf/2002.04615.pdf);
available in [Zenodo](https://doi.org/10.5281/zenodo.3689582).

  - 'plis_BBO.dat'
  - 'plis_DECIGO.dat'