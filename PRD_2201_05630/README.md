# GW_turbulence/PRD_2201_05630

The run directories correspond to the runs in Table I of the paper A. Roper Pol,
C. Caprini, A. Neronov, D. Semikoz, *The gravitational wave signal from primordial magnetic fields in the Pulsar
Timing Array frequency band,* [arXiv:2201.05630](https://arxiv.org/abs/2201.05630) (2022).

The reference for the data sets is: *Datasets for "The gravitational wave signal from primordial magnetic fields in the Pulsar
Timing Array frequency band,"* v2021.12.15, [Zenodo](https://zenodo.org/record/5782752) (2022).

## Abstract

The NANOGrav, Parkes, European, and International pulsar timing array (PTA) collaborations
have reported evidence for a common-spectrum process that can potentially correspond to a stochastic
gravitational wave background (SGWB) in the 1–100 nHz frequency range.
We consider the scenario in which this signal is produced by magnetohydrodynamic (MHD) turbulence in the early
universe, induced by a non-helical primordial magnetic field at the energy scale corresponding to the
quark confinement phase transition. We perform MHD simulations to study the dynamical evolution of the
magnetic field and compute the resulting SGWB. We show that the SGWB output from
the simulations can be very well approximated by assuming that the magnetic anisotropic stress is
constant in time, over a time interval related to the eddy turnover time. The analytical spectrum
that we derive under this assumption features a change of slope at a frequency corresponding to
the GW source duration that we confirm with the numerical simulations. We compare the SGWB
signal with the PTA data to constrain the temperature scale at which the SGWB is sourced, as well
as the amplitude and characteristic scale of the initial magnetic field. We find that the generation
temperature is constrained to be in the 1–200 MeV range, the magnetic field amplitude must be
\> 1% of the radiation energy density at that time, and the magnetic field characteristic scale is
constrained to be \> 10% of the horizon scale. We show that the turbulent decay of this magnetic
field will lead to a field at recombination that can help to alleviate the Hubble tension and can be
tested by measurements in the voids of the Large Scale Structure with gamma-ray telescopes like
the Cherenkov Telescope Array.

## Tables of runs

<p align="center">
<img src="https://github.com/AlbertoRoper/GW_turbulence/blob/master/PRD_2201_05630/tableI.png" width="800">

## Figures

<p align="center">  
<img src="https://github.com/AlbertoRoper/GW_turbulence/blob/master/PRD_2201_05630/plots/png/Omega_GW_runsA.png" width="500">
<img src="https://github.com/AlbertoRoper/GW_turbulence/blob/master/PRD_2201_05630/plots/png/Omega_GW_runsB.png" width="500">
<img src="https://github.com/AlbertoRoper/GW_turbulence/blob/master/PRD_2201_05630/plots/png/Omega_GW_runsC.png" width="500">
<img src="https://github.com/AlbertoRoper/GW_turbulence/blob/master/PRD_2201_05630/plots/png/Omega_GW_runsD.png" width="500">
<img src="https://github.com/AlbertoRoper/GW_turbulence/blob/master/PRD_2201_05630/plots/png/Omega_GW_runsE.png" width="500">
 
**Fig. 1**: Simulated GW spectra <img src="https://render.githubusercontent.com/render/math?math=h^2 \Omega_{\rm {GW}}^0 (k)">
of runs A–E (dots) compared to the analytical model developed in subsection (II D) assuming constant magnetic stresses:
Eq. (18) at <img src="https://render.githubusercontent.com/render/math?math=t_{\rm {fin}}"> (thin gray lines) and its envelope Eq. (24)
(dash-dotted black lines).
The maximal values over one oscillation of the numerical outputs at each wave number
are shown in different colors for runs with different domain sizes, and are combined to show the GW spectra from
sub-horizon scales up to the scales where the inertial range is developed.
Runs A1 to E1 (blue dots) are computing the smallest scales resolved, up to a Nyquist wave number of
<img src="https://render.githubusercontent.com/render/math?math=k_{\rm {Ny}} {\cal H}_*^{-1}"> = 126.
The spectra are shown in terms of <img src="https://render.githubusercontent.com/render/math?math=k {\cal H}_*^{-1}">
and compensated by <img src="https://render.githubusercontent.com/render/math?math=(g_*/10)^{-{1\over3}}">,
such that they can be scaled to the specific value of the comoving Hubble rate
at the time of generation <img src="https://render.githubusercontent.com/render/math?math={\cal H}_*">
and to different values of <img src="https://render.githubusercontent.com/render/math?math=g_*">.

<p align="center">
<img src="https://github.com/AlbertoRoper/GW_turbulence/blob/master/PRD_2201_05630/plots/png/dte_dtf_fit.png" width="500">
<img src="https://github.com/AlbertoRoper/GW_turbulence/blob/master/PRD_2201_05630/plots/png/dte_ratio_fit.png" width="500">
<img src="https://github.com/AlbertoRoper/GW_turbulence/blob/master/PRD_2201_05630/plots/png/Omega_GW_runs_model.png" width="500">

**Fig. 2**: *Upper Panel (left)*: the points represent the eddy
turnover times of the simulations <img src="https://render.githubusercontent.com/render/math?math=\delta t_{\rm e}">
and the corresponding values of <img src="https://render.githubusercontent.com/render/math?math=\delta t_{\rm {fin}}">
obtained by fitting the break into <img src="https://render.githubusercontent.com/render/math?math=k^3">
(occurring at <img src="https://render.githubusercontent.com/render/math?math=1/\delta t_{\rm {fin}}"> according to the analytical model).
The dashed line represents the fit of Eq. (36).
*Upper Panel (right)*: the points represent the ratio <img src="https://render.githubusercontent.com/render/math?math={\cal G}">
between the numerical and the analytical SGWB amplitudes at the peak.
The dashed line represents the fit of Eq. (37).
*Bottom panel*: the SGWBs computed with the analytical model given in Eq. (24) (dot-dashed lines), and from
the adjusted model given in Eq. (38) (solid lines), are compared to the results of the MHD simulations (colored points).
The compensated model uses the empirical fits shown in the upper and middle panel.
  
<p align="center">
<img src="https://github.com/AlbertoRoper/GW_turbulence/blob/master/PRD_2201_05630/plots/png/PTA_Om_vs_beta.png" width="500">
<img src="https://github.com/AlbertoRoper/GW_turbulence/blob/master/PRD_2201_05630/plots/png/Om_PTA_f.png" width="500">
 
**Fig. 3**: *Left*: 1<img src="https://render.githubusercontent.com/render/math?math=\sigma"> and
2<img src="https://render.githubusercontent.com/render/math?math=\sigma"> contours of the amplitude
<img src="https://render.githubusercontent.com/render/math?math=h^2 \Omega_{\rm {yr}}"> vs. slope
<img src="https://render.githubusercontent.com/render/math?math=\beta"> (cf. Eqs. (43) and (45)) derived
from the NANOGrav data set for both the broken PL (blue) and single PL (green) fits, and from the PPTA (red),
EPTA (purple), and IPTA (black) data sets for the single PL fit.
The grey shaded area shows the slopes <img src="https://render.githubusercontent.com/render/math?math=\beta \in (1, 3)">
characteristic of the SGWB produced by primordial MHD turbulence below the spectral peak, cf. subsection (II D 1).
*Right*: Shaded regions: range of the SGWB spectra <img src="https://render.githubusercontent.com/render/math?math=h^2 \Omega_{\rm {GW}}^0 (f)">
of Eqs. (43) and (45), corresponding to the 2<img src="https://render.githubusercontent.com/render/math?math=\sigma">
contours given in the upper panel.
The vertical line shows the reference frequency <img src="https://render.githubusercontent.com/render/math?math=f_{\rm {yr}}">.
Dashed lines: 2<img src="https://render.githubusercontent.com/render/math?math=\sigma"> maximum amplitude at each frequency —
such that larger amplitudes are, in principle, excluded by the PTA observations.
  
<p align="center">
<img src="https://github.com/AlbertoRoper/GW_turbulence/blob/master/PRD_2201_05630/plots/png/OmGW_compatible_T100.png" width="500">
<img src="https://github.com/AlbertoRoper/GW_turbulence/blob/master/PRD_2201_05630/plots/png/OmGW_compatible_T150.png" width="500">

**Fig. 4**:  For <img src="https://render.githubusercontent.com/render/math?math=T_*"> = 150 MeV and
<img src="https://render.githubusercontent.com/render/math?math=g_*"> = 15 in the upper
panel, and <img src="https://render.githubusercontent.com/render/math?math=T_*"> = 100 MeV and
<img src="https://render.githubusercontent.com/render/math?math=g_*"> = 10 in the lower
panel, we show the upper boundary (solid lines) and the lower boundary (dashed lines) of the regions compatible
with the PTA data at 2<img src="https://render.githubusercontent.com/render/math?math=\sigma">.
To be compatible with NANOGrav with broken PL, each SGWB spectrum must lie in the region within the blue solid and dashed
lines; with NANOGrav with single PL, within the green solid and dashed lines; with PPTA, within the red solid
and dashed lines; with EPTA, within the purple solid and dashed lines; with IPTA, within the black solid and dashed lines.
The shaded areas correspond to the range of allowed values <img src="https://render.githubusercontent.com/render/math?math=h^2 \Omega_{\rm {GW}}^0 (f)">
of Eqs. (43) and (45), restricted to the range of slopes of interest for a MHD-produced SGWB, i.e.
<img src="https://render.githubusercontent.com/render/math?math=\beta \in (1, 3)">.
The magnetic field characteristic scale is bound to <img src="https://render.githubusercontent.com/render/math?math=k_* \geq 2\pi {\cal H}_*">
and the magnetic energy densities to <img src="https://render.githubusercontent.com/render/math?math=\Omega^*_{\rm M} \leq"> 0.1.
The vertical lines show the upper bound of the PTA frequency subset to which we restrain the analysis:
<img src="https://render.githubusercontent.com/render/math?math=f \simeq 1.25 \times 10^{-8}"> Hz for the
single PL cases (dot-dashed line), and <img src="https://render.githubusercontent.com/render/math?math=f \simeq 9 \times 10^{-9}"> Hz
for the NANOGrav broken PL case (dashed line).
  
<p align="center">
<img src="https://github.com/AlbertoRoper/GW_turbulence/blob/master/PRD_2201_05630/plots/png/pars_compatible_T200.png" width="330">
<img src="https://github.com/AlbertoRoper/GW_turbulence/blob/master/PRD_2201_05630/plots/png/pars_compatible_T150.png" width="330">
<img src="https://github.com/AlbertoRoper/GW_turbulence/blob/master/PRD_2201_05630/plots/png/pars_compatible_T100.png" width="330">
<img src="https://github.com/AlbertoRoper/GW_turbulence/blob/master/PRD_2201_05630/plots/png/pars_compatible_T50.png" width="330">
<img src="https://github.com/AlbertoRoper/GW_turbulence/blob/master/PRD_2201_05630/plots/png/pars_compatible_T20.png" width="330">
<img src="https://github.com/AlbertoRoper/GW_turbulence/blob/master/PRD_2201_05630/plots/png/pars_compatible_T10.png" width="330">
<img src="https://github.com/AlbertoRoper/GW_turbulence/blob/master/PRD_2201_05630/plots/png/pars_compatible_T5.png" width="330">
<img src="https://github.com/AlbertoRoper/GW_turbulence/blob/master/PRD_2201_05630/plots/png/pars_compatible_T2.png" width="330">
<img src="https://github.com/AlbertoRoper/GW_turbulence/blob/master/PRD_2201_05630/plots/png/pars_compatible_T1.png" width="330">


**Fig. 5**:  For different values of <img src="https://render.githubusercontent.com/render/math?math=T_* \in (1, 200)"> MeV,
we show the allowed regions in the <img src="https://render.githubusercontent.com/render/math?math=(k_*, \Omega_{\rm M}^*)">
parameter space, derived as described in the main text from the 2<img src="https://render.githubusercontent.com/render/math?math=\sigma">
results of NANOGrav using the broken PL (blue) and single PL (green) fits, and from the 2<img src="https://render.githubusercontent.com/render/math?math=\sigma">
results of EPTA (purple), PPTA (red), and IPTA (black) using the single PL fits.
The vertical and horizontal dot-dashed lines show the physical limits
<img src="https://render.githubusercontent.com/render/math?math=k_* \geq 2 \pi {\cal H}_*"> and
<img src="https://render.githubusercontent.com/render/math?math=\Omega_{\rm M}^* \leq 0.1"> respectively: the allowed
parameter region lies within the rectangle.
The wave number of the largest processed eddies <img src="https://render.githubusercontent.com/render/math?math=k_*|_{\rm {LPE}}"> 
is also shown (dash-dotted diagonal line).
