# GW_turbulence/horndeski

The run directories correspond to the runs in Table I of the paper Y.He, A. Roper Pol, A. Brandenburg
*Modified propagation of gravitational waves from the early radiation era,* submitted to JCAP (2022).

The reference for the data sets is: *Datasets of “Modified propagation of gravitational
waves from the early radiation era", doi: 10.5281/zenodo.7408601 (v2022.12.07),
[Zenodo](https://zenodo.org/record/7408601).

## Abstract

We study the propagation of cosmological gravitational wave (GW) backgrounds
produced in the early radiation era until the present day in modified theories of gravity.
Compared to general relativity (GR), we study the effects that Horndeski parameters, such
as the run rate of the effective Planck mass αM and the tensor speed excess $\alpha_{\rm T}$, have on the
present-day GW spectrum using the WKB estimate, which provides an analytical description
but fails at superhorizon scales, and numerical simulations that allows us to go beyond the
WKB approximation. We show that αT makes relatively insignificant changes to the GR
solution, especially taking into account the constraints on its value from GW observations by
the LIGO-Virgo collaboration, while αM can introduce modifications to the spectral slopes of
the GW energy spectrum in the low-frequency regime depending on the considered time evo-
lution parameterization of αM. This effect is additional to the damping or growth occurring
equally at all scales that can be predicted by the WKB approximation. We discuss the obser-
vational implications in light of the recent observations by pulsar timing array collaborations
and future detectors such as SKA, LISA, and DECIGO.

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
 
<p align="center">
<img src="https://github.com/AlbertoRoper/GW_turbulence/blob/master/PRD_2201_05630/plots/png/fbr_T_interp.png" width="500">

**Fig. 6**:  Range of frequencies at which the <img src="https://render.githubusercontent.com/render/math?math=f^3">
to <img src="https://render.githubusercontent.com/render/math?math=f^1"> break,
<img src="https://render.githubusercontent.com/render/math?math=f_{\rm {br}}">, occurs, for the parameters
<img src="https://render.githubusercontent.com/render/math?math=(k_*, \Omega_{\rm M}^*)"> compatible with
the results of each of the PTA collaboration, for different <img src="https://render.githubusercontent.com/render/math?math=T_*">
(cf. Fig. 5), in the limit <img src="https://render.githubusercontent.com/render/math?math=\Omega_{\rm M}^* \leq 0.1">.
 
<p align="center">
<img src="https://github.com/AlbertoRoper/GW_turbulence/blob/master/PRD_2201_05630/plots/png/B_vs_l.png" width="700">

**Fig. 7**: Region in the magnetic field parameter space, given by its comoving amplitude B and characteristic scale
l, compatible with the observations of the different PTA collaborations: in blue, NANOGrav with broken PL; in
green, NANOGrav with single PL; in red, PPTA; in purple, EPTA; in black, IPTA.
The parameter space region accessible to LISA is shown in light blue.
The horizontal dot-dashed lines show the bounds <img src="https://render.githubusercontent.com/render/math?math=\Omega_{\rm M}^* \leq 1">,
and <img src="https://render.githubusercontent.com/render/math?math=\Omega_{\rm M}^* \leq 0.1"> from Nucleosynthesis.
The black dash-dotted diagonal lines show: i) the magnetic field amplitude when the
characteristic scale corresponds to the largest processed eddies at the QCD phase transition
 <img src="https://render.githubusercontent.com/render/math?math=l_*|_{\rm {LPE}}"> (cf. Eq. (52));
ii) the magnetic field amplitude reached at recombination (cf. Eq. (53)).
The dash-dotted red and brown lines show the evolutionary paths of the extremities of the parameter space region compatible with the PTA results up to
recombination, following compressible (red) and incompressible (brown) MHD free decay.
The solid red and brown lines indicate the evolutionary paths of an initial field with
<img src="https://render.githubusercontent.com/render/math?math=k_* = 2 \pi {\cal H}_*"> and
<img src="https://render.githubusercontent.com/render/math?math=\Omega_{\rm M}^* = 0.1"> at
<img src="https://render.githubusercontent.com/render/math?math=T_*"> = 100 MeV and
<img src="https://render.githubusercontent.com/render/math?math=g_*"> = 10 (right red dot), and at
<img src="https://render.githubusercontent.com/render/math?math=T_*"> = 150 MeV and
<img src="https://render.githubusercontent.com/render/math?math=g_*"> = 15 (left red dot).
The green line indicates the upper limit B <img src="https://render.githubusercontent.com/render/math?math=\lesssim"> 0.1 nG,
and the range <img src="https://render.githubusercontent.com/render/math?math=B_{\rm {rec}} \in (0.013, 0.1)"> nG,
proposed to alleviate the Hubble tension, both derived in ref. [62] from CMB
constraints on the baryon clumping.
The black solid diagonal lines show the Fermi Large Area Telescope (LAT) lower
bound on the intergalactic magnetic field from timing of the blazar signal (darker gray area), and from the search of
extended emission (lighter gray area) [2].
The blue line shows the expected sensitivity of CTA [63]. At larger scales,
the upper bound from Faraday rotation (FR) is shown [99], and the blue shaded region indicates the observations
of UHECR from the Perseus-Pisces supercluster [100, 101]. Note that the latter constraints refer to present time
magnetic field strength and characteristic scale, and they have been cut to avoid intersecting the evolutionary paths
from the QCD phase transition up to recombination in the plot, for clarity.
