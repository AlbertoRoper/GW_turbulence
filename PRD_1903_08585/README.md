# GW_turbulence/PRD_1903_08585

The run directories correspond to the runs in Table 1 of the paper A. Roper Pol, S. Mandal, A. Brandenburg,
T. Kahniashvili, & A. Kosowsky, *Numerical Simulations of Gravitational Waves from Early-Universe
Turbulence,* Phys. Rev. D **102**, 083512 (2020), [arXiv:1903.08585](https://arxiv.org/abs/1903.08585),
[DOI:10.1103/PhysRevD.102.083512](https://doi.org/10.1103/PhysRevD.102.083512).

The reference for the data sets is: *Datasets for "Numerical Simulations of Gravitational Waves from
Early-Universe Turbulence,"* v2020.02.28, [Zenodo](https://zenodo.org/record/3692072) (2020).

## Abstract

We perform direct numerical simulations of magnetohydrodynamic turbulence in the early universe and numerically compute the resulting stochastic background of gravitational waves and relic
magnetic fields. These simulations do not make the simplifying assumptions of earlier analytic
work. If the turbulence is assumed to have an energy-carrying scale that is about a hundredth
of the Hubble radius at the time of generation, as expected in a first-order phase transition, the
peak of gravitational wave power will be in the mHz frequency range for a signal produced at the
electroweak scale. The efficiency of gravitational wave (GW) production varies significantly with
how the turbulence is driven. Detectability of turbulence at the electroweak scale by the planned
Laser Interferometer Space Antenna (LISA) requires anywhere from 0.5% to 10% of the thermal
plasma energy density to be in plasma motions or magnetic fields, depending on the model of the
driving process. Our results predict a new universal form below the spectral peak frequency that
is shallower than previously thought. This implies larger values of the GW energy spectra in the
low-frequency range. This extends the range where turbulence is detectable with LISA to lower
frequencies, corresponding to higher energy scales than the assumed energy-carrying scale.

## Table of runs

<p align="center">
  <img src="https://github.com/AlbertoRoper/GW_turbulence/blob/master/PRD_1903_08585/tableI.png" width="700">

## Figures

<p align="center">

<!--![Figure 1](https://github.com/AlbertoRoper/GW_turbulence/blob/master/PRD_1903_08585/plots/EGW_EM_vs_k_ini2.png)-->
  
<img src="https://github.com/AlbertoRoper/GW_turbulence/blob/master/PRD_1903_08585/plots/EGW_EM_vs_k_ini2.png" width="500">
 
**Fig. 1**: Magnetic and GW energy spectra for run ini2 averaged over late times (t > 1.1), after the GW spectrum has
started to fluctuate around a steady state.

<!--![Figure 3](https://github.com/AlbertoRoper/GW_turbulence/blob/master/PRD_1903_08585/plots/EGW_vs_kt.png)-->

<p align="center">
<img src="https://github.com/AlbertoRoper/GW_turbulence/blob/master/PRD_1903_08585/plots/EGW_vs_kt.png" width="500">

**Fig. 3**: GW spectral energy versus time for four values of  <img src="https://render.githubusercontent.com/render/math?math=k">,
demonstrating the <img src="https://render.githubusercontent.com/render/math?math=k^2"> scaling at early times for run ini2.

<!--![Figure 4](https://github.com/AlbertoRoper/GW_turbulence/blob/master/PRD_1903_08585/plots/OmGW_vs_f_ini.png)-->
<!--![Figure 4](https://github.com/AlbertoRoper/GW_turbulence/blob/master/PRD_1903_08585/plots/hc_vs_f_ini.png)-->
  
<p align="center">
<img src="https://github.com/AlbertoRoper/GW_turbulence/blob/master/PRD_1903_08585/plots/OmGW_vs_f_ini.png" width="500">
<img src="https://github.com/AlbertoRoper/GW_turbulence/blob/master/PRD_1903_08585/plots/hc_vs_f_ini.png" width="500">

**Fig. 4**: Spectra of <img src="https://render.githubusercontent.com/render/math?math=h^2">, <img src="https://render.githubusercontent.com/render/math?math=h_0^2 \Omega_{\rm GW} (f)"> and <img src="https://render.githubusercontent.com/render/math?math=h_{\rm c} (f)">,
evaluated at the present time, along with the LISA power law sensitivity curve
(green dot-dashed line) to a stochastic GW background assuming four years of mission, and a threshold signal-to-noise
ratio of 10.

<p align="center">
<!--![Figure 5](https://github.com/AlbertoRoper/GW_turbulence/blob/master/PRD_1903_08585/plots/OmM_vs_t.png)-->
<!--![Figure 5](https://github.com/AlbertoRoper/GW_turbulence/blob/master/PRD_1903_08585/plots/OmGW_vs_t.png)-->
  
<img src="https://github.com/AlbertoRoper/GW_turbulence/blob/master/PRD_1903_08585/plots/OmM_vs_t.png" width="500">
<img src="https://github.com/AlbertoRoper/GW_turbulence/blob/master/PRD_1903_08585/plots/OmGW_vs_t.png" width="500">

**Fig. 5**: Evolution of <img src="https://render.githubusercontent.com/render/math?math=\Omega_{\rm M, K}"> (top)
and <img src="https://render.githubusercontent.com/render/math?math=\Omega_{\rm GW}"> (bottom) for runs
with initial energy (ini1–3) and runs where energy is driven through monochromatic forcing (hel1–2 and ac1).
Note that the energy densities are normalized with the radiation energy
density at the time of generation.

<!--![Figure 6](https://github.com/AlbertoRoper/GW_turbulence/blob/master/PRD_1903_08585/plots/OmGW_vs_f_driven.png)-->
  
<p align="center">
<img src="https://github.com/AlbertoRoper/GW_turbulence/blob/master/PRD_1903_08585/plots/OmGW_vs_f_driven.png" width="500">
<img src="https://github.com/AlbertoRoper/GW_turbulence/blob/master/PRD_1903_08585/plots/hc_vs_f_driven.png" width="500">

**Fig. 6**: Similar to Fig. 4, but for runs hel1–3, noh1, and ac1.
  
<p align="center">
<img src="https://github.com/AlbertoRoper/GW_turbulence/blob/master/PRD_1903_08585/plots/OmGW_vs_OmMK.png" width="500">
  
**Fig. 7**: <img src="https://render.githubusercontent.com/render/math?math=\Omega^{\rm sat}_{\rm GW}">
versus <img src="https://render.githubusercontent.com/render/math?math=\Omega^{\rm sat}_{\rm M, K}">.
The quadratic dependence inferred from the +2 slope of the lines holds within runs of the
same type. Note that runs ini3 (N = 10) and hel4 (N = 1000)
in green have different stirring scales than the rest of the runs
(N = 100).

