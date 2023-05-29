# cosmoGW

* Author: Alberto Roper Pol (alberto.roperpol@unige.ch)
* Created: 01/11/2021
* Latest stable release: [*v2.0.0 (30/05/2022)*](https://zenodo.org/record/6045844)

cosmoGW (updated version of *GW_turbulence*) is a project to produce results related to the production of cosmological gravitational wave backgrounds from different sources in the early Universe (especially from MHD turbulence).

Some of the routines use results from large-scale numerical simulations using the [Pencil Code](https://github.com/pencil-code);
see [Pencil Code Collaboration], *The Pencil Code, a modular MPI code for partial differential equations and particles: multipurpose and multiuser-maintained,* J. Open Source Softw. **6**, 2807 (2021), [arXiv:2009.08231](https://arxiv.org/abs/2009.08231), [DOI:10.21105/joss.02807](https://joss.theoj.org/papers/10.21105/joss.02807).

Other routines are part of cosmoGW, including postprocessing calculations, lower-scale numerical simulations, and plotting routines.

If you use any of the cosmoGW results, please cite this [repository](https://zenodo.org/record/6045844), the relevant reference/s listed in the routines, and the [Pencil Code paper](https://joss.theoj.org/papers/10.21105/joss.02807) and [code](https://github.com/pencil-code) when the results are based on Pencil Code simulations. I would also love to hear about your work, so don't dout to send me an email and tell me about your project!

If you have any issues, comments, suggestions, or you are just interested in discussing any of the presented results, you are more than welcome to contact me: alberto.roperpol@unige.ch.

## Routines

The main routines of cosmoGW are:

* [**cosmoGW.py**](cosmoGW.py): functions relevant for cosmological stochastic gravitational wave backgrounds (SGWB).
* [**cosmoMF.py**](cosmoMF.py): functions relevant for cosmological magnetic fields: bounds from different experiments, observations or projected sensitivities, and expectations from theory, among others.
* [**cosmology.py**](cosmology.py): functions relevant for cosmological calculations, including a Friedmann equations solver (see tutorial on Friedmann equations in [cosmology.ipnyb](cosmology/cosmology.ipynb)) that can generate the solution files being read in some PEncil Code simulations (see tutorial [cosmology_PC.ipnyb](cosmology/cosmology_PC.ipynb)).
* [**horndeski.py**](horndeski.py): functions relevant for GW production in the context of general theories of modified gravity.
* [**interferometry.py**](interferometry.py): functions to compute the response and sensitivity functions of interferometer space-based GW detectors (e.g., LISA and Taiji) to the detection of SGWBs (see tutorial on LISA interferometry in [interferometry.ipynb](interferometry/interferometry.ipynb)) energy density and polarization, including the space-based network LISA-Taiji to detect polarization.
* [**pta.py**](pta.py): functions used in the analysis of observations by pulsar timing array (PTA) collaborations: NANOGrav, PPTA, EPTA, and IPTA.
* [**reading.py**](reading.py): functions to read the output files of a specific set of runs (project) of the Pencil Code.
* [**run.py**](run.py): contains the class **run**, used to store all the variables computed in the Pencil Code and in cosmoGW from the Pencil Code solitions. It includes functions to initialize and postprocess the results of a set of runs.
* [**spectra.py**](spectra.py): contains description for specific spectral templates, postprocessing routines for numerical spectra, and other mathematical routines.

Some data files are available in cosmoGW that are useful in some of the projects:
* [**cosmology**](cosmology): includes files relevant for the cosmological evolution of the Universe and contains a tutorial on solving Friedmann equations.
* [**interferometry**](interferometry): includes files relevant for space-based GW interferometry calculations and contains a tutorial on computing the response functions, sensitivities and power law sensitivities to SGWB energy density and polarization.
* [**detector_sensitivity**](detector_sensitivity): includes the sensitivity of various detectors (ground-based, space-based, and pulsar timing arrays, among others), see the [README](detector_sensitivity/README.md) file for info and references.

## Projects

Each specific project is contained in a separate directory and (usually) corresponds to a publication. They all include a Jupyter notebook that allows to reproduce the results and plots of the publication. Note that (obviously) there are other authors involved in most of the collected projects so I am not the only one to credit when using these results! The Python routine [dirs.py](dirs.py) returns a dictionary linking common names of the run to their specific directories. The dictionary can be used to directly read the simulations of a specific project (see Jupyter notebooks contained in each project directory).

* The run directories and results in [**PRD_1903_08585**](PRD_1903_08585) (datasets also available in [Zenodo](https://zenodo.org/record/3692072)) correspond to:
> **A. Roper Pol, S. Mandal, A. Brandenburg, T. Kahniashvili, A. Kosowsky,** *Numerical Simulations of Gravitational Waves from Early-Universe Turbulence,* Phys. Rev. D **102**, 083512 (2020), [arXiv:1903.08585](https://arxiv.org/abs/1903.08585),
[DOI:10.1103/PhysRevD.102.083512](https://doi.org/10.1103/PhysRevD.102.083512).

* The run directories and results in [**PRR_2011_05556**](PRR_2011_05556) (datasets are also available in [Zenodo](https://zenodo.org/record/4256906)) correspond to:
> **T. Kahniashvili, A. Brandenburg, G. Gogoberidze, S. Mandal, A. Roper Pol,** *Circular polarization of gravitational waves from early-Universe helical turbulence,* Phys. Rev. Research **3**, 013193 (2021), [arXiv:2011.05556](https://arxiv.org/abs/2011.05556),
[DOI:10.1103/PhysRevResearch.3.013193](https://journals.aps.org/prresearch/abstract/10.1103/PhysRevResearch.3.013193).

* The run directories and results in [**JCAP_2107_05356**](JCAP_2107_05356) (datasets are also available in [Zenodo](https://zenodo.org/record/5525504)) correspond to:
> **A. Roper Pol, S. Mandal, A. Brandenburg, T. Kahniashvili,** *Polarization of gravitational waves from helical MHD turbulent sources,*
*JCAP* **04** (2022), 019, [arXiv:2107.05356](https://arxiv.org/abs/2107.05356).

* The run directories and results in [**PRD_2110_14456**](PRD_2110_14456) (datasets are also available in [Zenodo](https://zenodo.org/record/5603013)) correspond to:
> **Y. He, A. Roper Pol, A. Brandenburg,** *Leading-order nonlinear gravitational waves from reheating magnetogeneses,* *submitted to
Phys. Rev. D*, [arXiv:2110.14456](https://arxiv.org/abs/2110.14456) (2021).

* The run directories and results in [**PRD_2201_05630**](PRD_2201_05630) (datasets are also available in [Zenodo](https://zenodo.org/record/5782752)) correspond to:
> **A. Roper Pol, C. Caprini, A. Neronov, D. Semikoz,** *The gravitational wave signal from primordial magnetic fields in the Pulsar
Timing Array frequency band,* Phys. Rev. D **105**, 123502 (2022), [arXiv:2201.05630](https://arxiv.org/abs/2201.05630).

* The run directories and results in [**horndeski**](horndeski) (datasets are also available in
[Zenodo](https://zenodo.org/record/7408601)) correspond to:
> **Y. He, A. Roper Pol, A. Brandenburg,** *Modified propagation of gravitational waves from the early radiation era,*
*in press, JCAP* (2023), [arXiv:2212.06082](https://arxiv.org/abs/2212.06082).
