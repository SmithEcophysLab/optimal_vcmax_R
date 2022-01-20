# optimal_vcmax_R

## Repository Description
This repository contains the R functions necessary for calculating optimal vcmax in R as first described 
in Smith et al. (2019) Global photosynthetic capacity is optimized to the environment.
*Ecology Letters* 22(3): 506-517. doi: 10.1111/ele.13210. [link](https://onlinelibrary.wiley.com/doi/10.1111/ele.13210).

## Summary of main files and folders
The script [calc_optimal_vcmax.R](calc_optimal_vcmax.R) contains code for the optimal 
vcmax function.'
The function will calculate optimal vcmax, optimal jmax, and the optimal jmax/vcmax ratio.
The inputs required are temperature, PAR, VPD, elevation, and an estimate for the
quantum efficiency of photosynthetic electron transport, and the curvature of the light
response curve of photosynthetic electron transport.

Th folder [functions](functions) contains the functions necessary to run the 
[calc_optimal_vcmax.R](calc_optimal_vcmax.R) script.

The script [test_calc_optimal_vcmax.R](test_calc_optimal_vcmax.R) will test the functions.
It is suggested that users run this script first to ensure the function will run properly.

All function descriptions, including parameter descriptions, can be found in the script files.

## Model Inputs
- pathway: photosynthetic pathway, either "C3" or "C4"
- tg_c: acclimated temperature (degC)
- z: elevation (m)
- vpdo: vapor pressure deficit at sea level (kPa)
- cao: atmospheric CO2 at sea level (umol mol-1)
- oao: atmospheric O2 at sea level (ppm)
- paro: photosynthetically active radiation at sea level (µmol m-2 s-1)
- q0_resp: yes or no, use the q0 response curve calculation
- q0_int: intercept for the q0 response curve calculation
- q0: quantum efficiency of photosynthetic electron transport (mol/mol)
- theta: curvature of the light response of electron transport (unitless)
- chi: leaf intercellular to atmospheric CO2 ratio (ci/ca) (unitless),  defaults to "NA"
- f: fraction of year in growing season
- lma: leaf mass area (g m-2), defaults to "NA"

## Current DOI badge
[![DOI](https://zenodo.org/badge/156727566.svg)](https://zenodo.org/badge/latestdoi/156727566)

## Contact
Any questions or issues can be submitted via GitHub or directed to Nick Smith
(nick.smith@ttu.edu).
