# optimal_vcmax_R
This repository contains the R functions necessary for calculating optimal vcmax in R as 
in Smith et al. Photosynthetic capacity is optimized to the environment. Ecology Letters.

The script [calc_optimal_vcmax.R](calc_optimal_vcmax.R) contains code for the optimal 
vcmax function.'
The function will calculate optimal vcmax, optimal jmax, and the optimal jmax/vcmax ratio.
The inputs required are temperature, PAR, VPD, elevation, and an estimate for the
quantum efficiency of photosynthetic electron transport, and the curvature of the light
response curve of photosynthetic electron transport.

The script [test_calc_optimal_vcmax.R](test_calc_optimal_vcmax.R) will test the functions.
It is suggested that users run this script first to ensure the function will run properly.

Any questions or issues can be submitted via GitHub or directed to Nick Smith
(nick.smith@ttu.edu).
