# calculate optimal vcmax as in Smith et al. Photosynthetic capacity is optimized to the 
# environment. Ecology Letters. 

########################
## variable key
########################

# tg_c: acclimated temperature (degC)
# z: elevation (m)
# vpdo: vapor pressure deficit at sea level (kPa)
# cao: atmospheric CO2 at sea level (umol mol-1)
# paro: photosynthetically active radiation at sea level (µmol m-2 s-1)
# q0: quantum efficiency of photosynthetic electron transport (mol/mol)
# theta: curvature of the light response of electron transport (unitless)
# R: universal gas constant (J mol-1 K-1)
# tg_K: acclimated temperature (K)
# to: temperature optimum for vcmax (degC)
# patm: atmospheric pressure (Pa)
# ca: atmospheric CO2 at z (Pa)
# km: Michaelis-Menten constant for Rubisco (Pa)
# gammastar: CO2 compensation point (Pa)
# chi: leaf intercellular to atmospheric CO2 ratio (ci/ca) (unitless)
# vpd: vapor pressure deficit at z (kPa)
# par: photosynthetically active radiation at z (µmol m-2 s-1)
# ci: leaf intercellular CO2 concentation (Pa)
# m: CO2 limiation of electron transport rate limited photosynthesis (Pa)
# mc: CO2 limiation of Rubisco carboxylation rate limited photosynthesis (Pa)
# c: constant describing cost of maintaining electron transport (unitless)
# omega: omega term from Smith et al.
# omega_star: omega_star term from Smith et al.
# vcmax_star: maximum rate of Rubisco carboxylation without temperature correction (µmol m-2 s-1) 
# vcmax_prime: optimal maximum rate of Rubisco carboxylation at tg (µmol m-2 s-1)
# jvrat: ratio of the maximum rate of electron transport to the maximum rate of Rubisco carbocylation at tg (unitless)
# jmax_prime: optimal maximum rate of electron transport at tg (µmol m-2 s-1)

# All equation numbers refer to equations presented in Smith et al., 2019 

# libraries
# install.packages('R.utils')
library(R.utils)

# load necessary functions
sourceDirectory('functions')

calc_optimal_vcmax <- function(tg_c = 25, z = 0, vpdo = 1, cao = 400, paro = 800, q0 = 0.257, theta = 0.85){
	
	# constants
	R <- 8.314
	c <- 0.05336251

	# environmental terms
	patm <- calc_patm(z)
	par <- calc_par(paro, z)
	vpd <- calc_vpd(tg_c, z, vpdo)
	ca <- cao * 1e-6 * patm
	to <- (0.44 * tg_c + 24.92)              # Eq. 21, note: intercept differs due to use of °C
	tg_K <- tg_c + 273.15
	
	# K and Gamma* model terms
    km <- calc_km_pa(tg_c, z)                # Eq. 3
	gammastar <- calc_gammastar_pa(tg_c, z)   

	# Coordination and least-cost hypothesis model terms
	chi <- calc_chi(tg_c, z, vpdo, cao)              # Eq. 1
	ci <- chi * ca # Pa
	mc <- ((ci - gammastar) / (ci + km))             # Eq. 6
	m <- ((ci - gammastar)/(ci + (2 * gammastar)))   # Eq. 8
	omega <- calc_omega(theta = theta, c = c, m = m) # Eq. S4
	omega_star <- (1 + (omega) - sqrt((1 + (omega))^2 - (4 * theta * omega)))  # Eq. 18
	
	# calculate vcmax and jmax	
	vcmax_star <- ((q0 * par * m) / mc) * (omega_star / (8 * theta))	          # Eq. 19
	vcmax_prime <- vcmax_star * calc_tresp_mult(tg_c, tg_c, tref = to)         # Eq. 20
	jvrat <- ((8 * theta * mc * omega) / (m * omega_star))     # Eq. 15 / Eq. 19
	jmax_prime <- jvrat * vcmax_prime
	
	# equivalently
	# jmax_prime <- q0 * par * omega * calc_tresp_mult(tg_c, tg_c, tref = to)    # Eq. 15 * Eq. 20
	# jvrat <- jmax_prime/vcmax_prime
	
	# output
	results <- as.data.frame(cbind(tg_c, z, vpdo, cao, paro, q0, theta, c, to, par, patm, ca, vpd, chi, ci, km, 
	                               gammastar, omega, m, mc, omega_star, vcmax_star, vcmax_prime, jvrat, jmax_prime))
	
	colnames(results) <- c('tg_c', 'z', 'vpdo', 'cao', 'paro', 'q0', 'theta', 'c', 'to', 'par', 'patm', 'ca', 'vpd', 'chi', 'ci', 'km', 
	                       'gammastar', 'omega', 'm', 'mc', 'omega_star', 'vcmax_star', 'vcmax_prime', 'jvrat', 'jmax_prime')
	
	results	
	
}

