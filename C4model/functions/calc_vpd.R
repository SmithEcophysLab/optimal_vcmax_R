# calculate vpd in kPa, from Wang (2016) New Phytologist

calc_vpd <- function(temp, z, vpdo){
	
	vpdo_pa <- vpdo * 1000 # VPD at sea level (assumed to be the case for CRU) (Pa)
	
	po <- 101325 # Atmospheric pressure at sea level (Pa)
	
	patm <- calc_patm(z) # actual atmospheric pressure (Pa)
	
	es <- 610.8 * exp((17.27 * temp) / (237.3 + temp)) # saturation vapor pressure (Pa) http://cronklab.wikidot.com/calculation-of-vapour-pressure-deficit
	
	eao <- -(vpdo_pa - es) # actual vapor pressure at sea level (Pa)
	
	vpd <- es - eao * (patm / po) # VPD at z (Pa)
	
	vpd_kPa <- vpd / 1000 # VPD at z (kPa)
	
	return(vpd_kPa)
}




