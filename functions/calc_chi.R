# calculate the ratio of ci/ca

calc_chi = function(temp, z, vpdo, cao, beta = 146){ # temp in °C, z in m, vpd in kPa, ca in ppm
	
	patm = calc_patm(z) 
	
	vpd = calc_vpd(temp, z, vpdo) 
	
	vpd_pa = vpd * 1000
	ca_pa = cao * 1e-6 * patm
	
	km = calc_km_pa(temp, z) 
	gammastar = calc_gammastar_pa(temp, z)
	nstar = calc_nstar(temp, z)
	
	xi = sqrt((beta * (km + gammastar)) / (1.6 * nstar))
	
	chi = (gammastar / ca_pa) + (1 - (gammastar / ca_pa)) * (xi / (xi + sqrt(vpd_pa)))
	
  chi	
}

given_chi <- function(chi){
  ifelse(!is.na(chi), "yes", "no")
}

return_chi <- function(given_chi, chi, temp, z, vpdo, cao, beta){
  ifelse(given_chi == "yes",
         return(chi),
         return(calc_chi(temp = temp, z = z, vpdo = vpdo, cao = cao, beta = beta)))
}
