# calculate the ratio of ci/ca

calc_chi = function(temp, z, vpdo, cao){ # temp in °C, z in m, vpd in kPa, ca in ppm
	
	beta = 146 
	
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


